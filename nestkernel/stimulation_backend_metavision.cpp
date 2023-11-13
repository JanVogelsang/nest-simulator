/*
 *  stimulation_backend_metavision.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// C++ includes:
#include <fstream>
#include <iostream>
#include <vector>

// Includes from nestkernel:
#include "kernel_manager.h"
#include "nest_time.h"
#include "nest_types.h"
#include "stimulation_backend.h"
#include "stimulation_device.h"

#include "stimulation_backend_metavision.h"

#ifdef HAVE_METAVISION

namespace nest
{

StimulationBackendMetavision::StimulationBackendMetavision()
  : enrolled_( false )
  , pcc_filled_( false )
{
}

StimulationBackendMetavision::~StimulationBackendMetavision() noexcept
{
}


void
StimulationBackendMetavision::initialize()
{
  Metavision::AvailableSourcesList camera_list = Metavision::Camera::list_online_sources();
  size_t cam_idx = 0;
  for ( const auto &cameras_per_type : camera_list )
  {
      for ( const std::string& camera_serial : cameras_per_type.second )
      {
          cameras.push_back( Metavision::Camera::from_serial( camera_serial ) );
          // cameras.push_back( Metavision::Camera::from_file( "/home/vogelsang1/CLionProjects/metavision-sample/data/data.raw" ) );
          cameras[ cam_idx ].cd().add_callback( [=]( const Metavision::EventCD *begin, const Metavision::EventCD *end ) { collect_incoming_spikes( begin, end, cam_idx ); });
          cameras[ cam_idx ].add_runtime_error_callback( []( const Metavision::CameraException &e){ std::cout << "Metavision SDK reported an exception from one of the cameras: " << e.what() << std::endl; } );
          camera_resolutions.emplace_back( cameras[ cam_idx ].geometry().width(), cameras[ cam_idx ].geometry().height() );
          cameras[ cam_idx ].start();
          ++cam_idx;
      }
  }

  camera_start_indices.push_back( 0 );
  for ( cam_idx = 0; cam_idx != cameras.size() - 1; ++cam_idx )
  {
      camera_start_indices.push_back( camera_resolutions[ cam_idx - 1 ].first * camera_resolutions[ cam_idx - 1 ].second + camera_start_indices[ cam_idx - 1 ] );
  }

  // two stimulation devices per camera pixel
  devices_pcc_ = std::vector< StimulationDevice* >( camera_start_indices.back() + camera_resolutions.back().first * camera_resolutions.back().second );
  devices_ncc_ = std::vector< StimulationDevice* >( camera_start_indices.back() + camera_resolutions.back().first * camera_resolutions.back().second );
}

void
StimulationBackendMetavision::finalize()
{
  std::vector< StimulationDevice* >().swap( devices_pcc_ );
  std::vector< StimulationDevice* >().swap( devices_ncc_ );

  for ( Metavision::Camera& cam : cameras )
  {
      cam.stop();
  }
}

void
StimulationBackendMetavision::enroll( StimulationDevice& device, const DictionaryDatum& params )
{
  std::vector< StimulationDevice* >& devices_ = pcc_filled_ ? devices_ncc_ : devices_pcc_;

  devices_[ next_index ] = &device;
  // Find next pixel for which no stimulation device has been registered yet, starting from the last set index.
  next_index = std::find( devices_.begin() + next_index + 1, devices_.end(), nullptr ) - devices_.begin();

  enrolled_ = true;
}


void
StimulationBackendMetavision::disenroll( StimulationDevice& device )
{
  auto device_it = std::find( devices_pcc_.begin(), devices_pcc_.end(), &device );
  if ( device_it != devices_pcc_.end() )
  {
    (*device_it) = nullptr;
    pcc_filled_ = false;
  } else
  {
    device_it = std::find( devices_ncc_.begin(), devices_ncc_.end(), &device );
    if ( device_it != devices_ncc_.end() )
    {
      (*device_it) = nullptr;
    }
  }
}

void
StimulationBackendMetavision::prepare()
{
  if ( not enrolled_ )
  {
    return;
  }
}

void
StimulationBackendMetavision::pre_run_hook()
{
#pragma omp single
  {
    // block until all cameras reached the current timestep
    for ( Metavision::Camera& cam : cameras )
    {
      const double current_timestep = kernel().simulation_manager.get_slice_origin().get_ms();
      while ( cam.is_running()
        and static_cast< double >( cam.get_last_timestamp() ) / 1000 < current_timestep )
      {
        std::this_thread::sleep_for( std::chrono::milliseconds( 1 ) );  // TODO JV: Remove sleep
      }
    }
  }

  // TODO JV: Replace omp single by parallel for (for each camera)
#pragma omp single
  {
    // Receive spike trains from all cameras until the current time has been reached
    for ( CameraSpikeEvent cse : spikes_ring_buffer[ current_time_index ] )
    {
        std::cout << "Event received: coordinates (" << cse.idx << "), t: " << cse.t << ", polarity: " << cse.positive_contrast_change << std::endl;
        // send spike to device for pixel
        if ( cse.positive_contrast_change )
        {
            // TODO JV: idx might be larger than size if not enough neurons registered
            devices_pcc_[ cse.idx ]->set_data_from_stimulation_backend( cse.t, cse.offset );
        }
        else {
            devices_ncc_[ cse.idx ]->set_data_from_stimulation_backend( cse.t, cse.offset );
        }
    }
    spikes_ring_buffer[ current_time_index ].clear();
    current_time_index = ( current_time_index + 1 ) % spikes_ring_buffer.size();
  }
}

void
StimulationBackendMetavision::post_run_hook()
{
}

void
StimulationBackendMetavision::cleanup()
{
}

void StimulationBackendMetavision::collect_incoming_spikes( const Metavision::EventCD *begin, const Metavision::EventCD *end, const size_t cam_idx )
{
  const double current_timestep = kernel().simulation_manager.get_slice_origin().get_ms();
  for ( const Metavision::EventCD *ev = begin; ev != end; ++ev )
  {
    double precise_time_us = static_cast< double >( ev->t );
    unsigned int t = Time::delay_ms_to_steps( precise_time_us / 1000 );
    double offset = ( precise_time_us - t * static_cast< double >( Time::resolution ) * 1000 ) / 1000;  // TODO JV: verify correctness
    bool positive_contrast_change = ev->p == 1;
    unsigned int pixel_idx = camera_start_indices[ cam_idx ] + camera_resolutions[ cam_idx ].first * ev->y + ev->x;

    // calculate how many simulation cycles in the future the event has to be sent
    const size_t cycles_in_future = ( t - origin ) / ( origin_to - origin );
    // Check if the spikes must be inserted into a position in the ring buffer that does not exist yet. In this case
    // the ring buffer must be resized first.
    if ( cycles_in_future >= spikes_ring_buffer.size() )
    {
        change_ring_buffer_size( cycles_in_future + 1 );
    }
    spikes_ring_buffer[ ( current_time_index + cycles_in_future ) % spikes_ring_buffer.size() ].push_back( CameraSpikeEvent{ pixel_idx, t, positive_contrast_change, offset } );
  }
}

void StimulationBackendMetavision::change_ring_buffer_size( const size_t new_size )
{
    size_t size_diff = new_size - spikes_ring_buffer.size();  // Only required (and valid) when size increases
    spikes_ring_buffer.resize( new_size );
    if ( new_size > spikes_ring_buffer.size() )
    {
        // The ring buffer must be extended at the end, however the end is not necessarily at the end of the vector.
        for ( size_t i = 0; i != std::max( size_diff, current_time_index ); ++i )
        {
            spikes_ring_buffer[ current_time_index + i + 1 ].swap( spikes_ring_buffer[ i ] );
        }
    }
}

} // namespace nest

#endif  // HAVE_METAVISION
