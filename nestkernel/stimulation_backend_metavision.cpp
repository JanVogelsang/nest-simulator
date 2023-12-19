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
#include <iostream>
#include <mpi.h>
#include <vector>

// Includes from nestkernel:
#include "kernel_manager.h"
#include "nest_time.h"
#include "nest_types.h"
#include "stimulation_device.h"

#include "stimulation_backend_metavision.h"

#ifdef HAVE_METAVISION

namespace nest
{

StimulationBackendMetavision::StimulationBackendMetavision()
  : next_index_( 0 )
  , priority_waiting_( false )
  , sim_target_time_( 0 )
{
}

StimulationBackendMetavision::~StimulationBackendMetavision() noexcept
{
}


void
StimulationBackendMetavision::initialize()
{
  sim_target_time_ = 0;

  auto comm = kernel().mpi_manager.get_communicator();

  // TODO JV: Find a way to let the user specify which input source to use (e.g. cam via serial number, file, etc.)
    if ( kernel().mpi_manager.get_rank() == 0 )
    {
      cameras_.push_back( Metavision::Camera::from_file( "/home/vogelsang1/vision/camera_data/data.raw" ) );
      cameras_.back().add_runtime_error_callback( [ & ]( const Metavision::CameraException& e )
        { std::cout << "Metavision SDK reported an exception from one of the cameras: " << e.what() << std::endl; } );
//      cameras_.push_back( Metavision::Camera::from_file( "/home/vogelsang1/vision/camera_data/data.raw" ) );
//      cameras_.back().add_runtime_error_callback( [ & ]( const Metavision::CameraException& e )
//        { std::cout << "Metavision SDK reported an exception from one of the cameras: " << e.what() << std::endl; } );
    }

//  Metavision::AvailableSourcesList camera_list = Metavision::Camera::list_online_sources();
//  for ( const auto& cameras_per_type : camera_list )
//  {
//    for ( const std::string& camera_serial : cameras_per_type.second )
//    {
//      cameras_.push_back( Metavision::Camera::from_serial( camera_serial ) );
//      cameras_.back().add_runtime_error_callback( [ & ]( const Metavision::CameraException& e )
//        { std::cout << "Metavision SDK reported an exception from one of the cameras: " << e.what() << std::endl; } );
//    }
//  }

  const size_t num_processes = kernel().mpi_manager.get_num_processes();
  const int num_local_cameras = cameras_.size();

  temporary_spikes_buffer_ = std::vector< std::vector< CameraSpikeEvent > >( num_processes );
  current_time_index_ = std::vector< size_t >( num_processes, 0 );
  std::vector< size_t > process_num_devices = std::vector< size_t >( num_processes );
  // communicate information about all found cameras to other processes and receive information from them
  // first get number of cameras for each process
  std::vector< int > num_cameras_per_process = std::vector< int >( num_processes );
  MPI_Allgather( &num_local_cameras, 1, MPI_INT, num_cameras_per_process.data(), 1, MPI_INT, comm );
  global_num_cameras_ = std::accumulate( num_cameras_per_process.begin(), num_cameras_per_process.end(), 0 );
  for ( size_t rank = 0; rank < num_processes; ++rank )
  {
    for ( int idx = 0; idx != num_cameras_per_process[ rank ]; ++idx )
    {
      // set the rank each camera resides on
      camera_ranks_.push_back( rank );
    }
  }

  if ( global_num_cameras_ > 0 )
  {
    // find local start index
    const size_t local_start_index = std::accumulate(
      num_cameras_per_process.begin(), num_cameras_per_process.begin() + kernel().mpi_manager.get_rank(), 0 );
    camera_num_pixels_.resize( global_num_cameras_ );
    process_devices_start_indices_.resize( num_processes );
    process_devices_end_indices_.resize( num_processes );
    camera_resolutions_.resize( cameras_.size() );
    camera_times_.resize( cameras_.size(), 0 );

    for ( size_t cam_idx = 0; cam_idx != cameras_.size(); ++cam_idx )
    {
      camera_resolutions_[ cam_idx ] =
        std::make_pair( cameras_[ cam_idx ].geometry().width(), cameras_[ cam_idx ].geometry().height() );
      camera_num_pixels_[ local_start_index + cam_idx ] =
        camera_resolutions_[ cam_idx ].first * camera_resolutions_[ cam_idx ].second;
    }

    // Communicate the total number of pixels of each camera to determine the maximum number of devices each process
    // might have to accept for each camera. The exact value is not known yet, as it depends on the distribution of
    // devices (one per pixel) on the processes. Some processes will own one more device compared to other processes if
    // the total number of pixels is not divisible by the number of ranks.
    std::vector< int > displacements = std::vector< int >( global_num_cameras_ );
    displacements[ 0 ] = 0;
    for ( size_t cam_idx = 1; cam_idx != camera_num_pixels_.size(); ++cam_idx )
    {
      displacements[ cam_idx ] = num_cameras_per_process[ cam_idx - 1 ];
    }
    MPI_Allgatherv( camera_num_pixels_.data() + local_start_index,
      num_local_cameras,
      MPI_SIZE_T,
      camera_num_pixels_.data(),
      num_cameras_per_process.data(),
      displacements.data(),
      MPI_SIZE_T,
      comm );

    camera_start_indices_.resize( cameras_.size() );
    if ( cameras_.size() > 0 )
    {
      camera_start_indices_[ 0 ] =
        std::accumulate( camera_num_pixels_.begin(), camera_num_pixels_.begin() + local_start_index, 0 );
      for ( size_t cam_idx = 1; cam_idx != cameras_.size(); ++cam_idx )
      {
        camera_start_indices_[ cam_idx ] =
          camera_start_indices_[ cam_idx - 1 ] + camera_num_pixels_[ local_start_index + cam_idx ];
      }
    }

    size_t max_num_devices =
      ( std::accumulate( camera_num_pixels_.begin(), camera_num_pixels_.end(), 0 ) - 1 ) / num_processes + 1;
    devices_.resize( max_num_devices );

    spikes_ring_buffer_.resize( num_processes, std::vector< std::vector< CameraSpikeEvent > >( 1 ) );
  }
}

void
StimulationBackendMetavision::finalize()
{
  std::vector< StimulationDevice* >().swap( devices_ );
  next_index_ = 0;
  std::vector< Metavision::Camera >().swap( cameras_ );
  std::vector< Metavision::CallbackId >().swap( camera_callbacks_ );
  std::vector< std::vector< std::vector< CameraSpikeEvent > > >().swap( spikes_ring_buffer_ );
  std::vector< std::vector< CameraSpikeEvent > >().swap( temporary_spikes_buffer_ );
  std::vector< size_t >().swap( camera_start_indices_ );
  std::vector< std::pair< size_t, size_t > >().swap( camera_resolutions_ );
  std::vector< double >().swap( camera_times_ );
  std::vector< size_t >().swap( camera_ranks_ );
  std::vector< size_t >().swap( camera_num_pixels_ );
  std::vector< size_t >().swap( process_devices_start_indices_ );
  std::vector< size_t >().swap( process_devices_end_indices_ );
  std::vector< size_t >().swap( current_time_index_ );
}

void
StimulationBackendMetavision::enroll( const Node* node, StimulationDevice& device, const DictionaryDatum& params )
{
  if ( next_index_ >= devices_.size() )
  {
    throw KernelException( "Metavision stimulation backend is not accepting any more stimulation devices!" );
  }
  // Check if device has been enrolled already to prevent enrolling the same device multiple times.
  // auto device_it = std::find( devices_.begin(), devices_.begin() + next_index_, &device );
  // if ( device_it == devices_.begin() + next_index_ )
  // {
  devices_[ next_index_ ] = &device;
  // Find next pixel for which no stimulation device has been registered yet, starting from the last set index.
  next_index_ = std::find( devices_.begin() + next_index_ + 1, devices_.end(), nullptr ) - devices_.begin();
  // }
}

void
StimulationBackendMetavision::disenroll( const Node* node, StimulationDevice& device )
{
  auto device_it = std::find( devices_.begin(), devices_.begin() + next_index_, &device );
  if ( device_it != devices_.end() )
  {
    ( *device_it ) = nullptr;
    next_index_ = device_it - devices_.begin();
  }
}

void
StimulationBackendMetavision::prepare()
{
}

void
StimulationBackendMetavision::wait_for_camera_( const size_t cam_idx, std::unique_lock< std::mutex >& lock )
{
  // If the camera didn't reach simulation time yet, we continue waiting. Make sure this is only checked while the mutex
  // is locked by us, otherwise race conditions can occur when accessing the camera_times_.
  while ( camera_times_[ cam_idx ] < sim_target_time_ )
  {
    // wait for notification of camera that more spikes have been sent and the camera time has been advanced
    camera_busy_cv_.wait( lock );
  }
}

void
StimulationBackendMetavision::pre_run_hook()
{
  // Get the total number of devices on each process and communicate to everyone, so that all processes managing cameras
  // can send the events for all pixels depending on the number of devices. Device-pixel mapping is performed round-
  // robin, but some processes could in theory own more devices than others, depending on how the devices were created.
  // In the average case, when a whole population is created per camera, the round-robin mapping will match the round-
  // robin distribution of the devices on the processes. I.e., the neuron creation order (represented by global node
  // ids) matches the pixel order (row-wise).
  // communicate the number of devices per process
  auto comm = kernel().mpi_manager.get_communicator();
  MPI_Allgather( &next_index_, 1, MPI_SIZE_T, process_devices_start_indices_.data(), 1, MPI_SIZE_T, comm );
  // construct the global starting index for the first device of each rank based on the number of devices on each rank

  process_devices_end_indices_[ 0 ] = process_devices_start_indices_[ 0 ];
  for ( size_t idx = 1; idx != process_devices_start_indices_.size(); ++idx )
  {
    process_devices_start_indices_[ idx ] += process_devices_start_indices_[ idx - 1 ];
    process_devices_end_indices_[ idx ] = process_devices_start_indices_[ idx ];
  }
  // shift array by one index to the right
  std::rotate( process_devices_start_indices_.begin(),
    process_devices_start_indices_.end() - 1,
    process_devices_start_indices_.end() );
  process_devices_start_indices_[ 0 ] = 0;

  // start listening for incoming camera spikes as soon as all buffers have been set up and the simulation starts
  for ( size_t cam_idx = 0; cam_idx != cameras_.size(); ++cam_idx )
  {
    camera_callbacks_.push_back( cameras_[ cam_idx ].cd().add_callback(
      [ =, this ]( const Metavision::EventCD* begin, const Metavision::EventCD* end )
      { collect_incoming_spikes_( begin, end, cam_idx ); } ) );
    cameras_[ cam_idx ].start();
  }

  // manually run post_step_hook once, to initially load spikes from the camera
  post_step_hook();
}

void
StimulationBackendMetavision::post_step_hook()
{
  auto comm = kernel().mpi_manager.get_communicator();
  // First wait for the camera to finish producing and outputting spikes for the current time slice, and then iterate
  // over all currently produced spikes to insert them into the ring buffer. After processing all temporary spikes the
  // number of spikes produced in the current slice are communicated to all other processes to then exchange all spikes
  // in a subsequent step.
  const long current_timestep = kernel().simulation_manager.get_slice_origin().get_steps();
  const double current_timestep_double = Time::delay_steps_to_ms( current_timestep );
  // get end time of next time slice
  sim_target_time_ = Time::delay_steps_to_ms( current_timestep + kernel().simulation_manager.get_to_step() );

  const double time_slice_length = Time::delay_steps_to_ms( kernel().connection_manager.get_min_delay() );
  const size_t num_processes = kernel().mpi_manager.get_num_processes();

  std::unique_lock< std::mutex > lock;

#pragma omp master
  {
    // prevent the camera from taking over the mutex again
    priority_waiting_ = true;

    lock = std::unique_lock( camera_busy_mutex_ );  // implicitly locks the mutex

    // block until all cameras_ reached the current timestep
    for ( size_t cam_idx = 0; cam_idx != cameras_.size(); ++cam_idx )
    {
      // only wait for new camera spikes if the camera is still turned on
      if ( cameras_[ cam_idx ].is_running() )
      {
        wait_for_camera_( cam_idx, lock );
      }
    }
    priority_waiting_ = false;
  }

// Parallel sorting of temporarily buffered spikes into time slices
#pragma omp for
  for ( size_t rank = 0; rank != num_processes; ++rank )
  {
    for ( const CameraSpikeEvent& cse : temporary_spikes_buffer_[ rank ] )
    {
      // convert time from microsecond (long) to milliseconds (double)
      double t = static_cast< double >( cse.t ) / 1000;

      // calculate how many simulation cycles in the future the event has to be sent
      const size_t cycles_in_future = ( t - current_timestep_double ) / time_slice_length;

      // Check if the spikes must be inserted into a position in the ring buffer that does not exist yet. In this
      // case the ring buffer must be resized first.
      if ( cycles_in_future >= spikes_ring_buffer_[ rank ].size() )
      {
        change_ring_buffer_size_( rank, cycles_in_future + 1 );
      }
      spikes_ring_buffer_[ rank ]
                         [ ( current_time_index_[ rank ] + cycles_in_future ) % spikes_ring_buffer_[ rank ].size() ]
                           .push_back( cse );
    }
    temporary_spikes_buffer_[ rank ].clear();
  } // implicit barrier
#pragma omp master
  {
    lock.unlock();
    // Notify camera thread that spikes can be collected again now, as the temporary_spikes_buffer_ is no longer in use.
    camera_busy_cv_.notify_all();

    // All spikes need to be exchanged. In order to do so, we must first communicate the number of spikes for each rank,
    // to then calculate the receive buffer sizes on each rank and allocate the memory.
    // first communicate the number of spikes to expect
    std::vector< int > send_sizes = std::vector< int >( num_processes );
    std::vector< int > recv_sizes = std::vector< int >( num_processes );
    for ( size_t rank = 0; rank != num_processes; ++rank )
    {
      send_sizes[ rank ] = spikes_ring_buffer_[ rank ][ current_time_index_[ rank ] ].size();
    }

    MPI_Alltoall( send_sizes.data(), 1, MPI_INT, recv_sizes.data(), 1, MPI_INT, comm );
    // now allocate memory
    const size_t send_total_num = std::accumulate( send_sizes.begin(), send_sizes.end(), 0 );
    const size_t recv_total_num = std::accumulate( recv_sizes.begin(), recv_sizes.end(), 0 );
    std::vector< CameraSpikeEvent > send_buffer = std::vector< CameraSpikeEvent >( send_total_num );
    std::vector< CameraSpikeEvent > recv_buffer = std::vector< CameraSpikeEvent >( recv_total_num );
    std::vector< int > send_displacements = std::vector< int >( num_processes );
    std::vector< int > recv_displacements = std::vector< int >( num_processes );
    // collocate all spikes in the send buffer
    // This solution (based on implementation in mpi_manager) is very dirty and not portable, but efficient.
    const size_t effective_size_multiplier = sizeof( CameraSpikeEvent ) / sizeof( size_t );
    size_t current_idx = 0;
    for ( size_t rank = 0; rank != num_processes; ++rank )
    {
      send_sizes[ rank ] *= effective_size_multiplier;
      recv_sizes[ rank ] *= effective_size_multiplier;

      send_displacements[ rank ] = current_idx;
      assert( send_buffer.size() - current_idx >= spikes_ring_buffer_.at( rank ).at( current_time_index_.at( rank ) ).size() );
      std::copy( spikes_ring_buffer_.at( rank ).at( current_time_index_.at( rank ) ).begin(),
        spikes_ring_buffer_.at( rank ).at( current_time_index_.at( rank ) ).end(),
        send_buffer.begin() + current_idx );
      current_idx += spikes_ring_buffer_.at( rank ).at( current_time_index_.at( rank ) ).size();
    }
    recv_displacements[ 0 ] = 0;
    for ( size_t rank = 1; rank != num_processes; ++rank )
    {
      recv_displacements[ rank ] = recv_displacements[ rank - 1 ] + recv_sizes[ rank - 1 ];
    }
    // now communicate the actual spike data
    MPI_Alltoallv( send_buffer.data(),
      send_sizes.data(),
      send_displacements.data(),
      MPI_SIZE_T,
      recv_buffer.data(),
      recv_sizes.data(),
      recv_displacements.data(),
      MPI_SIZE_T,
      comm );
    for ( CameraSpikeEvent cse : recv_buffer )
    {
      // convert time from microsecond (long) to milliseconds (double)
      double t = static_cast< double >( cse.t ) / 1000;

      // send spike to device for pixel
      std::vector< double > spike { t, cse.positive_contrast_change ? 1. : -1. };

      devices_.at( cse.idx )->set_data_from_stimulation_backend( spike );
    }
    for ( size_t rank = 0; rank != num_processes; ++rank )
    {
      spikes_ring_buffer_[ rank ][ current_time_index_[ rank ] ].clear();
      current_time_index_[ rank ] = ( current_time_index_[ rank ] + 1 ) % spikes_ring_buffer_[ rank ].size();
    }
  }
}

void
StimulationBackendMetavision::post_run_hook()
{
  // stop listening for incoming camera spikes when the simulation ends
  for ( size_t cam_idx = 0; cam_idx != cameras_.size(); ++cam_idx )
  {
    cameras_[ cam_idx ].stop();
    cameras_[ cam_idx ].cd().remove_callback( camera_callbacks_[ cam_idx ] );
  }
  camera_callbacks_.clear();
}

void
StimulationBackendMetavision::cleanup()
{
}

void
StimulationBackendMetavision::collect_incoming_spikes_( const Metavision::EventCD* begin,
  const Metavision::EventCD* end,
  const size_t cam_idx )
{
  {
    std::unique_lock< std::mutex > l( camera_busy_mutex_ );
    camera_busy_cv_.wait( l, [ & ]() { return !priority_waiting_ or sim_target_time_ > camera_times_[ cam_idx ]; } );

    // we can discard all spikes that will never be processes by the simulation anymore
    const long simulation_end_time_us =
      static_cast< long >( kernel().simulation_manager.run_duration().get_ms() * 1000 );
    // copy all events into temporary buffer
    for ( const Metavision::EventCD* ev = begin; ev != end; ++ev )
    {
      const long time_us = ev->t;
      if ( time_us > simulation_end_time_us )
      {
        break;
      }

      const bool positive_contrast_change = ev->p == 1;
      const unsigned int pixel_idx =
        camera_start_indices_[ cam_idx ] + camera_resolutions_[ cam_idx ].first * ev->y + ev->x;

      // add +1 to pixel_idx to not find the first value "not less than", but the first "strictly greater than"
      const auto target_rank_it =
        std::lower_bound( process_devices_start_indices_.begin(), process_devices_start_indices_.end(), pixel_idx + 1 );
      const size_t rank = ( target_rank_it - process_devices_start_indices_.begin() ) - 1;
      if ( pixel_idx < process_devices_end_indices_[ rank ] )
      {
        const size_t rank_local_pixel_idx = pixel_idx - *( target_rank_it - 1 );
        temporary_spikes_buffer_[ rank ].push_back(
          CameraSpikeEvent( rank_local_pixel_idx, time_us, positive_contrast_change ) );
      }
    }
    if ( begin != end )
    {
      camera_times_[ cam_idx ] = static_cast< double >( ( end - 1 )->t ) / 1000;
    }
  }
  camera_busy_cv_.notify_all();
}

void
StimulationBackendMetavision::change_ring_buffer_size_( const size_t rank, const size_t new_size )
{
  const size_t old_size = spikes_ring_buffer_[ rank ].size();
  spikes_ring_buffer_[ rank ].resize( new_size );

  if ( new_size > old_size )
  {
    // The ring buffer must be extended at the end, however the end is not necessarily at the end of the vector.
    for ( size_t i = 0; i != current_time_index_[ rank ]; ++i )
    {
      spikes_ring_buffer_[ rank ][ ( old_size + i ) % spikes_ring_buffer_[ rank ].size() ].swap(
        spikes_ring_buffer_[ rank ][ i ] );
    }
  }
}

} // namespace nest

#endif // HAVE_METAVISION
