/*
 *  stimulation_backend_metavision.h
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

#ifndef STIMULATION_BACKEND_METAVISION_H
#define STIMULATION_BACKEND_METAVISION_H

#ifdef HAVE_METAVISION

#include "metavision/sdk/base/events/event_cd.h"
#include "metavision/sdk/driver/camera.h"
#include "nest_types.h"
#include "static_assert.h"
#include "stimulation_backend.h"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <syncstream>

/* BeginUserDocs: stimulation backend

Stimulation backend `metavision` - Receive stimulation parameters via Metavision SDK
##################################################################

The `metavision` stimulation backend collects data from event-based cameras_ using the Metavision SDK
and updates stimulation devices just before each run.

EndUserDocs */

namespace nest
{

struct CameraSpikeEvent
{
  unsigned int idx;                  // index of pixel in devices_
  long t : 31;                       // time in us
  bool positive_contrast_change : 1; // positive (1) or negative contrast change (0)

  CameraSpikeEvent() = default;
  CameraSpikeEvent( const size_t idx, const long t, const bool positive_contrast_change )
    : idx( idx )
    , t( t )
    , positive_contrast_change( positive_contrast_change )
  {
  }
};

//!< check exact bitsize
using success_event_size = StaticAssert< sizeof( CameraSpikeEvent ) == 8 >::success;


class StimulationDevice;

/**
 * A simple input backend METAVISION implementation.
 *
 * Communication protocol diagram:
 * \image html METAVISION_backend_protocol_of_communication.svg
 * General state machine diagram of NEST:
 * \image html METAVISION_backend_state_Nest.svg
 * Example of state machine diagram for the communication with NEST:
 * \image html METAVISION_backend_example_state_machine_communication_with_Nest.svg
 */
class StimulationBackendMetavision : public StimulationBackend
{
public:
  // TODO JV: Debugging only
  // std::osyncstream co;

  /**
   * InputBackend constructor
   *
   * The actual initialization is happening in InputBackend::initialize()
   */
  StimulationBackendMetavision();

  /**
   * InputBackend destructor
   *
   * The actual finalization is happening in InputBackend::finalize()
   */
  ~StimulationBackendMetavision() noexcept override;

  void initialize() override;

  void finalize() override;

  void enroll( const Node* node, StimulationDevice& device, const DictionaryDatum& params ) override;

  void disenroll( const Node* node, StimulationDevice& device ) override;

  void cleanup() override;

  void prepare() override;

  void pre_run_hook() override;

  void post_run_hook() override;

  void post_step_hook() override;

private:
  /**
   * True if there is at least one registered stimulation device.
   */
  bool enrolled_;

  /**
   * Maps (x, y) camera coordinates to stimulation devices.
   */
  std::vector< StimulationDevice* > devices_;
  /**
   * Precalculated start indices inside devices_ array for each camera index.
   */
  std::vector< size_t > camera_start_indices_;
  /**
   * List of camera resolutions (x, y).
   */
  std::vector< std::pair< size_t, size_t > > camera_resolutions_;
  /**
   * Index of the devices_ array with the position of the next pixel for which a stimulation device will be registered.
   */
  size_t next_index_;
  /**
   * List of cameras_ for which incoming spikes are received.
   */
  std::vector< Metavision::Camera > cameras_;
  /**
   * List of cameras_ for which incoming spikes are received.
   */
  std::vector< Metavision::CallbackId > camera_callbacks_;
  /**
   * Buffers spikes from the camera, with one 2D array per rank with time slots (i.e., one per update cycle) on the
   * first and the corresponding events on the second dimension.
   */
  std::vector< std::vector< std::vector< CameraSpikeEvent > > > spikes_ring_buffer_;
  /**
   * Temporary buffer for incoming spikes from the camera until the start of the next update cycle. One vector per
   * target rank.
   */
  std::vector< std::vector< CameraSpikeEvent > > temporary_spikes_buffer_;
  /**
   * Current position inside the ring buffer.
   */
  std::vector< size_t > current_time_index_;

  size_t global_num_cameras_;
  std::vector< size_t > camera_num_pixels_;
  std::vector< size_t > camera_ranks_;
  std::vector< size_t > process_devices_start_indices_;
  std::vector< size_t > process_devices_end_indices_;

  /**
   * Camera properties must not be accessed while the camera is busy (i.e., currently sending events).
   * The simulation should always have priority over accepting new spikes to prevent the camera constantly sending new
   * spikes and never freeing the mutex for long enough for the simulation thread to take over and actually process the
   * spikes. The simulation should only be blocked when all spikes have been processed and one must wait for the camera.
   */
  std::mutex camera_busy_mutex_;
  std::condition_variable camera_busy_cv_;
  /**
   * If priority_waiting_ is false, the main thread is currently busy simulating the network and the simulation backend
   * is not processing spikes from the camera. If the main thread is currently processing camera spikes
   * (priority_waiting_ == true), the camera thread should always wait for the simulation unless camera spikes for the
   * current time slice have not been added to the list of processable spikes yet (e.g. camera time is behind sim time).
   */
  std::atomic< bool > priority_waiting_;
  /**
   * Time in milliseconds when the latest camera spikes have been added to the queue. One entry per camera.
   */
  std::vector< double > camera_times_;
  /**
   * Time in milliseconds of the end of the current simulation time slice.
   */
  std::atomic< double > sim_target_time_;

  /**
   *
   */
  void
  collect_incoming_spikes_( const Metavision::EventCD* begin, const Metavision::EventCD* end, const size_t cam_idx );

  /**
   *
   */
  void change_ring_buffer_size_( const size_t rank, const size_t new_size );

  void wait_for_camera_( const size_t cam_idx, std::unique_lock< std::mutex >& lock );
};

} // namespace nest

#endif // HAVE_METAVISION

#endif /* #ifndef STIMULATION_BACKEND_METAVISION_H */
