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

#include "nest_types.h"
#include "stimulation_backend.h"
#include "static_assert.h"
#include "metavision/sdk/base/events/event_cd.h"
#include "metavision/sdk/driver/camera.h"

/* BeginUserDocs: stimulation backend

Stimulation backend `metavision` - Receive stimulation parameters via Metavision SDK
##################################################################

The `metavision` stimulation backend collects data from event-based cameras using the Metavision SDK
and updates stimulation devices just before each run.

EndUserDocs */

namespace nest
{

struct CameraSpikeEvent {
    unsigned int idx;               // index of pixel in devices_
    long t : NUM_BITS_LAG;          // time in steps
    bool positive_contrast_change : 1;  // positive (1) or negative contrast change (0)
    unsigned short offset;          // time offset relative to t in us
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

  void enroll( StimulationDevice& device, const DictionaryDatum& params ) override;

  void disenroll( StimulationDevice& device ) override;

  void cleanup() override;

  void prepare() override;

  void pre_run_hook() override;

  void post_run_hook() override;

private:
  /**
   * True if there is at least one registered stimulation device.
   */
  bool enrolled_;

  /**
   * Maps (x, y) camera coordinates to stimulation devices for positive and negative contrast changes (pcc and ncc).
   */
  std::vector< StimulationDevice* > devices_pcc_;
  std::vector< StimulationDevice* > devices_ncc_;
  bool pcc_filled_;
  /**
   * Precalculated start indices inside devices_ array for each camera index.
   */
  std::vector< size_t > camera_start_indices;
  /**
   * List of camera resolutions (x, y).
   */
  std::vector< std::pair< size_t, size_t > > camera_resolutions;
  /**
   * Index of the devices_ array with the position of the next pixel for which a stimulation device will be registered.
   */
  size_t next_index;
  /**
   * List of cameras for which incoming spikes are received.
   */
  std::vector< Metavision::Camera > cameras;
  /**
   * Buffer of incoming spikes from the camera, with one vector per time slot (i.e., one per update cycle).
   */
  std::vector< std::vector < CameraSpikeEvent > > spikes_ring_buffer;
  /**
   * Current position inside the ring buffer.
   */
  size_t current_time_index = 0;

  /**
   *
   */
  void collect_incoming_spikes( const Metavision::EventCD *begin, const Metavision::EventCD *end, const size_t cam_idx );

  /**
   *
   */
  void change_ring_buffer_size( const size_t new_size );
};

} // namespace nest

#endif  // HAVE_METAVISION

#endif /* #ifndef STIMULATION_BACKEND_METAVISION_H */
