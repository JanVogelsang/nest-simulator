/*
 *  metavision_spike_train_injector.h
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

#ifndef PRECISE_WEIGHTED_SPIKE_GENERATOR_H
#define PRECISE_WEIGHTED_SPIKE_GENERATOR_H


// C++ includes:
#include <vector>

// Includes from nestkernel:
#include "connection.h"
#include "device_node.h"
#include "event.h"
#include "nest_time.h"
#include "nest_types.h"
#include "stimulation_device.h"

namespace nest
{

/* BeginUserDocs: neuron, device, spike, generator

Short description
+++++++++++++++++

Generate spikes from an array with precise spike-times and weights

Description
+++++++++++

A spike generator can be used to generate weighted spikes at specific times which are given to the spike generator as an
array. Very similar to the ``spike_generator``, but always uses precise spikes and adds a weight to each spike.
Stimulation backends using precise spikes can be used with this generator.

Spike times are given in milliseconds as an array. The `spike_times` array must be sorted with the earliest spike first.
 All spike times must be strictly in the future. Trying to set a spike time in the past or at the current time step will
 cause a NEST error. Setting a spike time of 0.0 will also result in an error.

Multiple occurrences of the same time indicate that more than one event is to be generated at this particular time.

Additionally, `spike_weights` can be set. This is an array as well. It contains one weight value per spike time. The
 spikes are delivered with the respective weight multiplied with the weight of the connection.

The spike generator supports spike times that do not coincide with a time step, that is, are not falling on the grid
 defined by the simulation resolution.
Spike times will not be rounded but represented exactly as a combination of step and offset.


.. include:: ../models/stimulation_device.rst

spike_times
    List of spike times in ms.

spike_weights
    List of corresponding spike weights, the unit depends on the receiver.
    (e.g., nS for conductance-based neurons or pA for current based ones)

spike_offsets
    List of spike time offsets.


Set spike times from a stimulation backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spike times for this stimulation device can be updated with input coming from a stimulation backend.

Sends
+++++

SpikeEvent

EndUserDocs
*/
void register_precise_weighted_spike_generator( const std::string& name );


class precise_weighted_spike_generator : public Node, public StimulationDevice
{

public:
  precise_weighted_spike_generator();
  precise_weighted_spike_generator( const precise_weighted_spike_generator& );

  size_t send_test_event( Node&, size_t, synindex, bool ) override;
  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum& ) override;
  bool
  is_active( const Time& ) const override
  {
    return true;
  }
  bool
  is_off_grid() const override
  {
    return true;
  }

  StimulationDevice::Type get_type() const override;
  void set_data_from_stimulation_backend( std::vector< double >& input_spikes ) override;

private:
  void init_state_() override;
  void init_buffers_() override;
  void pre_run_hook() override;

  void update( Time const&, const long, const long ) override;

  void
  set_initialized_() final
  {
    StimulationDevice::set_initialized_( this );
  }
  Name
  get_element_type() const override
  {
    return names::stimulator;
  }

  // ------------------------------------------------------------

  struct State_
  {
    State_();
    size_t position_; //!< index of next spike to deliver
  };

  // ------------------------------------------------------------

  struct Parameters_
  {
    //! Spike time stamp as Time, rel to origin_
    std::vector< Time > spike_stamps_;

    //! Spike time offset
    std::vector< double > spike_offsets_;

    std::vector< double > spike_weights_; //!< Spike weights

    Parameters_(); //!< Sets default parameter values
    Parameters_( const Parameters_& ) = default;
    Parameters_& operator=( const Parameters_& ) = default;

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary

    /**
     * Set values from dictionary.
     * @note State is passed so that the position can be reset if the
     *       spike_times_ or spike_weights_ vector has been filled with
     *       new data, or if the origin was reset.
     */
    void set( const DictionaryDatum&, State_&, const Time&, const Time&, Node* node );

    /**
     * Insert spike time to arrays, throw BadProperty for invalid spike times.
     *
     * @param spike time, ms
     * @param origin
     * @param current simulation time
     */
    void assert_valid_spike_time_and_insert_( double, const Time&, const Time& );
  };

  // ------------------------------------------------------------

  Parameters_ P_;
  State_ S_;
};

inline size_t
precise_weighted_spike_generator::send_test_event( Node& target, size_t receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline void
precise_weighted_spike_generator::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  StimulationDevice::get_status( this, d );
}

inline void
nest::precise_weighted_spike_generator::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors

  // To detect "now" spikes and shift them, we need the origin. In case
  // it is set in this call, we need to extract it explicitly here.
  Time origin;
  double v;
  if ( updateValue< double >( d, names::origin, v ) )
  {
    origin = Time::ms( v );
  }
  else
  {
    origin = StimulationDevice::get_origin();
  }

  // throws if BadProperty
  ptmp.set( d, S_, origin, kernel().simulation_manager.get_time(), this );

  // We now know that ptmp is consistent. We do not write it back
  // to P_ before we are also sure that the properties to be set
  // in the parent class are internally consistent.
  StimulationDevice::set_status( this, d );

  // if we get here, temporary contains consistent set of properties
  P_ = ptmp;
}

inline StimulationDevice::Type
precise_weighted_spike_generator::get_type() const
{
  return StimulationDevice::Type::SPIKE_GENERATOR;
}

} // namespace nest

#endif /* #ifndef PRECISE_WEIGHTED_SPIKE_GENERATOR_H */
