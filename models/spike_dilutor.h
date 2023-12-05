/*
 *  spike_dilutor.h
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

#ifndef SPIKE_DILUTOR_H
#define SPIKE_DILUTOR_H

// Includes from nestkernel:
#include "connection.h"
#include "device_node.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "stimulation_device.h"

namespace nest
{

/* BeginUserDocs: device, generator

Short description
+++++++++++++++++

Repeat incoming spikes with a certain probability

Description
+++++++++++

The device repeats incoming spikes with a certain probability.
Targets will receive different spike trains.

In parallel simulations, a copy of the device is present on each process
and spikes are collected only from local sources.

.. admonition:: Deprecated model

   ``spike_dilutor`` is deprecated because it does not work with multiple threads.
   To create connections that transmit spikes with a given probability, use :doc:`bernoulli_synapse <bernoulli_synapse>`
   instead.

.. admonition:: Does not work with threads

   ``spike_dilutor`` only works in single-threaded simulations. It can be used with MPI-parallel simulations.

Parameters
++++++++++

p_copy
    Copy probability

Sends
+++++

SpikeEvent

See also
++++++++

mip_generator


Examples using this model
+++++++++++++++++++++++++

.. listexamples:: spike_dilutor

EndUserDocs */

void register_spike_dilutor( const std::string& name );

class spike_dilutor : public DeviceNode
{

public:
  spike_dilutor();
  spike_dilutor( const spike_dilutor& rhs );

  bool
  has_proxies() const override
  {
    return false;
  }

  bool
  local_receiver() const override
  {
    return true;
  }

  Name
  get_element_type() const override
  {
    return names::stimulator;
  }

  using Node::event_hook;
  using Node::handle;
  using Node::handles_test_event; // new

  size_t send_test_event( Node&, size_t, synindex, bool ) override;
  size_t handles_test_event( SpikeEvent&, size_t ) override;
  void handle( SpikeEvent& ) override;

  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum& ) override;

  void
  set_initialized_() final
  {
    device_.set_initialized_( this );
  }

private:
  void init_state_() override;
  void init_buffers_() override;
  void pre_run_hook() override;

  void update( Time const&, const long, const long ) override;

  void event_hook( DSSpikeEvent& ) override;

  // ------------------------------------------------------------

  /**
   * Store independent parameters of the model.
   */
  struct Parameters_
  {
    double p_copy_; //!< copy probability for each incoming spike

    Parameters_(); //!< Sets default parameter values
    Parameters_( const Parameters_& ) = default;
    Parameters_& operator=( const Parameters_& ) = default;

    void get( DictionaryDatum& ) const;             //!< Store current values in dictionary
    void set( const DictionaryDatum&, Node* node ); //!< Set values from dictionary
  };

  struct Buffers_
  {
    RingBuffer n_spikes_;
  };

  // ------------------------------------------------------------

  class DilutorStimulationDevice : public StimulationDevice
  {
    StimulationDevice::Type
    get_type() const override
    {
      return StimulationDevice::Type::SPIKE_GENERATOR;
    }

  public:
    using StimulationDevice::set_initialized_;
    void
    set_initialized_( const Node* node )
    {
      StimulationDevice::set_initialized_( node );
    }
  } device_;

  Parameters_ P_;
  Buffers_ B_;
};

inline size_t
spike_dilutor::send_test_event( Node& target, size_t receptor_type, synindex syn_id, bool )
{
  device_.enforce_single_syn_type( syn_id );

  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline size_t
spike_dilutor::handles_test_event( SpikeEvent&, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline void
spike_dilutor::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  device_.get_status( d );
}

inline void
spike_dilutor::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d, this );   // throws if BadProperty

  // We now know that ptmp is consistent. We do not write it back
  // to P_ before we are also sure that the properties to be set
  // in the parent class are internally consistent.
  device_.set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
}

} // namespace nest

#endif
