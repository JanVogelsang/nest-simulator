/*
 *  connection.h
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

#ifndef CONNECTION_H
#define CONNECTION_H

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection_label.h"
#include "connector_base_impl.h"
#include "delay_checker.h"
#include "event.h"
#include "kernel_manager.h"
#include "nest_names.h"
#include "nest_time.h"
#include "nest_timeconverter.h"
#include "nest_types.h"
#include "node.h"
#include "spikecounter.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"

namespace nest
{

class ConnectorModel;

/**
 * Base class for dummy nodes used in connection testing.
 *
 * This class provides a based for dummy node objects that
 * are used to test whether a connection can be established.
 * The base class provides empty implementations of all pure
 * virtual functions of class Node.
 *
 * Each connection class (i.e., each class derived from class
 * template Connection<T>), must derive a concrete ConnTestDummyNode
 * class that overrides method Node::handles_test_event() for all
 * event types that the connection supports.
 *
 * For details, see Kunkel et al, Front Neuroinform 8:78 (2014),
 * Sec 3.3.1. Note that the ConnTestDummyNode class is called
 * "check_helper" in the paper.
 *
 * @ingroup event_interface
 */
class ConnTestDummyNodeBase : public Node
{
  void
  pre_run_hook() override
  {
  }
  void
  update( const nest::Time&, const long, const long ) override
  {
  }
  void
  set_status( const DictionaryDatum& ) override
  {
  }
  void
  get_status( DictionaryDatum& ) const override
  {
  }
  void
  init_state_() override
  {
  }
  void
  init_buffers_() override
  {
  }
};

/**
 * Base class for representing connections.
 * It provides the mandatory properties receiver port and target,
 * as well as the functions get_status() and set_status()
 * to read and write them. A suitable Connector containing these
 * connections can be obtained from the template GenericConnector.
 *
 * \note Please note that the event received by the send() function is
 * a reference to a single object that is re-used by each Connection.
 * This means that the object must not be changed in the Connection,
 * or if needs to be changed, everything has to be reset after sending
 * (i.e. after Event::operator() has been called).
 */
class Connection
{
public:
  // this typedef may be overwritten in the derived connection classes in order
  // to attach a specific event type to this connection type, used in secondary
  // connections not used in primary connectors
  typedef SecondaryEvent EventType;

  Connection()
  {
    set_delay( 1.0 );
  }

  Connection( const Connection& rhs ) = default;
  Connection& operator=( const Connection& rhs ) = default;

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  };

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   *
   * @note Target and Rport cannot be changed after a connection has been
   * created.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Check syn_spec dictionary for parameters that are not allowed with the
   * given connection.
   *
   * Will issue warning or throw error if an illegal parameter is found. The
   * method does nothing if no illegal parameter is found.
   *
   * @note Classes requiring checks need to override the function with their own
   * implementation, as this base class implementation does not do anything.
   */
  void check_synapse_params( const DictionaryDatum& d ) const;

  /**
   * Calibrate the delay of this connection to the desired resolution.
   */
  void calibrate( const TimeConverter& );

  /**
   * Framework for STDP with predominantly axonal delays:
   * Correct this synapse and the corresponding previously sent spike
   * taking into account a new post-synaptic spike.
   */
  void correct_synapse_stdp_ax_delay( const double t_last_pre_spike,
    double* weight_revert,
    const double t_post_spike,
    const synindex syn_id,
    const CommonSynapseProperties&,
    Node* target );

  /**
   * Process a post-synaptic spike after it is backpropagated to the synapse.
   * @param t_syn The time the post-synaptic spike arrives at this connection
   */
  void process_post_synaptic_spike( const double t_syn, const CommonSynapseProperties& );

  /**
   * Get the total transmission delay.
   */
  double
  get_delay() const
  {
    // TODO JV (pt): Think about storing dendritic delay and assuming an axonal delay of 0. In derived connections
    //  (AxonalDelayConnections), both the dendritic and axonal delay would be stored and the total delay would be the
    //  sum instead of the other way around as it is the case right now.
    return Time::delay_steps_to_ms( delay_ );
  }

  /**
   * Get the proportion of the transmission delay attributed to the dendrite.
   */
  double
  get_dendritic_delay() const
  {
    return get_delay();
  }

  /**
   * Get the proportion of the transmission delay attributed to the axon.
   */
  double
  get_axonal_delay() const
  {
    throw IllegalConnection( "Connection does not support axonal delays." );
  }

  bool
  supports_axonal_delay() const
  {
    return false;
  }

  /**
   * Return the delay of the connection in steps
   */
  long
  get_delay_steps() const
  {
    return delay_;
  }

  /**
   * Set the delay of the connection
   */
  void
  set_delay( const double delay )
  {
    delay_ = Time::delay_ms_to_steps( delay );
  }

  /**
   * Set the delay of the connection in steps
   */
  void
  set_delay_steps( const long delay )
  {
    delay_ = delay;
  }

  long
  get_label() const
  {
    return UNLABELED_CONNECTION;
  }

  /**
   * Triggers an update of a synaptic weight. This function is needed for neuromodulated synaptic plasticity.
   */
  void trigger_update_weight( const thread,
    const std::vector< spikecounter >&,
    const double,
    const CommonSynapseProperties& );

  /**
   * Disables the connection.
   *
   * @see is_disabled
   */
  void
  disable()
  {
    // syn_id_delay_.disable();
  }

  /**
   * Returns a flag denoting if the connection is disabled.
   *
   * @see disable
   */
  bool
  is_disabled() const
  {
    return false;
    // return syn_id_delay_.is_disabled();
  }

protected:
  /* the order of the members below is critical as it influences the size of the object. Please leave unchanged as
     targetidentifierT target_;
     SynIdDelay syn_id_delay_;
  */
  // targetidentifierT target_;
  //! syn_id (9 bit), delay (21 bit) in timesteps of this connection and more_targets and disabled flags (each 1 bit)
  // SynIdDelay syn_id_delay_;
  double delay_;  // TODO JV: Only store delay in delaygroup in node
};

inline void
Connection::get_status( DictionaryDatum& d ) const
{
  def< double >( d, names::delay, get_delay() );
}

inline void
Connection::set_status( const DictionaryDatum& d, ConnectorModel& )
{
  double delay;
  if ( updateValue< double >( d, names::delay, delay ) )
  {
    kernel().connection_manager.get_delay_checker().assert_valid_delay_ms( delay );
    set_delay( delay );
  }
}

inline void
Connection::check_synapse_params( const DictionaryDatum& ) const
{
}

inline void
Connection::calibrate( const TimeConverter& tc )
{
  Time t = tc.from_old_steps( delay_ );
  delay_ = t.get_steps();

  if ( delay_ == 0 )
  {
    delay_ = 1;
  }
}

inline void
Connection::correct_synapse_stdp_ax_delay( const double,
  double*,
  const double,
  const synindex,
  const CommonSynapseProperties&,
  Node* )
{
  throw IllegalConnection( "Connection does not support correction in case of STDP with predominantly axonal delays." );
}

inline void
Connection::process_post_synaptic_spike( const double t_syn, const CommonSynapseProperties& )
{
  throw IllegalConnection( "Connection does not support updates that are triggered by a volume transmitter." );
}

inline void
Connection::trigger_update_weight( const thread,
  const std::vector< spikecounter >&,
  const double,
  const CommonSynapseProperties& )
{
  throw IllegalConnection( "Connection does not support updates that are triggered by a volume transmitter." );
}

} // namespace nest

#endif // CONNECTION_H
