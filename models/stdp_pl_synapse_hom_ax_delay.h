/*
 *  stdp_pl_synapse_hom_ax_delay.h
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

#ifndef STDP_PL_SYNAPSE_HOM_AX_DELAY_H
#define STDP_PL_SYNAPSE_HOM_AX_DELAY_H

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "axonal_delay_connection.h"
#include "connection.h"

namespace nest
{

/* BeginUserDocs: synapse, spike-timing-dependent plasticity

Short description
+++++++++++++++++

Synapse type for spike-timing dependent plasticity with power law

Description
+++++++++++

``stdp_pl_synapse`` is a connector to create synapses with spike time
dependent plasticity using homoegeneous parameters (as defined in [1]_).

Parameters
++++++++++

=========  ======  ====================================================
 tau_plus  ms      Time constant of STDP window, potentiation
                   (tau_minus defined in postsynaptic neuron)
 lambda    real    Learning rate
 alpha     real    Asymmetry parameter (scales depressing increments as
                   alpha*lambda)
 mu        real    Weight dependence exponent, potentiation
=========  ======  ====================================================

The parameters can only be set by SetDefaults and apply to all synapses of
the model.

.. warning::

   This synaptic plasticity rule does not take
   :ref:`precise spike timing <sim_precise_spike_times>` into
   account. When calculating the weight update, the precise spike time part
   of the timestamp is ignored.

References
++++++++++

.. [1] Morrison A, Aertsen A, Diesmann M. (2007) Spike-timing dependent
       plasticity in balanced random netrks. Neural Computation,
       19(6):1437-1467. DOI: https://doi.org/10.1162/neco.2007.19.6.1437

Transmits
+++++++++

SpikeEvent

See also
++++++++

stdp_synapse, tsodyks_synapse, static_synapse

EndUserDocs */

/**
 * Class containing the common properties for all synapses of type
 * stdp_pl_synapse_hom_ax_delay.
 */
class STDPPLHomAxDelayCommonProperties : public CommonSynapseProperties
{

public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  STDPPLHomAxDelayCommonProperties();

  /**
   * Get all properties and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  // data members common to all connections
  double tau_plus_;
  double tau_plus_inv_; //!< 1 / tau_plus for efficiency
  double lambda_;
  double alpha_;
  double mu_;
};


/**
 * Class representing an STDP connection with homogeneous parameters, i.e.
 * parameters are the same for all synapses.
 */
class stdp_pl_synapse_hom_ax_delay : public AxonalDelayConnection
{

public:
  typedef STDPPLHomAxDelayCommonProperties CommonPropertiesType;
  typedef Connection ConnectionBase;


  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  stdp_pl_synapse_hom_ax_delay();

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  stdp_pl_synapse_hom_ax_delay( const stdp_pl_synapse_hom_ax_delay& ) = default;
  stdp_pl_synapse_hom_ax_delay& operator=( const stdp_pl_synapse_hom_ax_delay& ) = default;

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay;
  using ConnectionBase::get_delay_steps;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   */
  void send( Event& e, thread t, const STDPPLHomAxDelayCommonProperties&, Node* target );

  /**
   * Process a post-synaptic spike after it is backpropagated to the synapse.
   * @param t_syn The time the post-synaptic spike arrives at this synapse
   */
  void process_post_synaptic_spike( double t_syn, const STDPPLHomAxDelayCommonProperties& cp );

  /**
   * Framework for STDP with predominantly axonal delays:
   * Correct this synapse and the corresponding previously sent spike
   * taking into account a new post-synaptic spike.
   */
  void correct_synapse_stdp_ax_delay( const double t_last_spike,
    double* weight_revert,
    const double t_post_spike,
    const synindex syn_id,
    const STDPPLHomAxDelayCommonProperties& cp,
    Node* target );

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport ) override
    {
      return invalid_port;
    }
  };

  /*
   * This function calls check_connection on the sender and checks if the
   * receiver accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for STDP
   * connections we have to call register_stdp_pl_connection on the target
   * neuron to inform the Archiver to collect spikes for this connection.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
   */
  void
  check_connection( Node& s, Node& t, const rport receptor_type, const synindex syn_id, const CommonPropertiesType& cp )
  {
    ConnTestDummyNode dummy_target;

    const double delay = get_delay();
    if ( get_axonal_delay() > delay )
    {
      throw BadProperty( "Axonal delay should not exceed total synaptic delay." );
    }
    t.register_stdp_connection( get_dendritic_delay(), syn_id );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

  void perform_facilitation( const double t_syn, const STDPPLHomAxDelayCommonProperties& cp );

private:
  double
  facilitate_( double w, double kplus, const STDPPLHomAxDelayCommonProperties& cp )
  {
    return w + cp.lambda_ * std::pow( w, cp.mu_ ) * kplus;
  }

  double
  depress_( double w, double kminus, const STDPPLHomAxDelayCommonProperties& cp )
  {
    double new_w = w - w * cp.lambda_ * cp.alpha_ * kminus;
    return new_w > 0.0 ? new_w : 0.0;
  }

  // data members of each connection
  double weight_;
  double Kplus_;
  double t_lastspike_;
};

//
// Implementation of class stdp_pl_synapse_hom_ax_delay.
//

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 */
inline void
stdp_pl_synapse_hom_ax_delay::send( Event& e, thread, const STDPPLHomAxDelayCommonProperties& cp, Node* target )
{
  const double axonal_delay = get_axonal_delay();
  const double dendritic_delay = get_dendritic_delay();

  const double t_spike = e.get_stamp().get_ms();

  // Todo JV (help): Move this to ConnectorBase. This would require a common base class for all STDP synapses.
  // Process all post-synaptic spikes which have to be processed before this pre-synaptic spike. The post-synaptic
  // spikes could not be processed before, as this connection might have an axonal delay lower than min_delay. Hence,
  // the order at which pre- and post-synaptic spikes have to be processed can only be known after the communication.
  // TODO JV: If the axonal delay is longer than min_delay, this might not be necessary
  auto [ K_minus, post_spikes_to_process ] = target->get_stdp_history(
    t_lastspike_ + axonal_delay, t_spike + axonal_delay, dendritic_delay, e.get_sender_spike_data().get_syn_id() );

  for ( const double t_post_spike_syn : post_spikes_to_process )
  {
    perform_facilitation( t_post_spike_syn, cp );
  }

  // Framework for STDP with predominantly axonal delays:
  // Store weight before depression for potential later correction
  const double weight_revert = weight_;

  // depression due to new pre-synaptic spike
  weight_ = depress_( weight_, K_minus, cp );

  // std::cout << std::setprecision( 17 ) << "Pre " << t_spike + axonal_delay << " " << K_minus << " " << weight_ <<
  // std::endl;

  e.set_weight( weight_ );
  e.set_delay_steps( get_delay_steps() );
  e();

  // TODO JV (pt): This does not take min_delay into account. The pre-synaptic spike could be safe already, due to the
  //  (up to min_delay-1) late delivery.
  if ( ( axonal_delay - dendritic_delay ) > kernel().connection_manager.get_stdp_eps() )
  {
    target->add_correction_entry_stdp_ax_delay(
      dynamic_cast< SpikeEvent& >( e ), t_lastspike_, weight_revert, dendritic_delay );
  }

  Kplus_ = Kplus_ * std::exp( ( t_lastspike_ - t_spike ) * cp.tau_plus_inv_ ) + 1.0;

  t_lastspike_ = t_spike;
}

inline void
stdp_pl_synapse_hom_ax_delay::perform_facilitation( const double t_syn, const STDPPLHomAxDelayCommonProperties& cp )
{
  // facilitation due to postsynaptic spike
  double minus_dt = t_lastspike_ + get_axonal_delay() - t_syn;
  assert( minus_dt < -1.0 * kernel().connection_manager.get_stdp_eps() );
  weight_ = facilitate_( weight_, Kplus_ * std::exp( minus_dt * cp.tau_plus_inv_ ), cp );

  // std::cout << std::setprecision( 17 ) << "Post " << start->t + dendritic_delay << " "
  //           << Kplus_ * std::exp( minus_dt * cp.tau_plus_inv_ ) << " " << weight_ << std::endl;
}

inline void
stdp_pl_synapse_hom_ax_delay::process_post_synaptic_spike( double t_syn, const STDPPLHomAxDelayCommonProperties& cp )
{
  // Todo JV (pt): Move this to ConnectorBase. This would require t_lastspike to be a common member of all STDP synapses.
  // Check if synapse has been updated to this point in time already and ignore the post-spike if that is the case.
  if ( t_lastspike_ + get_axonal_delay() + kernel().connection_manager.get_stdp_eps() < t_syn )
  {
    perform_facilitation( t_syn, cp );
  }
}

inline void
stdp_pl_synapse_hom_ax_delay::correct_synapse_stdp_ax_delay( const double t_last_spike,
  double* weight_revert,
  const double t_post_spike,
  const synindex syn_id,
  const STDPPLHomAxDelayCommonProperties& cp,
  Node* target )
{
  const double t_spike = t_lastspike_; // no new pre-synaptic spike since last send()
  const double wrong_weight = weight_; // incorrectly transmitted weight
  weight_ = *weight_revert;            // removes the last depressive step

  const double axonal_delay = get_axonal_delay();
  double dendritic_delay = get_delay() - axonal_delay;

  // facilitation due to new post-synaptic spike
  const double minus_dt = t_last_spike + axonal_delay - ( t_post_spike + dendritic_delay );

  double K_plus_revert;
  // Multiple pre-spikes before post-spike possible. Only facilitate if not facilitated already (only if first
  // correction for this post-spike)
  if ( minus_dt < -1.0 * kernel().connection_manager.get_stdp_eps() )
  {
    // Kplus value at t_last_spike_ needed
    K_plus_revert = ( Kplus_ - 1.0 ) / std::exp( ( t_last_spike - t_spike ) * cp.tau_plus_inv_ );
    weight_ = facilitate_( weight_, K_plus_revert * std::exp( minus_dt * cp.tau_plus_inv_ ), cp );

    // update weight_revert in case further correction will be required later
    *weight_revert = weight_;
  }

  // depression taking into account new post-synaptic spike
  const double K_minus = target->get_stdp_history( t_last_spike + axonal_delay, t_spike + axonal_delay, dendritic_delay, syn_id ).first;
  weight_ = depress_( weight_, K_minus, cp );

  // std::cout << std::setprecision( 17 ) << "Post " << t_post_spike + dendritic_delay << " " << K_plus_revert << " " <<
  // weight_
  // << std::endl;

  // send a correcting event to the target neuron
  SpikeEvent e;
  e.set_weight( weight_ - wrong_weight );
  e.set_delay_steps( get_delay_steps() );
  e.set_stamp( Time::ms_stamp( t_spike ) );
  e();
}

} // of namespace nest

#endif // of #ifndef STDP_PL_SYNAPSE_HOM_AX_DELAY_H