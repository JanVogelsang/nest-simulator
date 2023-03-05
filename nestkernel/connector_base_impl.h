/*
 *  connector_base_impl.h
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

#include "connector_base.h"

// Includes from nestkernel:
#include "connector_model.h"
#include "kernel_manager.h"

#ifndef CONNECTOR_BASE_IMPL_H
#define CONNECTOR_BASE_IMPL_H

namespace nest
{

template < typename ConnectionT >
void
Connector< ConnectionT >::send_weight_event( const thread tid,
  const index local_target_connection_id,
  Event& e,
  const CommonSynapseProperties& cp,
  Node* target )
{
  // If the pointer to the receiver node in the event is invalid,
  // the event was not sent, and a WeightRecorderEvent is therefore not created.
  if ( cp.get_weight_recorder() )
  {
    // Create new event to record the weight and copy relevant content.
    WeightRecorderEvent wr_e;
    wr_e.set_port( e.get_port() );
    wr_e.set_stamp( e.get_stamp() );
    wr_e.set_sender( e.get_sender() );
    wr_e.set_sender_node_id( sources_[ local_target_connection_id ].get_node_id() );
    wr_e.set_weight( e.get_weight() );
    wr_e.set_delay_steps( e.get_delay_steps() );
    // Set weight_recorder as receiver
    // index wr_node_id = cp.get_wr_node_id();
    // Node* wr_node = kernel().node_manager.get_node_or_proxy( wr_node_id, tid );
    // Put the node_id of the postsynaptic node as receiver node ID
    wr_e.set_receiver_node_id( target->get_node_id() );
    DeviceNode* wr_node =
      static_cast< DeviceNode* >( kernel().node_manager.get_node_or_proxy( cp.get_wr_node_id(), tid ) );
    wr_node->handle( wr_e );
  }
}

template < typename ConnectionT >
void
Connector< ConnectionT >::correct_synapse_stdp_ax_delay( const index local_target_connection_id,
  const double t_last_pre_spike,
  double* weight_revert,
  const double t_post_spike,
  const synindex syn_id,
  Node* target )
{
  typename ConnectionT::CommonPropertiesType const& cp = static_cast< GenericConnectorModel< ConnectionT >* >(
    kernel().model_manager.get_connection_models( target->get_thread() )[ syn_id_ ] )
                                                           ->get_common_properties();
  C_[ local_target_connection_id ].correct_synapse_stdp_ax_delay( t_last_pre_spike,
    weight_revert,
    t_post_spike,
    syn_id,
    get_dendritic_delay( local_target_connection_id ),
    cp,
    target );
}

template < typename ConnectionT >
std::pair< double, std::vector< double > >
Connector< ConnectionT >::get_stdp_history( const double last_pre_spike_time,
  const double pre_spike_time,
  const double dendritic_delay,
  const double tau_minus_inv,
  const std::deque< double >::const_iterator history_begin,
  const std::deque< double >::const_iterator history_end )
{
  const double eps = kernel().connection_manager.get_stdp_eps();

  std::vector< double > post_spikes_to_process;
  auto [ last_post_spike, Kminus ] = get_Kminus( dendritic_delay );
  for ( auto it = history_begin; it != history_end; ++it )
  {
    const double t_post = *it + dendritic_delay;
    // check if post-synaptic spike should be processed after the pre-synaptic one
    if ( t_post > pre_spike_time + eps )
    {
      break;
    }

    // if pre- and post-synaptic spikes occur at same time, don't add post-synaptic spike to trace yet
    if ( std::abs( t_post - pre_spike_time ) < eps )
    {
      post_spikes_to_process.push_back( t_post );
      break;
    }

    // only add 1 to trace if post-synaptic spike happened strictly before pre-synaptic spike
    Kminus = Kminus * std::exp( ( last_post_spike - t_post ) * tau_minus_inv ) + 1;
    last_post_spike = t_post;

    // check if post-synaptic spike has not been processed by the connection yet (influences the trace in any case)
    if ( t_post > last_pre_spike_time + eps )
    {
      post_spikes_to_process.push_back( t_post );
    }
  }
  Kminus *= std::exp( ( last_post_spike - pre_spike_time ) * tau_minus_inv );

  return { Kminus, post_spikes_to_process };
}

} // of namespace nest

#endif
