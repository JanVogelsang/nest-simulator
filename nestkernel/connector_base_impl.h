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
  Event& e,
  const CommonSynapseProperties& cp,
  const index target_node_id )
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
    // wr_e.set_sender_node_id( sources_[ local_target_connection_id ] );
    wr_e.set_weight( e.get_weight() );
    wr_e.set_delay_steps( e.get_delay_steps() );
    // Set weight_recorder as receiver
    // index wr_node_id = cp.get_wr_node_id();
    // Node* wr_node = kernel().node_manager.get_node_or_proxy( wr_node_id, tid );
    // Put the node_id of the postsynaptic node as receiver node ID
    wr_e.set_receiver_node_id( target_node_id );
    DeviceNode* wr_node =
      static_cast< DeviceNode* >( kernel().node_manager.get_node_or_proxy( cp.get_wr_node_id(), tid ) );
    wr_node->handle( wr_e );
  }
}

template < typename ConnectionT >
void
DendriticDelayConnector< ConnectionT >::update_stdp_connections( const index node_id,
  const double post_spike_time_syn,
  const delay dendritic_delay,
  const ConnectorModel* cm )
{
  typename ConnectionT::CommonPropertiesType const& cp =
    static_cast< const GenericConnectorModel< ConnectionT >* >( cm )->get_common_properties();

  const double eps = kernel().connection_manager.get_stdp_eps();
  // TODO JV (pt): Consider moving the dendritic delay to the regions, sort the regions by dendritic delay after all
  //  regions are known, and do a binary search here instead of a dictionary lookup. Will remove one indirection and
  //  lookup in small arrays is usually faster than in dictionaries.
  const auto delay_region_it = dendritic_delay_region_indices_.find( dendritic_delay );
  // check if there are connections with given dendritic delay
  if ( delay_region_it != dendritic_delay_region_indices_.end() )
  {
    for ( ConnectionT& conn : dendritic_delay_regions_[ delay_region_it->second ].connections )
    {
      // Check if synapse has been updated to this point in time already and ignore the post-spike if that is the case
      const double last_presynaptic_spike_time_ms = conn.get_last_presynaptic_spike();
      if ( last_presynaptic_spike_time_ms + eps < post_spike_time_syn )
      {
        conn.process_post_synaptic_spike( post_spike_time_syn, cp );
        // TODO JV (pt): Post-synaptic spikes can be sent as weight recorder events as well now as the weight is always
        //  correct.
        // send_weight_event( tid, std::distance( begin, it ), e, cp, target );
      }
    }
  }
}

template < typename ConnectionT >
void
DendriticDelayConnector< ConnectionT >::send_weight_event( const thread tid,
  Event& e,
  const CommonSynapseProperties& cp,
  const index target_node_id )
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
    // wr_e.set_sender_node_id( sources_[ local_target_connection_id ] );  // TODO JV (pt)
    wr_e.set_sender_node_id( DISABLED_NODE_ID );
    wr_e.set_weight( e.get_weight() );
    wr_e.set_delay_steps( e.get_delay_steps() );
    // Set weight_recorder as receiver
    // index wr_node_id = cp.get_wr_node_id();
    // Node* wr_node = kernel().node_manager.get_node_or_proxy( wr_node_id, tid );
    // Put the node_id of the postsynaptic node as receiver node ID
    wr_e.set_receiver_node_id( target_node_id );
    DeviceNode* wr_node =
      static_cast< DeviceNode* >( kernel().node_manager.get_node_or_proxy( cp.get_wr_node_id(), tid ) );
    wr_node->handle( wr_e );
  }
}

template < typename ConnectionT >
double
DendriticDelayConnector< ConnectionT >::get_trace( const double pre_spike_time,
  const size_t dendritic_delay_id,
  const double tau_minus_inv,
  const std::deque< double >& history )
{
  // get the post-synaptic trace in synapse time, i.e. the trace of the neuron "dendritic delay"-milliseconds ago
  const double eps = kernel().connection_manager.get_stdp_eps();
  double last_post_spike = dendritic_delay_regions_[ dendritic_delay_id ].last_post_spike;
  double Kminus = dendritic_delay_regions_[ dendritic_delay_id ].Kminus;
  const double dendritic_delay =
    Time::delay_steps_to_ms( dendritic_delay_regions_[ dendritic_delay_id ].dendritic_delay );

  for ( auto it = history.cbegin(); it != history.cend(); ++it )
  {
    const double t_post = *it + dendritic_delay;

    // skip post-synaptic spikes which have been processed by corresponding dendritic delay group already
    if ( t_post < last_post_spike or std::abs( t_post - last_post_spike ) < eps )
    {
      continue;
    }

    // only update trace until finding a spike that occurs at the same time or after the current pre-synaptic spike
    if ( t_post > pre_spike_time or std::abs( t_post - pre_spike_time ) < eps )
    {
      break;
    }

    // only add 1 to trace if post-synaptic spike happened strictly before pre-synaptic spike
    Kminus = Kminus * std::exp( ( last_post_spike - t_post ) * tau_minus_inv ) + 1;
    last_post_spike = t_post;
  }
  Kminus *= std::exp( ( last_post_spike - pre_spike_time ) * tau_minus_inv );

  return Kminus;
}

} // of namespace nest

#endif
