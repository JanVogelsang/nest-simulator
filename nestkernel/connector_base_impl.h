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
inline void
Connector< ConnectionT >::update_stdp_connections( const double post_spike_time_syn,
  const delay dendritic_delay,
  const ConnectorModel* cm )
{
  typename ConnectionT::CommonPropertiesType const& cp =
    static_cast< const GenericConnectorModel< ConnectionT >* >( cm )->get_common_properties();

  const double eps = kernel().connection_manager.get_stdp_eps();
  const auto begin = C_.begin();
  // TODO JV (pt): Verify for precise spike times
  const auto delay_region_it = std::lower_bound( dendritic_delay_regions_.cbegin(), dendritic_delay_regions_.cend(), dendritic_delay, []( const DelayRegion& lhs, const delay d ) -> const bool { return lhs.dendritic_delay < d; });
  if ( delay_region_it != dendritic_delay_regions_.end() and delay_region_it->dendritic_delay == dendritic_delay ) // check if there are connections with given dendritic delay
  {
    const auto end = begin + delay_region_it->end;
    for ( auto it = begin + delay_region_it->start; it != end; ++it )
    {
      // Check if synapse has been updated to this point in time already and ignore the post-spike if that is the case
      const double last_presynaptic_spike_time_ms = it->get_last_presynaptic_spike();
      if ( last_presynaptic_spike_time_ms + eps < post_spike_time_syn )
      {
        it->process_post_synaptic_spike( post_spike_time_syn, cp );
        // TODO JV (pt): Post-synaptic spikes can be sent as weight recorder events as well now as the weight is always
        //  correct now.
        // send_weight_event( tid, std::distance( begin, it ), e, cp, target );
      }
    }
  }
}

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
    wr_e.set_sender_node_id( sources_[ local_target_connection_id ] );
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
double
Connector< ConnectionT >::get_trace( const double pre_spike_time,
  const double dendritic_delay,
  const double tau_minus_inv,
  const std::deque< double >::const_iterator history_begin,
  const std::deque< double >::const_iterator history_end )
{
  const double eps = kernel().connection_manager.get_stdp_eps();

  auto [ last_post_spike, Kminus ] = get_Kminus( dendritic_delay );
  for ( auto it = history_begin; it != history_end; ++it )
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

template < typename ConnectionT >
void
Connector< ConnectionT >::prepare_connections( const thread tid, const index target_lid )
{
  // TODO JV: Benchmark if using std::array would be faster here (stack vs. heap)
  std::vector< ConnectionT > temp_connections;
  std::vector< index > temp_sources;
  temp_connections.reserve( C_.size() );
  temp_sources.reserve( sources_.size() );
  for ( const auto& region : connection_indices_by_delay_ )
  {
    index region_start = static_cast< index >( std::distance( temp_sources.cbegin(), temp_sources.cend() ) );
    for ( index idx : region.second )
    {
#ifdef USE_ADJACENCY_LIST
      // Add an entry to the adjacency list with the new (sorted) index
      kernel().connection_manager.add_adjacency_list_target(
        tid, syn_id_, sources_[ idx ], target_lid, temp_connections.size(), axonal_delays_[ idx ] );
#endif

      temp_connections.push_back( C_.at( idx ) );
      temp_sources.push_back( sources_.at( idx ) );
    }
    dendritic_delay_regions_.emplace_back( region.first, region_start, static_cast< index >( std::distance( temp_sources.cbegin(), temp_sources.cend() ) ) );
  }
  C_.swap( temp_connections );
  sources_.swap( temp_sources );
  std::map< delay, std::vector< index > >().swap( connection_indices_by_delay_ );
#ifdef USE_ADJACENCY_LIST
  std::vector< delay >().swap( axonal_delays_ );
#endif
  std::sort( dendritic_delay_regions_.begin(), dendritic_delay_regions_.end(), []( const DelayRegion& lhs, const DelayRegion& rhs ) -> const bool { return lhs.dendritic_delay < rhs.dendritic_delay; } );
}

} // of namespace nest

#endif
