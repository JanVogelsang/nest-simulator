/*
 *  event_delivery_manager_impl.h
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

#ifndef EVENT_DELIVERY_MANAGER_IMPL_H
#define EVENT_DELIVERY_MANAGER_IMPL_H

#include "event_delivery_manager.h"

// Includes from nestkernel:
#include "connection_manager_impl.h"
#include "kernel_manager.h"

namespace nest
{

template < class EventT >
inline void
EventDeliveryManager::send_local_( Node& source, EventT& e, const long lag )
{
  assert( not source.has_proxies() );
  e.set_stamp( kernel().simulation_manager.get_slice_origin() + Time::step( lag + 1 ) );
  e.set_sender( source );
  const thread t = source.get_thread();
  const index ldid = source.get_local_device_id();
  kernel().connection_manager.send_from_device( t, ldid, e );
}

inline void
EventDeliveryManager::send_local_( Node& source, SecondaryEvent& e )
{
  assert( not source.has_proxies() );
  e.set_stamp( kernel().simulation_manager.get_slice_origin() + Time::step( 1 ) );
  e.set_sender( source );
  const thread t = source.get_thread();
  const index ldid = source.get_local_device_id();
  kernel().connection_manager.send_from_device( t, ldid, e );
}

template < class EventT >
inline void
EventDeliveryManager::send( Node& source, EventT& e, const long lag )
{
  send_local_( source, e, lag );
}

template <>
inline void
EventDeliveryManager::send< SpikeEvent >( Node& source, SpikeEvent& e, const long lag )
{
  const thread tid = source.get_thread();
  e.set_sender_node_id( source.get_node_id() );
  if ( source.has_proxies() )
  {
    local_spike_counter_[ tid ] += e.get_multiplicity();

    e.set_stamp( kernel().simulation_manager.get_slice_origin() + Time::step( lag + 1 ) );
    e.set_sender( source );

    if ( source.is_off_grid() )
    {
      send_off_grid_remote( tid, e, lag );
    }
    else
    {
      send_remote( tid, e, lag );
    }
    kernel().connection_manager.send_to_devices( tid, source.get_thread_lid(), e );
  }
  else
  {
    send_local_( source, e, lag );
  }
}

template <>
inline void
EventDeliveryManager::send< DSSpikeEvent >( Node& source, DSSpikeEvent& e, const long lag )
{
  e.set_sender_node_id( source.get_node_id() );
  send_local_( source, e, lag );
}

inline void
EventDeliveryManager::send_remote( thread tid, SpikeEvent& e, const long lag )
{
  // Put the spike in a buffer for the remote machines
  const index lid = kernel().vp_manager.node_id_to_lid( e.get_sender().get_node_id() );
  const std::vector< Target >& targets = kernel().connection_manager.get_remote_targets_of_local_node( tid, lid );

  for ( std::vector< Target >::const_iterator it = targets.begin(); it != targets.end(); ++it )
  {
    const thread assigned_tid = ( *it ).get_rank() / kernel().vp_manager.get_num_assigned_ranks_per_thread();

    // Unroll spike multiplicity as plastic synapses only handle individual spikes.
    for ( int i = 0; i < e.get_multiplicity(); ++i ) // TODO JV (pt): Remove multiplicity
    {
      emitted_spikes_register_[ tid ][ assigned_tid ][ lag ].push_back( *it );
    }
  }
}

inline void
EventDeliveryManager::send_off_grid_remote( thread tid, SpikeEvent& e, const long lag )
{
  // Put the spike in a buffer for the remote machines
  const index lid = kernel().vp_manager.node_id_to_lid( e.get_sender().get_node_id() );
  const std::vector< Target >& targets = kernel().connection_manager.get_remote_targets_of_local_node( tid, lid );

  for ( std::vector< Target >::const_iterator it = targets.begin(); it != targets.end(); ++it )
  {
    const thread assigned_tid = ( *it ).get_rank() / kernel().vp_manager.get_num_assigned_ranks_per_thread();

    // Unroll spike multiplicity as plastic synapses only handle individual spikes.
    for ( int i = 0; i < e.get_multiplicity(); ++i )
    {
      off_grid_emitted_spike_register_[ tid ][ assigned_tid ][ lag ].push_back( OffGridTarget( *it, e.get_offset() ) );
    }
  }
}

inline void
EventDeliveryManager::send_secondary( Node& source, SecondaryEvent& e )
{
  const thread tid = kernel().vp_manager.get_thread_id();
  const index lid = kernel().vp_manager.node_id_to_lid( source.get_node_id() );

  if ( source.has_proxies() )
  {

    // We need to consider every synapse type this event supports to
    // make sure also labeled and connection created by CopyModel are
    // considered.
    const std::set< synindex >& supported_syn_ids = e.get_supported_syn_ids();
    for ( const auto& syn_id : supported_syn_ids )
    {
      const std::vector< size_t >& positions =
        kernel().connection_manager.get_secondary_send_buffer_positions( tid, lid, syn_id );

      for ( size_t i = 0; i < positions.size(); ++i )
      {
        std::vector< unsigned int >::iterator it = send_buffer_secondary_events_.begin() + positions[ i ];
        e >> it;
      }
    }
    kernel().connection_manager.send_to_devices( tid, source.get_thread_lid(), e );
  }
  else
  {
    send_local_( source, e );
  }
}

inline size_t
EventDeliveryManager::write_toggle() const
{
  return kernel().simulation_manager.get_slice() % 2;
}

#ifdef USE_ADJACENCY_LIST
inline void
EventDeliveryManager::deliver_to_adjacency_list( const thread tid,
  const index adjacency_list_index,
  SpikeEvent& se,
  const std::vector< ConnectorModel* >& cm )
{
  auto [ adjacency_list_it, adjacency_list_end ] = kernel().connection_manager.get_targets( tid, adjacency_list_index );
  for ( ; adjacency_list_it != adjacency_list_end; ++adjacency_list_it )
  {
    const index local_target_node_id = adjacency_list_it->local_target_node_id;
    const index local_target_connection_id = adjacency_list_it->local_target_connection_id;
    const synindex syn_id = adjacency_list_it->syn_id;
    se.set_sender_node_id_info( tid, syn_id, local_target_node_id, local_target_connection_id );
    Node* target_node = kernel().node_manager.thread_lid_to_node( tid, local_target_node_id );
    target_node->deliver_event( tid, syn_id, local_target_connection_id, cm, se );
    // TODO JV: Let node handle delivery differently when delay is already provided (needs a way to know if total or
    //  only axonal-dendritic delay is sent)
  }
}
#endif

} // of namespace nest

#endif
