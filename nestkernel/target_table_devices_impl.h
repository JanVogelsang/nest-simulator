/*
 *  target_table_devices_impl.h
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

#ifndef TARGET_TABLE_DEVICES_IMPL_H
#define TARGET_TABLE_DEVICES_IMPL_H

// Includes from nestkernel:
#include "connector_base.h"
#include "kernel_manager.h"
#include "node.h"
#include "node_impl.h"
#include "target_table_devices.h"
#include "vp_manager_impl.h"

namespace nest
{

inline void
TargetTableDevices::add_connection_from_device( Node& source,
  Node& target,
  const index local_target_connection_id,
  const thread tid,
  const synindex syn_id )
{
  const index ldid = source.get_local_device_id();
  assert( ldid != invalid_index );
  assert( ldid < targets_from_devices_[ tid ].size() );

  targets_from_devices_[ tid ][ ldid ].push_back( Target(
    target.get_thread(), kernel().vp_manager.get_vp(), syn_id, target.get_thread_lid(), local_target_connection_id ) );

  // store node ID of sending device
  sending_devices_node_ids_[ tid ][ ldid ] = source.get_node_id();
}

inline void
TargetTableDevices::add_connection_to_device( Node& source,
  Node& target,
  const index local_target_connection_id,
  const thread tid,
  const synindex syn_id )
{
  const index source_lid = source.get_thread_lid();
  // Check if the source node has not been registered as source to devices yet
  if ( not targets_to_devices_[ tid ].count( source_lid ) )
  {
    targets_to_devices_[ tid ].emplace( source_lid, std::vector< Target >() );
  }

  targets_to_devices_[ tid ][ source_lid ].push_back( Target(
    target.get_thread(), kernel().vp_manager.get_vp(), syn_id, target.get_thread_lid(), local_target_connection_id ) );
}

// inline void
// nest::TargetTableDevices::get_synapse_status_to_device( const thread tid,
//   const index source_node_id,
//   const synindex syn_id,
//   DictionaryDatum& dict,
//   const index lcid ) const
//{
//   const index lid = kernel().vp_manager.node_id_to_lid( source_node_id );
//   if ( targets_to_devices_[ tid ][ lid ][ syn_id ] )
//   {
//     targets_to_devices_[ tid ][ lid ][ syn_id ]->get_synapse_status( tid, lcid, dict );
//   }
// }

// inline void
// nest::TargetTableDevices::set_synapse_status_to_device( const thread tid,
//   const index source_node_id,
//   const synindex syn_id,
//   ConnectorModel& cm,
//   const DictionaryDatum& dict,
//   const index lcid )
//{
//   const index lid = kernel().vp_manager.node_id_to_lid( source_node_id );
//   if ( targets_to_devices_[ tid ][ lid ][ syn_id ] )
//   {
//     targets_to_devices_[ tid ][ lid ][ syn_id ]->set_synapse_status( lcid, dict, cm );
//   }
// }

template < typename EventT >
inline void
TargetTableDevices::send_from_device( const thread tid, const index ldid, EventT& e )
{
  const std::vector< ConnectorModel* > cm = kernel().model_manager.get_connection_models( tid );
  for ( std::vector< Target >::iterator it = targets_from_devices_[ tid ][ ldid ].begin();
        it != targets_from_devices_[ tid ][ ldid ].end();
        ++it )
  {
    Node* target_node = kernel().node_manager.thread_lid_to_node( tid, it->get_local_target_node_id() );
    target_node->deliver_event_from_device( tid, it->get_syn_id(), it->get_local_target_connection_id(), cm, e );
  }
}

template < typename EventT >
inline void
TargetTableDevices::send_to_devices( const thread tid, const index source_node_lid, EventT& e )
{
  const std::vector< ConnectorModel* > cm = kernel().model_manager.get_connection_models( tid );
  for ( std::vector< Target >::iterator it = targets_to_devices_[ tid ][ source_node_lid ].begin();
        it != targets_to_devices_[ tid ][ source_node_lid ].end();
        ++it )
  {
    DeviceNode* target_node =
      static_cast< DeviceNode* >( kernel().node_manager.thread_lid_to_node( tid, it->get_local_target_node_id() ) );
    target_node->deliver_event_to_device( tid, it->get_syn_id(), it->get_local_target_connection_id(), cm, e );
  }
}

}

#endif /* TARGET_TABLE_DEVICES_IMPL_H */
