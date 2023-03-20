/*
 *  node_impl.h
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

#ifndef NODE_IMPL_H
#define NODE_IMPL_H

#include "node.h"

// Includes from nestkernel:
#include "connector_base.h"
#include "connection_type_enum.h"

namespace nest
{

template < typename ConnectionT >
inline void
Node::check_connection( Node& source, const synindex syn_id, const rport receptor_type )
{
  // 1. does this connection support the event type sent by source
  // try to send event from source to dummy_target
  // this line might throw an exception
  typename ConnectionT::ConnTestDummyNode dummy_node;
  source.send_test_event( dummy_node, receptor_type, syn_id, true );

  // 2. does the target accept the event type sent by source
  // try to send event from source to target
  // this returns the port of the incoming connection
  // p must be stored in the base class connection
  // this line might throw an exception
  // TODO JV (pt): Add rport to neurons who need it, but not to all connections
  // target_.set_rport( source.send_test_event( target, receptor_type, get_syn_id(), false ) );

  // 3. do the events sent by source mean the same thing as they are
  // interpreted in target?
  // note that we here use a bitwise and operation (&), because we interpret
  // each bit in the signal type as a collection of individual flags
  if ( not( source.sends_signal() & receives_signal() ) )
  {
    throw IllegalConnection( "Source and target neuron are not compatible (e.g., spiking vs binary neuron)." );
  }
}

template < typename ConnectionT >
const index
Node::add_connection( Node& source_node,
  const synindex syn_id,
  ConnectionT& connection,
  const rport receptor_type,
  const bool is_primary,
  const ConnectionType connection_type,
  const delay axonal_delay,
  const delay dendritic_delay )
{
  // Check if the source of the connection is a device to add the connection to the corresponding container
  if ( connection_type == ConnectionType::CONNECT_FROM_DEVICE )
  {
    if ( not connections_from_devices_.at( syn_id ) )
    {
      // No homogeneous Connector with this syn_id exists, we need to create a new homogeneous Connector.
      connections_from_devices_.at( syn_id ) = std::make_unique< Connector< ConnectionT > >( syn_id );
    }
    Connector< ConnectionT >* vc = static_cast< Connector< ConnectionT >* >( connections_from_devices_.at( syn_id ).get() );
    return vc->add_device_connection( connection, source_node.get_node_id() );
  }
  else
  {
    if ( not connections_.at( syn_id ) )
    {
      // No homogeneous Connector with this syn_id exists, we need to create a new homogeneous Connector.
      connections_.at( syn_id ) = std::make_unique< Connector< ConnectionT > >( syn_id );
    }
    Connector< ConnectionT >* vc = static_cast< Connector< ConnectionT >* >( connections_.at( syn_id ).get() );
    if ( connection_type == ConnectionType::CONNECT_TO_DEVICE )
    {
      return vc->add_device_connection( connection, source_node.get_node_id() );
    }
    else{
      vc->add_connection( connection, source_node.get_node_id(), axonal_delay, dendritic_delay );
      return invalid_index;  // TODO JV (pt): This index should never be used as it will change after sorting
    }
  }
}

template < typename EventT >
inline void
Node::deliver_event_from_device( const thread tid,
  const synindex syn_id,
  const index local_target_connection_id,
  const delay dendritic_delay,
  const std::vector< ConnectorModel* >& cm,
  EventT& e )
{
  // Send the event to the connection over which this event is transmitted to the node. The connection modifies the
  // event by adding a weight and optionally updates its internal state as well.
  connections_from_devices_[ syn_id ]->send( tid, local_target_connection_id, 0, dendritic_delay, cm, e, this );

  // TODO JV (pt): Optionally, the rport can be set here (somehow). For example by just handing it as a parameter to
  //  handle, or just handing the entire local connection id to the handle function.

  handle( e );
}

}

#endif // NODE_IMPL_H
