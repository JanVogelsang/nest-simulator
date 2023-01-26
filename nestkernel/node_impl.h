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

#ifndef NEST_NODE_IMPL_H
#define NEST_NODE_IMPL_H

#include "node.h"

// Includes from nestkernel:
#include "connector_base.h"

namespace nest
{

template < typename ConnectionT >
void
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
void
Node::add_connection( Node& source_node,
  const synindex syn_id,
  ConnectionT& connection,
  const rport receptor_type,
  const bool is_primary,
  typename ConnectionT::CommonPropertiesType const& cp )
{
  // Add to connection table
  ConnectorBase* connector = connections_[ syn_id ];

  if ( not connector )
  {
    // No homogeneous Connector with this syn_id exists, we need to create a new homogeneous Connector.
    connections_[ syn_id ] = new Connector< ConnectionT >( syn_id );
  }

  assert( connector );

  const Source src( source_node.get_node_id(), is_primary );
  Connector< ConnectionT >* vc = static_cast< Connector< ConnectionT >* >( connector );
  vc->add_connection( connection, src );

}

}

#endif // NEST_NODE_IMPL_H
