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
#include "connection_type_enum.h"
#include "connector_base.h"

namespace nest
{

template < typename ConnectionT >
inline void
Node::check_connection( Node& source,
  ConnectionT& connection,
  const synindex syn_id,
  const rport receptor_type,
  const delay total_delay,
  const typename ConnectionT::CommonPropertiesType& cp )
{
  // Does the target accept the event type sent by source
  // try to send event from source to target
  // this returns the port of the incoming connection
  // p must be stored in the base class connection
  // this line might throw an exception
  // TODO JV (pt): Add rport to neurons who need it, but not to all connections
  // set_rport();
  source.send_test_event( *this, receptor_type, syn_id );

  // Do the events sent by source mean the same thing as they are interpreted in target?
  // Note that we here use a bitwise and operation (&), because we interpret each bit in the signal type as a collection
  // of individual flags
  if ( not( source.sends_signal() & receives_signal() ) )
  {
    throw IllegalConnection( "Source and target neuron are not compatible (e.g., spiking vs binary neuron)." );
  }
}

template < typename ConnectionT >
const std::pair< index, size_t >
Node::add_connection( Node& source_node,
  const synindex syn_id,
  ConnectionT& connection,
  const rport receptor_type,
  const bool is_primary,
  const ConnectionType connection_type,
  const delay axonal_delay,
  const delay dendritic_delay )
{
  // TODO JV (pt): The whole axonal delay idea needs to be rethought when adding connections, as it depends on the
  //  synapse type if axonal delay has to be provided or not.
  // Check if the source of the connection is a device to add the connection to the corresponding container
  if ( connection_type == ConnectionType::CONNECT_FROM_DEVICE )
  {
    if ( not connections_from_devices_.at( syn_id ) )
    {
      // No homogeneous Connector with this syn_id exists, we need to create a new homogeneous Connector.
      connections_from_devices_.at( syn_id ) = std::make_unique< Connector< ConnectionT > >( syn_id, 1 );
    }
    Connector< ConnectionT >* vc =
      static_cast< Connector< ConnectionT >* >( connections_from_devices_.at( syn_id ).get() );
    return { vc->add_connection( connection ), 0 };
  }
  else
  {
    //    if ( not connections_.at( syn_id ) )
    //    {
    //      // No homogeneous Connector with this syn_id exists, we need to create a new homogeneous Connector.
    //      connections_.at( syn_id ) = std::make_unique< Connector< ConnectionT > >( syn_id );
    //    }
    //    Connector< ConnectionT >* vc = static_cast< Connector< ConnectionT >* >( connections_.at( syn_id ).get() );
    //    if ( connection_type == ConnectionType::CONNECT_TO_DEVICE )
    //    {
    //      return vc->add_connection( connection, source_node.get_node_id() );
    //    }
    //    else
    //    {
    //      vc->add_connection( connection, source_node.get_node_id(), axonal_delay, dendritic_delay );
    //      return invalid_index; // TODO JV (pt): This index should never be used as it will change after sorting
    //    }

    if ( connection_type == ConnectionType::CONNECT_TO_DEVICE )
    {
      if ( not connections_.at( syn_id ) )
      {
        // No homogeneous Connector with this syn_id exists, we need to create a new homogeneous Connector.
        connections_.at( syn_id ) = std::make_unique< Connector< ConnectionT > >( syn_id, 1 );
      }
      Connector< ConnectionT >* vc = static_cast< Connector< ConnectionT >* >( connections_.at( syn_id ).get() );
      return { vc->add_connection( connection ), 0 };
    }
    else
    {
      const ConnectorModel* cm = kernel().model_manager.get_connection_models( thread_ )[ syn_id ];
      if ( not connections_.at( syn_id ) )
      {
        // TODO JV: Benchmark this both with and without reserve
        size_t num;
        if ( syn_id == 55 )
        {
          num = 2250;
        }
        else
        {
          num = 9000;
        }
        if ( kernel().connection_manager.reserve_connections )
        {
          num = 1;
        }

        if ( cm->requires_postponed_delivery() )
        {
          connections_.at( syn_id ) = std::make_unique< DendriticDelayConnector< ConnectionT > >( syn_id, num );
        }
        else
        {
          connections_.at( syn_id ) = std::make_unique< Connector< ConnectionT > >( syn_id, num );
        }
      }
      if ( cm->requires_postponed_delivery() )
      {
        register_stdp_connection( axonal_delay, dendritic_delay, syn_id );
        DendriticDelayConnector< ConnectionT >* vc =
          static_cast< DendriticDelayConnector< ConnectionT >* >( connections_.at( syn_id ).get() );
        return vc->add_connection( connection, dendritic_delay );
      }
      else
      {
        Connector< ConnectionT >* vc = static_cast< Connector< ConnectionT >* >( connections_.at( syn_id ).get() );
        return { vc->add_connection( connection ), 0 };
      }
    }
  }
}

template <>
inline void
Node::deliver_event_from_device< SpikeEvent >( const thread tid,
  const synindex syn_id,
  const index local_target_connection_id,
  const delay d,
  const ConnectorModel* cm,
  SpikeEvent& e )
{
  connections_from_devices_[ syn_id ]->send( tid, node_id_, local_target_connection_id, d, cm, e );
  static_cast< DeviceNode& >( e.get_sender() ).event_hook( e );
  // TODO JV (pt): Make this cleaner, as only needed for poisson generators probably
  if ( e.get_multiplicity() )
  {
    handle( e );
  }
}

template <>
inline void
Node::deliver_event_from_device< CurrentEvent >( const thread tid,
  const synindex syn_id,
  const index local_target_connection_id,
  const delay d,
  const ConnectorModel* cm,
  CurrentEvent& e )
{
  connections_from_devices_[ syn_id ]->send( tid, node_id_, local_target_connection_id, d, cm, e );
  static_cast< DeviceNode& >( e.get_sender() ).event_hook( e );
  handle( e );
}

template < typename EventT >
inline void
Node::deliver_event_from_device( const thread tid,
  const synindex syn_id,
  const index local_target_connection_id,
  const delay d,
  const ConnectorModel* cm,
  EventT& e )
{
  assert( false ); // TODO JV (pt)
}

}

#endif // NODE_IMPL_H
