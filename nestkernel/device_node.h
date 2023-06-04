/*
 *  device_node.h
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

#ifndef DEVICE_NODE_H
#define DEVICE_NODE_H

// Includes from nestkernel:
#include "node.h"

namespace nest
{

/**
 * Base class for device objects.
 */
class DeviceNode : public Node
{

public:
  DeviceNode()
    : Node()
    , local_device_id_( invalid_index )
  {
  }

  DeviceNode( DeviceNode const& dn )
    : Node( dn )
    , local_device_id_( invalid_index )
  {
  }

  void set_local_device_id( const index ldid ) override;
  index get_local_device_id() const override;

  // TODO JV (pt): The regular deliver_event function has to be empty for devices, so "remote" spikes won't be sent to
  //  devices again, after they have been sent explicitly to local devices when the event was initially emitted.
  //  A cleaner way would be to move devices out of NodeCollections entirely and have explicit containers for devices
  //  which are separate from regular nodes.
  template < typename EventT >
  void
  deliver_event_to_device( const thread tid,
    const synindex syn_id,
    const index local_target_connection_id,
    const ConnectorModel* cm,
    EventT& e )
  {
    // Send the event to the connection over which this event is transmitted to the node. The connection modifies the
    // event by adding a weight and optionally updates its internal state as well.
    connections_[ syn_id ]->send( tid, local_target_connection_id, 0, cm, e, this );

    // TODO JV (pt): Optionally, the rport can be set here (somehow). For example by just handing it as a parameter to
    //  handle, or just handing the entire local connection id to the handle function.

    handle( e );
  }

#ifdef TIMER_DETAILED
  void
  deliver_event( const synindex, const index, const ConnectorModel*, const Time lag, const delay d, const double offset, const delay min_delay, Stopwatch&, Stopwatch&, Stopwatch& ) override
#else
  void
  deliver_event( const synindex, const index, const ConnectorModel*, const Time lag, const delay d, const double offset, const delay min_delay ) override
#endif
  {
  }

  virtual void event_hook( SpikeEvent& e ) {};
  virtual void event_hook( CurrentEvent& e ) {};

protected:
  index local_device_id_;
};

inline void
DeviceNode::set_local_device_id( const index ldid )
{
  local_device_id_ = ldid;
}

inline index
DeviceNode::get_local_device_id() const
{
  return local_device_id_;
}

} // namespace

#endif /* #ifndef DEVICE_NODE_H */
