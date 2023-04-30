/*
 *  event.cpp
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

#include "event.h"

// Includes from nestkernel:
#include "kernel_manager.h"
#include "node.h"

namespace nest
{
Event::Event()
  : sender_node_id_( 0 ) // initializing to 0 as this is an unsigned type
                         // node ID 0 is network, can never send an event, so
                         // this is safe
  , spike_data_()
  , sender_( nullptr )
  , p_( -1 )
  , d_( 1 )
  , stamp_( Time::step( 0 ) )
  , stamp_steps_( 0 )
  , offset_( 0.0 )
  , w_( 0.0 )
{
}

index
Event::retrieve_sender_node_id_from_source_table() const
{
  if ( sender_node_id_ > 0 )
  {
    return sender_node_id_;
  }
  else
  {
    const index node_id = kernel().connection_manager.get_source_node_id(
      spike_data_.tid, spike_data_.syn_id, spike_data_.local_node_id, spike_data_.local_connection_id );
    return node_id;
  }
}

} // namespace nest
