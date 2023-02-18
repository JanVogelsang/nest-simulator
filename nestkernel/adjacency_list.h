/*
 *  adjacency_list.h
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

#ifndef NEST_ADJACENCY_LIST_H
#define NEST_ADJACENCY_LIST_H

// C++ includes:
#include <map>
#include <vector>

// Includes from libnestutil:
#include "stopwatch.h"

// Includes from nestkernel:
#include "nest_time.h"
#include "nest_types.h"
#include "static_assert.h"

namespace nest
{
class Node;
class Event;

/**
 *
 */
struct AdjacencyListTarget
{
  unsigned short local_target_node_id;  // : NUM_BITS_LOCAL_NODE_ID;
  unsigned short local_target_connection_id;  // : NUM_BITS_LOCAL_CONNECTION_ID;
  // delay needs 21 bits, if more than 4 bytes are needed here, change from unsigned int to size_t

  AdjacencyListTarget(const index local_target_node_id, const index local_target_connection_id) : local_target_node_id(local_target_node_id), local_target_connection_id(local_target_connection_id) {}
};

//! check legal size
using success_syn_id_delay_data_size = StaticAssert< sizeof( AdjacencyListTarget ) == 4 >::success;

/**
 * Maps incoming spikes to thread-local target neurons and the corresponding node-local synapse over which the spike
 * reaches the neuron.
 */
class AdjacencyList
{
  /**
   * Intermediate structure to map source neurons to the corresponding index in the adjacency_list for each thread.
   * Deleted after communication of targets to pre-synaptic processes.
   */
  std::vector< std::map< index, index > > sources_;

  /**
   * Stores all targets for each source neuron with targets on this process.
   * Two dimensional object:
   *   - first dim: threads
   *   - second dim: source neuron
   *   - third dim: targets of source neuron
   */
  std::vector< std::vector< std::vector < AdjacencyListTarget > > > adjacency_list_;

public:
  void add_target( const thread tid, const index source_node_id, const index target_node_id, const index target_connection_id );

  void clear_sources()
  {
    std::vector< std::map< index, index > >().swap( sources_ );
  }
};

} // nest

#endif // NEST_ADJACENCY_LIST_H
