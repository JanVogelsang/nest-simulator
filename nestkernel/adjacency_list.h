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
#include <vector>

// Includes from libnestutil:
#include "stopwatch.h"

// Includes from nestkernel:
#include "nest_time.h"
#include "nest_types.h"

namespace nest
{
class Node;
class Event;

/**
 *
 */
class AdjacencyListTarget {
  unsigned short neuron_index_;
  unsigned short synapse_index_;
  // delay needs 21 bits
};

/**
 * Maps incoming spikes to thread-local target neurons and the corresponding neuron-local synapse over which the spike
 * reaches the neuron.
 */
class AdjacencyList
{
  /**
   * Stores all targets for each source neuron with targets on this process.
   * Two dimensional object:
   *   - first dim: threads
   *   - second dim: targets of source neuron
   */
  std::vector< std::vector< AdjacencyListTarget > > adjacency_list_;

};

} // nest

#endif // NEST_ADJACENCY_LIST_H
