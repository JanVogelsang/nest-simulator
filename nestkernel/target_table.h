/*
 *  target_table.h
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

#ifndef TARGET_TABLE_H
#define TARGET_TABLE_H

// C++ includes:
#include <cassert>
#include <iostream>
#include <map>
#include <vector>

// Includes from nestkernel:
#include "nest_types.h"
#include "spike_data.h"
#include "target.h"
#include "target_data.h"

namespace nest
{

/**
 * This data structure stores all targets of the local neurons.
 *
 * This is the presynaptic part of the connection infrastructure.
 */
class TargetTable
{
private:
  /**
   * Stores targets of local neurons
   *
   * Four dimensional objects:
   *   - first dim: threads
   *   - second dim: local neurons
   *   - third dim: synapse type  // TODO JV: Reconsider this
   *   - fourth dim: targets
   */
  // std::vector< std::vector< std::vector< std::vector< Target > > > > targets_;
  std::vector< std::vector< std::vector< Target > > > targets_;

public:
  /**
   * Initializes data structures.
   */
  void initialize();

  /**
   * Deletes data structure.
   */
  void finalize();

  /**
   * Adjusts targets_ to number of local nodes.
   */
  void prepare( const size_t tid );

  /**
   * Adds entry to targets_.
   */
  void add_target( const size_t tid, const size_t target_rank, const TargetData& target_data );

  /**
   * Returns all targets of a neuron. Used for filling
   * EventDeliveryManager::emitted_spikes_register_.
   */
  // const std::vector< std::vector< Target > >& get_targets( const size_t tid, const size_t lid ) const;
  const std::vector< Target >& get_targets( const size_t tid, const size_t lid ) const;

  /**
   * Clears all entries of targets_.
   */
  void clear( const size_t tid );
};

// inline const std::vector< std::vector< Target > >&
inline const std::vector< Target >&
TargetTable::get_targets( const size_t tid, const size_t lid ) const
{
  return targets_[ tid ][ lid ];
}

inline void
TargetTable::clear( const size_t tid )
{
  targets_[ tid ].clear();
}

} // namespace nest

#endif /* #ifndef TARGET_TABLE_H */
