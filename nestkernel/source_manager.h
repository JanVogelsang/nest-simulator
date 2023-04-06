/*
 *  source_manager.h
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

#ifndef SOURCE_MANAGER_H
#define SOURCE_MANAGER_H

// C++ includes:
#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

// Includes from libnestutil:
#include "manager_interface.h"

// Includes from nestkernel:
#include "mpi_manager.h"
#include "nest_types.h"
#include "per_thread_bool_indicator.h"
#include "source_table_position.h"
#include "spike_data.h"

namespace nest
{

class TargetData;

/**
 * Manages information about presynaptic neurons during postsynaptic connection creation, before the connection
 * information has been transferred to the presynaptic side. After all connections have been created, the information
 * stored by this structure in the neurons is transferred to the presynaptic side and can be cleared.
 */
class SourceManager : public ManagerInterface
{
private:
  /**
   * Whether the source table has been deleted on all nodes.
   */
  PerThreadBoolIndicator is_cleared_;

  //! Needed during readout of sources_.
  std::vector< SourceTablePosition > current_positions_;

public:
  SourceManager();
  ~SourceManager() override;

  void initialize() override;
  void finalize() override;
  void change_number_of_threads() override;
  void set_status( const DictionaryDatum& ) override;
  void get_status( DictionaryDatum& ) override;

  /**
   * Clears sources_.
   */
  void clear( const thread tid );

  /**
   * Returns true if sources_ has been cleared.
   */
  bool is_cleared() const;

#ifndef USE_ADJACENCY_LIST
  /**
   * Returns if the current target position points to a valid entry.
   */
  bool is_next_target_valid( const thread source_rank );

  /**
   * Returns the next target data, according to the current_positions_.
   */
  void get_next_target_data( const thread source_rank, TargetData& next_target_data );

  /**
   * Resets saved_positions_ to end of sources_.
   */
  void reset_entry_points();
#endif

  /**
   * Returns the source node ID corresponding to the connection ID of target node.
   */
  index get_node_id( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_connection_id ) const;

  /**
   * Determines maximal saved_positions_ after which it is safe to
   * delete sources during clean().
   */
  SourceTablePosition find_maximal_position() const;

  /**
   * Removes all entries marked as processed.
   */
  void clean( const thread tid );

  /**
   * Computes MPI buffer positions for unique combination of source
   * node ID and synapse type across all threads for all secondary
   * connections.
   */
  void compute_buffer_pos_for_unique_secondary_sources( const thread tid,
    std::map< index, size_t >& buffer_pos_of_source_node_id_syn_id_ );

  /**
   * Finds the first entry in sources_ at the given thread id and
   * synapse type that is equal to snode_id.
   */
  index find_first_source( const thread tid, const synindex syn_id, const index snode_id ) const;

  /**
   * Marks entry in sources_ at given position as disabled.
   */
  void disable_connection( const thread tid, const synindex syn_id, const index lcid );

  /**
   * Removes all entries from sources_ that are marked as disabled.
   */
  index remove_disabled_sources( const thread tid, const synindex syn_id );

  /**
   * Returns the number of unique node IDs for given thread id and
   * synapse type in sources_. This number corresponds to the number
   * of targets that need to be communicated during construction of
   * the presynaptic connection infrastructure.
   */
  // size_t num_unique_sources( const thread tid, const synindex syn_id ) const;

  /**
   * Marks that this thread has at least one target from the given source.
   */
  //  void
  //  add_source( const thread tid, const index snode_id )
  //  {
  //    if ( has_source_[ tid ].size() <= snode_id )
  //    {
  //      // Adds as many entries as required to cover all sources up until
  //      has_source_[ tid ].resize( snode_id + 1 ); // default initialized to false
  //    }
  //
  //    has_source_[ tid ][ snode_id ] = true;
  //  }

  /**
   * Encodes combination of node ID and synapse types as single
   * long number.
   */
  index pack_source_node_id_and_syn_id( const index source_node_id, const synindex syn_id ) const;
};

inline index
SourceManager::find_first_source( const thread tid, const synindex syn_id, const index snode_id ) const
{
  assert( false ); // TODO JV (pt): Structural plasticity

  /*// binary search in sorted sources
  const std::vector< Source >::const_iterator begin = sources_[ tid ][ syn_id ].begin();
  const std::vector< Source >::const_iterator end = sources_[ tid ][ syn_id ].end();
  std::vector< Source >::const_iterator it = std::lower_bound( begin, end, Source( snode_id, true ) );

  // source found by binary search could be disabled, iterate through sources until a valid one is found
  while ( it != end )
  {
    if ( it->get_node_id() == snode_id and not it->is_disabled() )
    {
      const index lcid = it - begin;
      return lcid;
    }
    ++it;
  }

  // no enabled entry with this snode ID found
  return invalid_index;*/
}

inline void
SourceManager::disable_connection( const thread tid, const synindex syn_id, const index lcid )
{
  assert( false ); // TODO JV (pt): Structural plasticity
  // disabling a source changes its node ID to 2^62 -1
  // source here
  /*assert( not sources_[ tid ][ syn_id ][ lcid ].is_disabled() );
  sources_[ tid ][ syn_id ][ lcid ].disable();*/
}

// inline size_t
// SourceManager::num_unique_sources( const thread tid, const synindex syn_id ) const
//{
//   return std::count( has_source_[ tid ].begin(), has_source_[ tid ].end(), true );
// }

inline index
SourceManager::pack_source_node_id_and_syn_id( const index source_node_id, const synindex syn_id ) const
{
  assert( source_node_id < 72057594037927936 ); // TODO JV (pt): What is this random number??
  assert( syn_id < invalid_synindex );
  // syn_id is maximally 256, so shifting node ID by 8 bits and storing
  // syn_id in the lowest 8 leads to a unique number
  return ( source_node_id << 8 ) + syn_id;
}

} // namespace nest

#endif // SOURCE_MANAGER_H
