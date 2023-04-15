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

#ifndef ADJACENCY_LIST_H
#define ADJACENCY_LIST_H
#ifdef USE_ADJACENCY_LIST

// C++ includes:
#include <cassert>
#include <map>
#include <vector>

// Includes from nestkernel:
#include "nest_types.h"
#include "static_assert.h"

namespace nest
{
/**
 *
 */
struct AdjacencyListTarget
{
  index local_target_node_id : NUM_BITS_LOCAL_NODE_ID;
  index local_target_connection_id : NUM_BITS_LOCAL_CONNECTION_ID;
  synindex syn_id : NUM_BITS_SYN_ID;
  delay axonal_delay : NUM_BITS_DELAY;
  // TODO JV (pt): Still some bits to spare here

  AdjacencyListTarget( const index local_target_node_id,
    const index local_target_connection_id,
    const synindex syn_id,
    const delay axonal_delay )
    : local_target_node_id( local_target_node_id )
    , local_target_connection_id( local_target_connection_id )
    , syn_id( syn_id )
    , axonal_delay( axonal_delay )
  {
  }
};

//! check legal size
using success_adjacency_list_target_size = StaticAssert< sizeof( AdjacencyListTarget ) == 8 >::success;

/**
 * Maps incoming spikes to thread-local target neurons and the corresponding node-local synapse over which the spike
 * reaches the neuron.
 */
class AdjacencyList
{
  /**
   * Stores all targets for each source neuron with targets on this process.
   * Three dimensional object:
   *   - first dim: threads
   *   - second dim: source neuron
   *   - third dim: targets of source neuron
   */
  // TODO JV (pt): Evaluate if adding local target node id as another dimension might make sense for a low number of VPs
  std::vector< std::vector< std::vector< AdjacencyListTarget > > > adjacency_list_;

  /**
   * Intermediate structure to map source neurons to the corresponding index in the adjacency_list for each thread.
   * Deleted after communication of targets to pre-synaptic processes.
   * Three dimensional object:
   *   - first dim: target threads
   *   - second dim: source ranks
   *   - third dim: map of source node id to index in adjacency_list
   */
  std::vector< std::vector< std::map< index, size_t > > > sources_;

  /**
   * A structure to map incoming spikes to the adjacency list entry for each thread on the postsynaptic side if spike
   * compression is enabled. Internally arranged in a 2d structure:
   * sources|(thread->adjacency list index)
   */
  std::vector< std::map< thread, index > > compressed_indices_;

  /**
   * Intermediate two-dimensional structure. For each source rank, map all source node ids to the "compressed" index in
   * the structure which stores all adjacency list indices, one for each thread. See compressed_indices_.
   */
  std::vector< std::map< index, size_t > > source_to_compressed_index_;

  /**
   * For each source rank, save the position of the next pair that is going to be retrieved. Required to communicate
   * target data for compressed spikes.
   */
  std::vector< std::map< index, size_t >::const_iterator > next_compressed_index_;

  /**
   * For each source rank, save the position of the next pair that is going to be retrieved. Required to communicate
   * target data.
   */
  std::vector< std::map< index, index >::const_iterator > next_source_index_;
  /**
   * For each source rank, save the current target thread position. Required to communicate target data.
   */
  std::vector< size_t > next_source_index_thread_;

public:
  void resize( const thread num_threads, const thread num_ranks, const bool compressed );

  void add_target( const thread tid,
    const synindex syn_id,
    const index source_node_id,
    const thread source_rank,
    const index target_node_id,
    const index target_connection_id,
    const delay axonal_delay );

  std::pair< std::vector< AdjacencyListTarget >::const_iterator, std::vector< AdjacencyListTarget >::const_iterator >
  get_iterators( const thread tid, const index adjacency_list_index ) const;

  const std::map< thread, index >& get_compressed_spike_data( const index idx ) const;

  void clear_sources();

  void finalize();

  size_t num_unique_sources( const thread tid ) const;

  void reset_entry_points( const thread num_ranks, const bool compressed );

  //! Return the current compressed target and increment the current compressed target index
  std::pair< index, size_t > get_next_compressed_target( const thread source_rank );

  //! Return the current target and increment the current target index
  std::tuple< std::pair< index, size_t >, thread > get_next_target( const thread source_rank );

  //! Return if the last target for this source rank has been processed already
  bool reached_last_target( const thread source_rank, const bool compressed ) const;

  //! Make sure the current compressed target index points to a valid target
  bool seek_next_compressed_target( const thread source_rank );

  //! Make sure the current target index points to a valid target
  bool seek_next_target( const thread source_rank );

  void prepare_compressed_targets( const thread num_threads, const thread num_ranks );

  void clear_compressed_indices();
};

inline void
AdjacencyList::resize( const thread num_threads, const thread num_ranks, const bool compressed )
{
  adjacency_list_.resize( num_threads );
  sources_.resize( num_threads, std::vector< std::map< index, size_t > >( num_ranks ) );
  if ( compressed )
  {
    source_to_compressed_index_.resize( num_ranks );
  }
  else
  {
    clear_compressed_indices();
  }
}

inline std::pair< std::vector< AdjacencyListTarget >::const_iterator,
  std::vector< AdjacencyListTarget >::const_iterator >
AdjacencyList::get_iterators( const thread tid, const index adjacency_list_index ) const
{
  assert( tid >= 0 );
  assert( static_cast< size_t >( tid ) < adjacency_list_.size() );
  assert( adjacency_list_index < adjacency_list_[ tid ].size() );

  return { adjacency_list_[ tid ][ adjacency_list_index ].cbegin(),
    adjacency_list_[ tid ][ adjacency_list_index ].cend() };
}

inline void
AdjacencyList::clear_sources()
{
  std::vector< std::vector< std::map< index, index > > >().swap( sources_ );
}

inline void
AdjacencyList::clear_compressed_indices()
{
  std::vector< std::map< index, size_t > >().swap( source_to_compressed_index_ );
}

inline void
AdjacencyList::finalize()
{
  std::vector< std::vector< std::vector< AdjacencyListTarget > > >().swap( adjacency_list_ );
  std::vector< std::vector< std::map< index, index > > >().swap( sources_ );
  std::vector< std::map< thread, index > >().swap( compressed_indices_ );
  std::vector< std::map< index, size_t > >().swap( source_to_compressed_index_ );
  std::vector< std::map< index, size_t >::const_iterator >().swap( next_compressed_index_ );
  std::vector< std::map< index, index >::const_iterator >().swap( next_source_index_ );
  std::vector< size_t >().swap( next_source_index_thread_ );
}

inline size_t
AdjacencyList::num_unique_sources( const thread tid ) const
{
  return adjacency_list_[ tid ].size();
}

inline const std::map< thread, index >&
AdjacencyList::get_compressed_spike_data( const index idx ) const
{
  return compressed_indices_[ idx ];
}

inline void
AdjacencyList::reset_entry_points( const thread num_ranks, const bool compressed )
{
  if ( compressed )
  {
    next_compressed_index_.clear();
    next_compressed_index_.reserve( num_ranks );
    for ( thread source_rank = 0; source_rank != num_ranks; ++source_rank )
    {
      next_compressed_index_.push_back( source_to_compressed_index_[ source_rank ].cbegin() );
    }
  }
  else
  {
    next_source_index_thread_.assign( num_ranks, 0 );
    next_source_index_.clear();
    next_source_index_.reserve( num_ranks );
    for ( thread source_rank = 0; source_rank != num_ranks; ++source_rank )
    {
      next_source_index_.push_back( sources_[ 0 ][ source_rank ].cbegin() );
    }
  }
}

inline bool
AdjacencyList::seek_next_compressed_target( const thread source_rank )
{
  // check if the last target has been reached already
  if ( reached_last_target( source_rank, true ) )
  {
    return false;
  }
  return true;
}

inline std::pair< index, size_t >
AdjacencyList::get_next_compressed_target( const thread source_rank )
{
  // return next pair of source and compressed index and afterward increase index
  return *( next_compressed_index_[ source_rank ]++ );
}

inline bool
AdjacencyList::reached_last_target( const thread source_rank, const bool compressed ) const
{
  if ( compressed )
  {
    return next_compressed_index_[ source_rank ] == source_to_compressed_index_[ source_rank ].cend();
  }
  else
  {
    return next_source_index_thread_[ source_rank ] == sources_.size();
  }
}

inline bool
AdjacencyList::seek_next_target( const thread source_rank )
{
  while (
    next_source_index_[ source_rank ] == sources_[ next_source_index_thread_[ source_rank ] ][ source_rank ].cend() )
  {
    ++next_source_index_thread_[ source_rank ];
    // check if the last target has been reached already (i.e., the next thread id equals the total number of threads)
    if ( reached_last_target( source_rank, false ) )
    {
      return false;
    }
    next_source_index_[ source_rank ] = sources_[ next_source_index_thread_[ source_rank ] ][ source_rank ].cbegin();
  }
  return true;
}

inline std::tuple< std::pair< index, size_t >, thread >
AdjacencyList::get_next_target( const thread source_rank )
{
  std::pair< index, index > next_target = *( next_source_index_[ source_rank ] );
  thread next_target_thread = next_source_index_thread_[ source_rank ];

  ++next_source_index_[ source_rank ];

  return { next_target, next_target_thread };
}

inline void
AdjacencyList::prepare_compressed_targets( const thread num_threads, const thread num_ranks )
{
  index source_node_id;
  size_t adjacency_list_index;
  for ( thread tid = 0; tid != num_threads; ++tid )
  {
    for ( thread source_rank = 0; source_rank != num_ranks; ++source_rank )
    {
      for ( auto& source : sources_[ tid ][ source_rank ] )
      {
        std::tie( source_node_id, adjacency_list_index ) = source;
        // first check if there has already been another thread with a connection from this source node
        const auto it = source_to_compressed_index_[ source_rank ].find( source_node_id );

        size_t source_node_idx;
        if ( it
          == source_to_compressed_index_[ source_rank ].end() ) // this source node id is not yet known to any thread
        {
          source_to_compressed_index_[ source_rank ][ source_node_id ] = compressed_indices_.size();
          source_node_idx = compressed_indices_.size();
          // TODO JV: Combine the following two lines to a single one
          compressed_indices_.emplace_back(); // append empty map
          compressed_indices_[ source_node_idx ].emplace( tid, adjacency_list_index );
        }
        else
        {
          compressed_indices_[ it->second ][ tid ] = adjacency_list_index;
        }
      }
    }
  }
}

inline void
AdjacencyList::add_target( const thread tid,
  const synindex syn_id,
  const index source_node_id,
  const thread source_rank,
  const index local_target_node_id,
  const index local_target_connection_id,
  const delay axonal_delay )
{
  assert( tid >= 0 );
  assert( static_cast< size_t >( tid ) < adjacency_list_.size() );

  assert( static_cast< size_t >( tid ) < sources_.size() );
  assert( static_cast< size_t >( source_rank ) < sources_[ tid ].size() );

  // without compression, the sources data structure has to be filled to map to adjacency list positions
  auto source_index = sources_[ tid ][ source_rank ].find( source_node_id );
  // Check if this is the first connection from this source node to any target node managed by this thread
  if ( source_index != sources_[ tid ][ source_rank ].end() ) // not the first connection
  {
    adjacency_list_[ tid ][ source_index->second ].emplace_back(
      local_target_node_id, local_target_connection_id, syn_id, axonal_delay );
  }
  else // actually the first connection
  {
    const index new_index = adjacency_list_[ tid ].size(); // set index for this source node id
    sources_[ tid ][ source_rank ][ source_node_id ] = new_index;
    adjacency_list_[ tid ].emplace_back( std::initializer_list< AdjacencyListTarget > {
      { local_target_node_id, local_target_connection_id, syn_id, axonal_delay } } );
  }
}

} // nest

#endif // USE_ADJACENCY_LIST
#endif // ADJACENCY_LIST_H
