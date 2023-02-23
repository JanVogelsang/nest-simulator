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

// TODO JV: Define only for debugging!
// #define USE_ADJACENCY_LIST ON
// TODO JV: Test with very low communication sizes, so multiple communication rounds are needed (both spikes and
//  targets)

#ifndef NEST_ADJACENCY_LIST_H
#define NEST_ADJACENCY_LIST_H
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
  delay partial_delay : NUM_BITS_DELAY;
  // TODO JV (pt): Still some bits to spare here

  AdjacencyListTarget( const index local_target_node_id,
    const index local_target_connection_id,
    const synindex syn_id,
    const delay partial_delay )
    : local_target_node_id( local_target_node_id )
    , local_target_connection_id( local_target_connection_id )
    , syn_id( syn_id )
    , partial_delay( partial_delay )
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
  std::vector< std::vector< std::vector< AdjacencyListTarget > > > adjacency_list_;

  /**
   * Intermediate structure to map source neurons to the corresponding index in the adjacency_list for each thread.
   * Deleted after communication of targets to pre-synaptic processes.
   * Two dimensional object:
   *   - first dim: threads
   *   - second dim: map of source node id to index in adjacency_list
   */
  std::vector< std::map< index, size_t > > sources_;

  /**
   * A structure to map incoming spikes to the adjacency list entry for each thread on the postsynaptic side if spike
   * compression is enabled. Internally arranged in a 2d structure:
   * sources|(thread->adjacency list index)
   */
  std::vector< std::map< thread, index > > compressed_indices_;

  /**
   * Intermediate structure. For each source node id, the "compressed" index in the structure which stores all adjacency
   * list indices for each source node, one for each thread. See compressed_indices_.
   */
  std::map< index, size_t > source_to_compressed_index_;

  /**
   * For each reading thread, save the position of the next pair that is going to be retrieved. Required to communicate
   * target data for compressed spikes.
   */
  std::vector< std::map< index, size_t >::const_iterator > next_compressed_index_;

  /**
   * For each reading thread, save the position of the next pair that is going to be retrieved. Required to communicate
   * target data.
   */
  std::vector< std::map< index, index >::const_iterator > next_source_index_;
  /**
   * For each reading thread, save the thread that owns the next pair that is going to be retrieved. Required to
   * communicate target data.
   */
  std::vector< size_t > next_source_index_thread_;

public:
  void
  resize( const thread num_threads )
  {
    sources_.resize( num_threads );
    adjacency_list_.resize( num_threads );
  }

  void add_target( const thread tid,
    const synindex syn_id,
    const index source_node_id,
    const index target_node_id,
    const index target_connection_id,
    const delay partial_delay,
    const bool prepare_for_compression );

  std::pair< std::vector< AdjacencyListTarget >::const_iterator, std::vector< AdjacencyListTarget >::const_iterator >
  get_iterators( const thread tid, const index adjacency_list_index ) const;

  const std::map< thread, index >& get_compressed_spike_data( const index idx ) const;

  void clear_sources( const thread tid );

  void finalize();

  size_t num_unique_sources( const thread tid ) const;

  void reset_entry_point( const thread num_threads );

  void reject_last_target_data( const thread tid );

  void no_targets_to_process( const thread tid );

  std::tuple< std::pair< index, size_t >, bool > get_next_compressed_target( const thread tid );

  std::tuple< std::pair< index, size_t >, thread, bool > get_next_target( const thread tid );

  void clear_compressed_indices();
};

inline std::pair< std::vector< AdjacencyListTarget >::const_iterator,
  std::vector< AdjacencyListTarget >::const_iterator >
AdjacencyList::get_iterators( const thread tid, const index adjacency_list_index ) const
{
  assert( tid >= 0 );
  assert( static_cast<size_t>(tid) < adjacency_list_.size() );
  assert( adjacency_list_index < adjacency_list_[ tid ].size() );

  return { adjacency_list_[ tid ][ adjacency_list_index ].cbegin(),
    adjacency_list_[ tid ][ adjacency_list_index ].cend() };
}

inline void
AdjacencyList::add_target( const thread tid,
  const synindex syn_id,
  const index source_node_id,
  const index local_target_node_id,
  const index local_target_connection_id,
  const delay partial_delay,
  const bool prepare_for_compression )
{
  assert( tid >= 0 );
  assert( static_cast<size_t>(tid) < adjacency_list_.size() );
  assert( static_cast<size_t>(tid) < sources_.size() );

  auto source_index = sources_[ tid ].find( source_node_id );

  // Check if this is the first connection from this source node to any target node managed by this thread
  if ( source_index != sources_[ tid ].end() ) // not the first connection
  {
    adjacency_list_[ tid ][ ( *source_index ).second ].emplace_back(
      local_target_node_id, local_target_connection_id, syn_id, partial_delay );
  }
  else  // actually the first connection
  {
    const index new_index = adjacency_list_[ tid ].size(); // set index for this source node id
    sources_[ tid ][ source_node_id ] = new_index;
    adjacency_list_[ tid ].emplace_back( std::initializer_list< AdjacencyListTarget > {
      { local_target_node_id, local_target_connection_id, syn_id, partial_delay } } );

    // if spike compression is enabled, fill the compression data structures as well
    if ( prepare_for_compression )
    {
      // first check if there has already been another thread with a connection from this source node
      const auto it = source_to_compressed_index_.find( source_node_id );

      if ( it == source_to_compressed_index_.end() ) // this source node id is not yet known to any thread
      {
        source_to_compressed_index_[ source_node_id ] = compressed_indices_.size();
        compressed_indices_.emplace_back(); // append empty map
      }
      // add the index to the adjacency list for this thread for the new source
      compressed_indices_[ it->second ][ tid ] = new_index;
    }
  }
}

inline void
AdjacencyList::clear_sources( const thread tid )
{
  std::map< index, index >().swap( sources_[ tid ] );
}

inline void
AdjacencyList::clear_compressed_indices()
{
  std::map< index, size_t >().swap( source_to_compressed_index_ );
}

inline void
AdjacencyList::finalize()
{
  std::vector< std::vector< std::vector< AdjacencyListTarget > > >().swap( adjacency_list_ );
  std::vector< std::map< index, index > >().swap( sources_ );
  std::vector< std::map< thread, index > >().swap( compressed_indices_ );
  std::map< index, size_t >().swap( source_to_compressed_index_ );
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
AdjacencyList::reset_entry_point( const thread num_threads )
{
  next_compressed_index_.assign( num_threads, source_to_compressed_index_.cbegin() );
  next_source_index_.assign( num_threads, sources_[ 0 ].cbegin() );
  next_source_index_thread_.assign( num_threads, 0 );
}

inline void
AdjacencyList::reject_last_target_data( const thread tid )
{
  // just update both indices, as this is less costly than another branch and won't hurt
  --next_compressed_index_[ tid ];
  --next_source_index_[ tid ];
}

inline void
AdjacencyList::no_targets_to_process( const thread tid )
{
  next_compressed_index_[ tid ] = source_to_compressed_index_.cend();
  next_source_index_thread_[ tid ] = sources_.size() - 1;
  next_source_index_[ tid ] = sources_[ next_source_index_thread_[ tid ] ].cbegin();
}

inline std::tuple< std::pair< index, size_t >, bool >
AdjacencyList::get_next_compressed_target( const thread tid )
{
  // check if the last target has been reached already
  if ( next_compressed_index_[ tid ] == source_to_compressed_index_.cend() )
  {
    return { { 0, 0 }, false };
  }

  // return next pair of source and compressed index and afterwards increase index
  return { *( next_compressed_index_[ tid ]++ ), true };
}

inline std::tuple< std::pair< index, size_t >, thread, bool >
AdjacencyList::get_next_target( const thread tid )
{
  // ensure valid iterator
  while ( next_source_index_[ tid ] == sources_[ next_source_index_thread_[ tid ] ].cend() )
  {
    ++next_source_index_thread_[ tid ];

    // check if the last target has been reached already
    if ( next_source_index_thread_[ tid ] == sources_.size() )
    {
      return { { 0, 0 }, -1, false };
    }

    next_source_index_[ tid ] = sources_[ next_source_index_thread_[ tid ] ].cbegin();
  }

  std::pair< index, index > next_target = *( next_source_index_[ tid ] );
  thread next_target_thread = next_source_index_thread_[ tid ];

  ++next_source_index_[ tid ];
  return { next_target, next_target_thread, true };
}

} // nest

#endif // USE_ADJACENCY_LIST
#endif // NEST_ADJACENCY_LIST_H
