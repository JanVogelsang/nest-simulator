/*
 *  source_manager.cpp
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

#include "source_manager.h"

// Includes from nestkernel:
#include "kernel_manager.h"
#include "vp_manager_impl.h"

namespace nest
{

SourceManager::SourceManager()
{
}

SourceManager::~SourceManager()
{
}

void
SourceManager::set_status( const DictionaryDatum& )
{
}

void
SourceManager::get_status( DictionaryDatum& )
{
}

void
SourceManager::change_number_of_threads()
{
  finalize();
  initialize();
}

void
SourceManager::initialize()
{
  is_cleared_.initialize( kernel().vp_manager.get_num_threads(), false );
  current_positions_.resize( kernel().mpi_manager.get_num_processes() );
}

void
SourceManager::finalize()
{
  // only clear sources when there were nodes added to the simulation already
  if ( kernel().node_manager.size() > 0 )
  {
#pragma omp parallel
    {
      const thread tid = kernel().vp_manager.get_thread_id();
      if ( is_cleared_[ tid ].is_false() )
      {
        clear( tid );
      }
    }
  }

  current_positions_.clear();
}

bool
SourceManager::is_cleared() const
{
  return is_cleared_.all_true();
}

SourceTablePosition
SourceManager::find_maximal_position() const
{
  SourceTablePosition max_position( -1, -1, -1, -1 );
  for ( thread source_rank = 0; source_rank < kernel().mpi_manager.get_num_processes(); ++source_rank )
  {
    if ( max_position < current_positions_[ source_rank ] )
    {
      max_position = current_positions_[ source_rank ];
    }
  }
  return max_position;
}

void
SourceManager::clean( const thread tid )
{
  // has_source_.clear();

  // Find maximal position in source table among threads to make sure
  // unprocessed entries are not removed. Given this maximal position,
  // we can safely delete all larger entries since they will not be
  // touched anymore.
  const SourceTablePosition max_position = find_maximal_position();
  const SparseNodeArray& thread_local_nodes = kernel().node_manager.get_local_nodes( tid );

  // If this thread corresponds to max_position's thread, we can only
  // delete part of the sources table, with indices larger than those
  // in max_position; if this thread is larger than max_positions's
  // thread, we can delete all sources; otherwise we do nothing.
  if ( max_position.tid == tid )
  {
    // clear all sources of all visited nodes
    for ( SparseNodeArray::const_iterator n = thread_local_nodes.begin() + max_position.local_target_node_id + 1;
          n != thread_local_nodes.end();
          ++n )
    {
      n->get_node()->clear_sources();
    }
    // don't remove any sources of the current node yet, they can be removed at a later point in time
  }
  else if ( max_position.tid < tid )
  {
    for ( SparseNodeArray::const_iterator n = thread_local_nodes.begin(); n != thread_local_nodes.end(); ++n )
    {
      n->get_node()->clear_sources();
    }
  }
  else
  {
    // do nothing
    assert( tid < max_position.tid );
  }
}

index
SourceManager::get_node_id( const thread tid,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_connection_id ) const
{
  if ( not kernel().connection_manager.get_keep_source_table() )
  {
    throw KernelException( "Cannot use SourceManager::get_node_id when get_keep_source_table is false" );
  }

  return kernel()
    .node_manager.thread_lid_to_node( tid, local_target_node_id )
    ->get_source( syn_id, local_connection_id );
}

#ifndef USE_ADJACENCY_LIST
void
SourceManager::reset_entry_points()
{
  // Since we read the source table backwards, we need to set saved values to the biggest possible value. These will be
  // used to initialize current_positions_ correctly upon calling restore_entry_point. However, this can only be done if
  // other values have valid values.

  const thread num_threads = kernel().vp_manager.get_num_threads();
  const thread num_ranks = kernel().mpi_manager.get_num_processes();
  if ( num_threads > 0 )
  {
    // The saved position is explicitly set to an invalid position with the thread id being set one value too high, so
    // the decrease call looks for the next valid position
    for ( thread source_rank = 0; source_rank != num_ranks; ++source_rank )
    {
      current_positions_[ source_rank ].tid = num_threads;
      current_positions_[ source_rank ].syn_id = -1;
      current_positions_[ source_rank ].local_target_node_id = -1;
      current_positions_[ source_rank ].local_target_connection_id = -1;
      current_positions_[ source_rank ].decrease( source_rank );
    }
  }
  else
  {
    for ( thread source_rank = 0; source_rank != num_ranks; ++source_rank )
    {
      current_positions_[ source_rank ].tid = -1;
      current_positions_[ source_rank ].syn_id = -1;
      current_positions_[ source_rank ].local_target_node_id = -1;
      current_positions_[ source_rank ].local_target_connection_id = -1;
    }
  }
}

bool
SourceManager::is_next_target_valid( const thread source_rank )
{
  return not current_positions_[ source_rank ].reached_end();
}

void
SourceManager::get_next_target_data( const thread source_rank, TargetData& next_target_data )
{
  SourceTablePosition& current_position = current_positions_[ source_rank ];
  index source_node_id = current_position.current_source;

  // if ( primary ) {
  next_target_data.set_is_primary( true );
  next_target_data.set_source_lid( kernel().vp_manager.node_id_to_lid( source_node_id ) );
  // we store the thread index of the target neuron, not our own tid!
  next_target_data.set_source_tid(
    kernel().vp_manager.vp_to_thread( kernel().vp_manager.node_id_to_vp( source_node_id ) ) );
  next_target_data.reset_marker();
  next_target_data.target_data.set_tid( current_position.tid );
  next_target_data.target_data.set_syn_id( current_position.syn_id );
  next_target_data.target_data.set_local_target_node_id( current_position.local_target_node_id );
  next_target_data.target_data.set_local_target_connection_id( current_position.local_target_connection_id );
  // TODO JV (pt): Secondary events. Secondary event can't be combined yet with the adjacency-list-based
  //  delivery which offers quite straightforward compression. In general, the adjacency list is not built
  //  for secondary events yet, so one would have to rethink the whole handling of secondary events here. In theory,
  //  there should not even be a strict difference between secondary and primary events, as they can be delivered the
  //  exact same way (afaik). Simply using the adjacency-list-based delivery with all sorts of events (from node to
  //  node) should actually work.
  // } else {  // -> secondary
  /*next_target_data.set_is_primary( false );

  // the source rank will write to the buffer position relative to the first position from the absolute position in
  // the receive buffer
  const size_t relative_recv_buffer_pos = kernel().connection_manager.get_secondary_recv_buffer_position(
                                            current_position.tid, current_position.syn_id, current_position.lcid )
    - kernel().mpi_manager.get_recv_displacement_secondary_events_in_int( source_rank );

  SecondaryTargetDataFields& secondary_fields = next_target_data.secondary_data;
  secondary_fields.set_recv_buffer_pos( relative_recv_buffer_pos );
  secondary_fields.set_syn_id( current_position.syn_id );
   }*/
  current_position.decrease( source_rank );
}
#endif

void
SourceManager::clear( const thread tid )
{
  // TODO JV: Make sure iteration over all nodes is efficient
  const SparseNodeArray& thread_local_nodes = kernel().node_manager.get_local_nodes( tid );

  for ( SparseNodeArray::const_iterator n = thread_local_nodes.begin(); n != thread_local_nodes.end(); ++n )
  {
    n->get_node()->clear_sources();
  }

  is_cleared_[ tid ].set_true();
}

index
SourceManager::remove_disabled_sources( const thread tid, const synindex syn_id )
{
  assert( false ); // TODO JV (pt): Structural plasticity

  /*if ( sources_[ tid ].size() <= syn_id )
  {
    return invalid_index;
  }

  std::vector< Source >& mysources = sources_[ tid ][ syn_id ];
  const index max_size = mysources.size();
  if ( max_size == 0 )
  {
    return invalid_index;
  }

  // lcid needs to be signed, to allow lcid >= 0 check in while loop
  // to fail; afterwards we can be certain that it is non-negative and
  // we can static_cast it to index
  long lcid = max_size - 1;
  while ( lcid >= 0 and mysources[ lcid ].is_disabled() )
  {
    --lcid;
  }
  ++lcid; // lcid marks first disabled source, but the while loop only
          // exits if lcid points at a not disabled element, hence we
          // need to increase it by one again
  mysources.erase( mysources.begin() + lcid, mysources.end() );
  if ( static_cast< index >( lcid ) == max_size )
  {
    return invalid_index;
  }
  return static_cast< index >( lcid );*/
}

void
SourceManager::compute_buffer_pos_for_unique_secondary_sources( const thread tid,
  std::map< index, size_t >& buffer_pos_of_source_node_id_syn_id )
{
  assert( false ); // TODO JV (pt): Secondary events

  /*
  // set of unique sources & synapse types, required to determine
  // secondary events MPI buffer positions
  // initialized and deleted by thread 0 in this method
  static std::set< std::pair< index, size_t > >* unique_secondary_source_node_id_syn_id;
#pragma omp single
  {
    unique_secondary_source_node_id_syn_id = new std::set< std::pair< index, size_t > >();
  }

  // collect all unique pairs of source node ID and synapse-type id
  // corresponding to continuous-data connections on this MPI rank;
  // using a set makes sure secondary events are not duplicated for
  // targets on the same process, but different threads
  for ( size_t syn_id = 0; syn_id < sources_[ tid ].size(); ++syn_id )
  {
    if ( not kernel().model_manager.get_connection_model( syn_id, tid ).is_primary() )
    {
      for ( std::vector< Source >::const_iterator source_cit = sources_[ tid ][ syn_id ].begin();
            source_cit != sources_[ tid ][ syn_id ].end();
            ++source_cit )
      {
#pragma omp critical
        {
          ( *unique_secondary_source_node_id_syn_id ).insert( std::make_pair( source_cit->get_node_id(), syn_id ) );
        }
      }
    }
  }
#pragma omp barrier

#pragma omp single
  {
    // compute receive buffer positions for all unique pairs of source
    // node ID and synapse-type id on this MPI rank
    std::vector< int > recv_counts_secondary_events_in_int_per_rank( kernel().mpi_manager.get_num_processes(), 0 );

    for ( std::set< std::pair< index, size_t > >::const_iterator cit =
            ( *unique_secondary_source_node_id_syn_id ).begin();
          cit != ( *unique_secondary_source_node_id_syn_id ).end();
          ++cit )
    {
      const thread source_rank = kernel().mpi_manager.get_process_id_of_node_id( cit->first );
      const size_t event_size = kernel().model_manager.get_secondary_event_prototype( cit->second, tid ).size();

      buffer_pos_of_source_node_id_syn_id.insert(
        std::make_pair( pack_source_node_id_and_syn_id( cit->first, cit->second ),
          recv_counts_secondary_events_in_int_per_rank[ source_rank ] ) );

      recv_counts_secondary_events_in_int_per_rank[ source_rank ] += event_size;
    }

    // each chunk needs to contain one additional int that can be used
    // to communicate whether waveform relaxation has converged
    for ( auto& recv_count : recv_counts_secondary_events_in_int_per_rank )
    {
      ++recv_count;
    }

    kernel().mpi_manager.set_recv_counts_secondary_events_in_int_per_rank(
      recv_counts_secondary_events_in_int_per_rank );
    delete unique_secondary_source_node_id_syn_id;
  } // of omp single
   */
}

}
