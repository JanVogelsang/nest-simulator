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
#include "mpi_manager_impl.h"
#include "target_data.h"
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
  assert( sizeof( Source ) == 8 );
  const thread num_threads = kernel().vp_manager.get_num_threads();
  is_cleared_.initialize( num_threads, false );
  current_positions_.resize( num_threads );
  saved_positions_.resize( num_threads );
  compressible_sources_.resize( num_threads );
  compressed_spike_data_map_.resize( num_threads );
  has_source_.resize( num_threads );

#pragma omp parallel
  {
    const thread tid = kernel().vp_manager.get_thread_id();
    compressible_sources_[ tid ].resize( 0 );
    compressed_spike_data_map_[ tid ].resize( 0 );
  } // of omp parallel
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
        compressible_sources_[ tid ].clear();
        compressed_spike_data_map_[ tid ].clear();
      }
    }
  }

  current_positions_.clear();
  saved_positions_.clear();
  compressible_sources_.clear();
  compressed_spike_data_map_.clear();
  has_source_.clear();
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
  for ( thread tid = 0; tid < kernel().vp_manager.get_num_threads(); ++tid )
  {
    if ( max_position < saved_positions_[ tid ] )
    {
      max_position = saved_positions_[ tid ];
    }
  }
  return max_position;
}

void
SourceManager::clean( const thread tid )
{
  has_source_.clear();

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

    /*// for the node for which sources are currently transferred, only remove part of the sources yet
    for ( synindex syn_id = max_position.syn_id; syn_id < kernel().model_manager.get_num_connection_models();
          ++syn_id )
    {
      const Node *current_node = kernel().node_manager.get_local_nodes( tid ).get_node_by_index(
    max_position.local_target_node_id ); if ( max_position.syn_id == syn_id )
      {
        // we need to add 2 to max_position.lcid since
        // max_position.lcid + 1 can contain a valid entry which we
        // do not want to delete.
        if ( max_position.local_target_connection_id + 2 < static_cast< long >( sources.size() ) )
        {
          current_node->erase_sources( max_position.local_target_connection_id + 2 );
        }
      }
      else
      {
        assert( max_position.syn_id < syn_id );
        current_node->clear_sources( syn_id );
      }
    }*/
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
    ->get_source( syn_id, local_connection_id )
    .get_node_id();
}


void
SourceManager::reset_entry_point( const thread tid )
{
  // Since we read the source table backwards, we need to set saved values to the biggest possible value. These will be
  // used to initialize current_positions_ correctly upon calling restore_entry_point. However, this can only be done if
  // other values have valid values.

  const thread num_threads = kernel().vp_manager.get_num_threads();
  if ( num_threads > 0 )
  {
    // The saved position is explicitly set to an invalid position with the thread id being set one value too high, so
    // the decrease call looks for the next valid position
    saved_positions_[ tid ].tid = num_threads;
    saved_positions_[ tid ].syn_id = -1;
    saved_positions_[ tid ].local_target_node_id = -1;
    saved_positions_[ tid ].local_target_connection_id = -1;
    saved_positions_[ tid ].decrease();
  }
  else
  {
    saved_positions_[ tid ].tid = -1;
    saved_positions_[ tid ].syn_id = -1;
    saved_positions_[ tid ].local_target_node_id = -1;
    saved_positions_[ tid ].local_target_connection_id = -1;
  }
}

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

void
SourceManager::reject_last_target_data( const thread tid )
{
  // The last target data returned by get_next_target_data() could not be inserted into MPI buffer due to overflow.
  // We hence need to correct the processed flag of the last entry
  current_positions_[ tid ].increase();
  kernel()
    .node_manager.get_local_nodes( current_positions_[ tid ].tid )
    .get_node_by_index( current_positions_[ tid ].local_target_node_id )
    ->get_source( current_positions_[ tid ].syn_id, current_positions_[ tid ].local_target_connection_id )
    .set_processed( false );
}

void
SourceManager::reset_processed_flags( const thread tid )
{
  // TODO JV (help): Make sure iteration over all nodes is efficient
  const SparseNodeArray& thread_local_nodes = kernel().node_manager.get_local_nodes( tid );

  for ( SparseNodeArray::const_iterator n = thread_local_nodes.begin(); n != thread_local_nodes.end(); ++n )
  {
    n->get_node()->reset_sources_processed_flags();
  }
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

bool
SourceManager::source_should_be_processed_( const thread rank_start, const thread rank_end, const Source& source ) const
{
  const thread source_rank = kernel().mpi_manager.get_process_id_of_node_id( source.get_node_id() );

  return not( source.is_processed()
    or source.is_disabled()
    // is this thread responsible for this part of the MPI buffer?
    or source_rank < rank_start or rank_end <= source_rank );
}

bool
SourceManager::populate_target_data_fields_( const SourceTablePosition& current_position,
  const Source& current_source,
  const thread source_rank,
  TargetData& next_target_data ) const
{
  const auto node_id = current_source.get_node_id();

  // set values of next_target_data
  next_target_data.set_source_lid( kernel().vp_manager.node_id_to_lid( node_id ) );
  next_target_data.set_source_tid( kernel().vp_manager.vp_to_thread( kernel().vp_manager.node_id_to_vp( node_id ) ) );
  next_target_data.reset_marker();

  if ( current_source.is_primary() ) // primary connection, i.e., chemical synapses
  {
    next_target_data.set_is_primary( true );

    TargetDataFields& target_fields = next_target_data.target_data;
    target_fields.set_syn_id( current_position.syn_id );
    if ( kernel().connection_manager.use_compressed_spikes() )
    {
      // WARNING: we set the tid field here to zero just to make sure
      // it has a defined value; however, this value is _not_ used
      // anywhere when using compressed spikes
      target_fields.set_tid( 0 );
      auto it_idx = compressed_spike_data_map_.at( current_position.tid )
                      .at( current_position.syn_id )
                      .find( current_source.get_node_id() );
      if ( it_idx != compressed_spike_data_map_.at( current_position.tid ).at( current_position.syn_id ).end() )
      {
        // WARNING: no matter how tempting, do not try to remove this
        // entry from the compressed_spike_data_map_; if the MPI buffer
        // is already full, this entry will need to be communicated the
        // next MPI comm round, which, naturally, is not possible if it
        // has been removed
        // TODO JV: Spike compression
        // target_fields.set_lcid( it_idx->second );
      }
      else // another thread is responsible for communicating this compressed source
      {
        return false;
      }
    }
    else
    {
      // we store the thread index of the source table, not our own tid!
      target_fields.set_tid( current_position.tid );
      target_fields.set_local_target_node_id( current_position.local_target_node_id );
      target_fields.set_local_target_connection_id( current_position.local_target_connection_id );
    }
  }
  else // secondary connection, e.g., gap junctions
  {
    assert( false ); // TODO JV (pt): Secondary events
    /*next_target_data.set_is_primary( false );

    // the source rank will write to the buffer position relative to the first position from the absolute position in
    // the receive buffer
    const size_t relative_recv_buffer_pos = kernel().connection_manager.get_secondary_recv_buffer_position(
                                              current_position.tid, current_position.syn_id, current_position.lcid )
      - kernel().mpi_manager.get_recv_displacement_secondary_events_in_int( source_rank );

    SecondaryTargetDataFields& secondary_fields = next_target_data.secondary_data;
    secondary_fields.set_recv_buffer_pos( relative_recv_buffer_pos );
    secondary_fields.set_syn_id( current_position.syn_id );*/
  }

  return true;
}

bool
SourceManager::get_next_target_data( const thread tid,
  const thread rank_start,
  const thread rank_end,
  thread& source_rank,
  TargetData& next_target_data )
{
  SourceTablePosition& current_position = current_positions_[ tid ];

  // we stay in this loop either until we can return a valid TargetData object or we have reached the end of all
  // sources tables
  while ( not current_position.reached_end() )
  {

    // the current position contains an entry, so we retrieve it
    Source& current_source = kernel()
                               .node_manager.get_local_nodes( current_position.tid )
                               .get_node_by_index( current_position.local_target_node_id )
                               ->get_source( current_position.syn_id, current_position.local_target_connection_id );

    if ( not source_should_be_processed_( rank_start, rank_end, current_source ) )
    {
      current_position.decrease();
      continue;
    }

    // reaching this means we found an entry that should be communicated via MPI, so we prepare to return the relevant
    // data
    // set the source rank
    source_rank = kernel().mpi_manager.get_process_id_of_node_id( current_source.get_node_id() );

    if ( not populate_target_data_fields_( current_position, current_source, source_rank, next_target_data ) )
    {
      current_position.decrease();
      continue;
    }

    // we are about to return a valid entry, so mark it as processed
    current_source.set_processed( true );

    current_position.decrease();
    return true; // found a valid entry
  }
  return false; // reached the end of all sources tables
}

void
SourceManager::resize_compressible_sources()
{
  for ( thread tid = 0; tid < static_cast< thread >( compressible_sources_.size() ); ++tid )
  {
    compressible_sources_[ tid ].clear();
    compressible_sources_[ tid ].resize(
      kernel().model_manager.get_num_connection_models(), std::map< index, SpikeData >() );
  }
}

void
SourceManager::collect_compressible_sources( const thread tid )
{
  assert( false ); // TODO JV: Spike compression
  /*
  for ( synindex syn_id = 0; syn_id < sources_[ tid ].size(); ++syn_id )
  {
    index local_target_connection_id = 0;
    auto& syn_sources = sources_[ tid ][ syn_id ];

    for ( auto node_id )
      while ( lcid < syn_sources.size() )
      {
        const index old_source_node_id = syn_sources[ lcid ].get_node_id();
        const std::pair< index, SpikeData > source_node_id_to_spike_data =
          std::make_pair( old_source_node_id, SpikeData( tid, syn_id, lcid, 0 ) );
        compressible_sources_[ tid ][ syn_id ].insert( source_node_id_to_spike_data );

        // find next source with different node_id (assumes sorted sources)
        ++lcid;
        while ( ( lcid < syn_sources.size() ) and ( syn_sources[ lcid ].get_node_id() == old_source_node_id ) )
        {
          ++lcid;
        }
      }
  }*/
}

void
SourceManager::fill_compressed_spike_data(
  std::vector< std::vector< std::vector< SpikeData > > >& compressed_spike_data )
{
  assert( false ); // TODO JV: Spike compression
  /*
  compressed_spike_data.clear();
  compressed_spike_data.resize( kernel().model_manager.get_num_connection_models() );

  for ( thread tid = 0; tid < static_cast< thread >( compressible_sources_.size() ); ++tid )
  {
    compressed_spike_data_map_[ tid ].clear();
    compressed_spike_data_map_[ tid ].resize(
      kernel().model_manager.get_num_connection_models(), std::map< index, size_t >() );
  }

  // pseudo-random thread selector to balance memory usage across
  // threads of compressed_spike_data_map_
  size_t thread_idx = 0;

  // for each local thread and each synapse type we will populate this
  // vector with spike data containing information about all process
  // local targets
  std::vector< SpikeData > spike_data;

  for ( thread tid = 0; tid < static_cast< thread >( compressible_sources_.size() ); ++tid )
  {
    for ( synindex syn_id = 0; syn_id < compressible_sources_[ tid ].size(); ++syn_id )
    {
      for ( auto it = compressible_sources_[ tid ][ syn_id ].begin();
            it != compressible_sources_[ tid ][ syn_id ].end(); )
      {
        spike_data.clear();

        // add target position on this thread
        spike_data.push_back( it->second );

        // add target positions on all other threads
        for ( thread other_tid = tid + 1; other_tid < static_cast< thread >( compressible_sources_.size() );
              ++other_tid )
        {
          auto other_it = compressible_sources_[ other_tid ][ syn_id ].find( it->first );
          if ( other_it != compressible_sources_[ other_tid ][ syn_id ].end() )
          {
            spike_data.push_back( other_it->second );
            compressible_sources_[ other_tid ][ syn_id ].erase( other_it );
          }
        }

        // WARNING: store source-node-id -> process-global-synapse
        // association in compressed_spike_data_map on a
        // pseudo-randomly selected thread which houses targets for
        // this source; this tries to balance memory usage of this
        // data structure across threads
        const thread responsible_tid = spike_data[ thread_idx % spike_data.size() ].get_tid();
        ++thread_idx;

        compressed_spike_data_map_[ responsible_tid ][ syn_id ].insert(
          std::make_pair( it->first, compressed_spike_data[ syn_id ].size() ) );
        compressed_spike_data[ syn_id ].push_back( spike_data );

        it = compressible_sources_[ tid ][ syn_id ].erase( it );
      }
      compressible_sources_[ tid ][ syn_id ].clear();
    }
  }*/
}

}
