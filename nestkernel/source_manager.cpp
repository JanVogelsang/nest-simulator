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
#include "connector_model.h"
#include "kernel_manager.h"
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
  const thread num_threads = kernel().vp_manager.get_num_threads();
  is_cleared_.initialize( num_threads, false );
  current_positions_.resize( num_threads );
  has_source_.resize( num_threads );

  const size_t num_connection_models = kernel().model_manager.get_num_connection_models();
  for( thread tid = 0; tid < num_threads; ++tid )
  {
    has_source_[ tid ].resize( num_connection_models );
  }
  compressed_spike_data_map_.resize( num_connection_models );

//  if ( kernel().connection_manager.use_compressed_spikes() )
//  {
//    current_positions_.resize( num_threads, compressed_spike_data_map_[ 0 ].cbegin() );
//  }
//  else
//  {
//    for( thread tid = 0; tid < num_threads; ++tid )
//    {
//      current_positions_.emplace_back( has_source_[ tid ][ 0 ].cbegin() );
//    }
//  }
}

void
SourceManager::finalize()
{
  // only clear sources when there were nodes added to the simulation already
  if ( kernel().node_manager.size() > 0 ){
#pragma omp parallel
  {
    const thread tid = kernel().vp_manager.get_thread_id();
    if ( is_cleared_[ tid ].is_false() )
    {
      clear( tid );
    }
  }
  }

  std::vector< SourcePosition >().swap( current_positions_ );
  std::vector< std::map< index, size_t > >().swap( compressed_spike_data_map_ );
  std::vector< std::vector < std::vector < bool > > >().swap( has_source_ );
}

void SourceManager::resize_sources()
{
  const size_t num_connection_models = kernel().model_manager.get_num_connection_models();
  for ( thread tid = 0; tid < kernel().vp_manager.get_num_threads(); ++tid )
  {
    has_source_[ tid ].resize( num_connection_models );
  }
  compressed_spike_data_map_.resize( num_connection_models );
}

bool
SourceManager::is_cleared() const
{
  return is_cleared_.all_true();
}

thread
SourceManager::find_minimal_position() const
{
  const thread num_threads = kernel().vp_manager.get_num_threads();
  thread min_thread = num_threads;
  for ( thread tid = 0; tid < num_threads; ++tid )
  {
    if ( min_thread > current_positions_[ tid ].tid )
    {
      min_thread = current_positions_[ tid ].tid;
    }
  }
  return min_thread;
}

void
SourceManager::clean( const thread tid )
{
  if ( kernel().connection_manager.use_compressed_spikes() )
  {
    // Find the minimal thread number for which all sources have been communicated already. We can safely delete the
    // source information for threads smaller than the minimal thread
    const thread min_thread = find_minimal_position();

    if ( min_thread > tid )
    {
      clear_sources( tid );
    }
    else
    {
      // do nothing
      assert( min_thread <= tid );
    }
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

void
SourceManager::reset_entry_point( const thread tid )
{
  current_positions_[ tid ].tid = 0;
  current_positions_[ tid ].syn_id = 0;
  if ( kernel().connection_manager.use_compressed_spikes() )
  {
    current_positions_[ tid ].c_it = compressed_spike_data_map_[ 0 ].cbegin();
  }
  else
  {
    current_positions_[ tid ].it = has_source_[ 0 ][ 0 ].cbegin();
  }
  ensure_valid_source_position( tid );
}

bool
SourceManager::ensure_valid_source_position( const thread tid )
{
  SourcePosition& source_position = current_positions_[ tid ];

  if ( kernel().connection_manager.use_compressed_spikes() )
  {
    while ( source_position.c_it == compressed_spike_data_map_[ source_position.syn_id ].cend() )
    {
      ++source_position.syn_id;
      if ( source_position.syn_id == kernel().model_manager.get_num_connection_models() )
      {
        return false;
      }
      source_position.c_it = compressed_spike_data_map_[ source_position.syn_id ].cbegin();
    }
  }
  else
  {
    // TODO JV (pt): This is not ideal yet, refactor please
    if ( not seek_next_possible_source_position( tid ) )
    {
      return false;
    }
    while ( not *source_position.it )
    {
      ++source_position.it;
      // we just make sure the current position is potentially valid, but might still not result in a value of "true"
      // for this source node id
      if ( not seek_next_possible_source_position( tid ) )
      {
        return false;
      }
    }
  }
  return true;
}

bool
SourceManager::seek_next_possible_source_position( const thread tid )
{
  SourcePosition& source_position = current_positions_[ tid ];
  while ( source_position.it == has_source_[ source_position.tid ][ source_position.syn_id ].cend() )
  {
    ++source_position.syn_id;
    if ( source_position.syn_id == kernel().model_manager.get_num_connection_models() )
    {
      ++source_position.tid;
      if ( source_position.tid == kernel().vp_manager.get_num_threads() )
      {
        return false;
      }
      source_position.syn_id = 0;
    }
    source_position.it = has_source_[ source_position.tid ][ source_position.syn_id ].cbegin();
  }
  return true;
}

void
SourceManager::return_to_previous_valid_source_position( const thread tid )
{
  SourcePosition& source_position = current_positions_[ tid ];

  if ( kernel().connection_manager.use_compressed_spikes() )
  {
    while ( ++source_position.c_it == compressed_spike_data_map_[ source_position.syn_id ].cbegin() )
    {
      --source_position.syn_id;
      source_position.c_it = compressed_spike_data_map_[ source_position.syn_id ].cend();
    }
    --source_position.c_it;
  }
  else
  {
    do
    {
      while ( source_position.it == has_source_[ source_position.tid ][ source_position.syn_id ].cbegin() )
      {
        if ( source_position.syn_id == 0 )
        {
          --source_position.tid;
          source_position.syn_id = kernel().model_manager.get_num_connection_models();
        }
        --source_position.syn_id;

        source_position.it = has_source_[ source_position.tid ][ source_position.syn_id ].cend();
      }
      --source_position.it;
    } while ( not *source_position.it );
  }
}

void
SourceManager::clear( const thread tid )
{
  // TODO JV: Make sure iteration over all nodes is efficient
  const SparseNodeArray& thread_local_nodes = kernel().node_manager.get_local_nodes( tid );

  for ( SparseNodeArray::const_iterator n = thread_local_nodes.begin(); n != thread_local_nodes.end(); ++n )
  {
    n->get_node()->remove_disabled_connections();
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

bool
SourceManager::get_next_target_data( const thread tid,
  const thread rank_start,
  const thread rank_end,
  thread& source_rank,
  const std::vector< ConnectorModel* >& cm,
  TargetData& next_target_data )
{
  ensure_valid_source_position( tid );

  SourcePosition& current_position = current_positions_[ tid ];

  const bool is_primary = cm[ current_position.syn_id ]->is_primary();
  next_target_data.set_is_primary( is_primary );
  next_target_data.set_syn_id( current_position.syn_id );
  next_target_data.reset_marker();

  index next_source;
  if ( kernel().connection_manager.use_compressed_spikes() )
  {
    next_source = current_position.c_it->first;
    source_rank = kernel().mpi_manager.get_process_id_of_node_id( next_source );
    if ( rank_start <= source_rank and source_rank < rank_end )
    {
      if ( is_primary ) // primary connection, i.e., chemical synapses
      {
        // target thread id us used to communicate index in compressed_spike_data
        next_target_data.set_compressed_index( current_position.c_it->second );
        ++current_position.c_it;
      }
      else // secondary connection, e.g., gap junctions
      {
        assert( false ); // TODO JV (pt): Secondary events
        // the source rank will write to the buffer position relative to the first position from the absolute position in
        // the receive buffer
        /*const size_t relative_recv_buffer_pos = kernel().connection_manager.get_secondary_recv_buffer_position(
                                                  current_position.tid, current_position.syn_id, current_position.lcid )
          - kernel().mpi_manager.get_recv_displacement_secondary_events_in_int( source_rank );

        SecondaryTargetDataFields& secondary_fields = next_target_data.secondary_data;
        secondary_fields.set_recv_buffer_pos( relative_recv_buffer_pos );
        secondary_fields.set_syn_id( current_position.syn_id );*/
      }
    }
  }
  else
  {
    next_source = current_position.it - has_source_[ current_position.tid ][ current_position.syn_id ].cbegin();
    source_rank = kernel().mpi_manager.get_process_id_of_node_id( next_source );
    next_target_data.set_target_tid( current_position.tid );
    ++current_position.it;
  }

  // TODO JV (pt): Structural plasticity: Handle deactivation of sources
  next_target_data.set_source_lid( kernel().vp_manager.node_id_to_lid( next_source ) );
  next_target_data.set_source_tid( kernel().vp_manager.vp_to_thread( kernel().vp_manager.node_id_to_vp( next_source ) ) );

  return false; // reached the end of all sources
}

void
SourceManager::fill_compressed_spike_data( std::vector < std::vector< std::vector< thread > > >& compressed_spike_data )
{
  const thread num_threads = has_source_.size();
  for ( thread tid = 0; tid < num_threads; ++tid )
  {
    for ( synindex syn_id = 0; syn_id < has_source_[ tid ].size(); ++syn_id )
    {
      std::vector<bool>& sources_mask = has_source_[ tid ][ syn_id ];
      for ( index source_node_id = 0; source_node_id < sources_mask.size(); ++source_node_id )
      {
        // TODO JV (pt): Somehow get rid of this branch-mania. This might not be as bad as it looks, though, as this could
        //  potentially be predictable by the compiler.
        if ( sources_mask[ source_node_id ] )
        {
          const auto compressed_spike_data_index = compressed_spike_data_map_[ syn_id ].find( source_node_id );
          if ( compressed_spike_data_index != compressed_spike_data_map_[ syn_id ].end() )
          {
            compressed_spike_data[ syn_id ][ (*compressed_spike_data_index).second ].push_back( tid );
          }
          else
          {
            compressed_spike_data_map_[ syn_id ][ source_node_id ] = compressed_spike_data.size();
            compressed_spike_data[ syn_id ].emplace_back( std::initializer_list< thread >{ tid } );
          }
        }
      }
    }
  }
}

}
