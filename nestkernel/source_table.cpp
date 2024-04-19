/*
 *  source_table.cpp
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

// C++ includes:
#include <iostream>

// Includes from nestkernel:
#include "connection_manager.h"
#include "connection_manager_impl.h"
#include "kernel_manager.h"
#include "mpi_manager_impl.h"
#include "source_table.h"
#include "vp_manager_impl.h"

nest::SourceTable::SourceTable()
{
}

nest::SourceTable::~SourceTable()
{
}

void
nest::SourceTable::initialize()
{
  assert( sizeof( Source ) == 8 );
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  sources_.resize( num_threads );
  is_cleared_.initialize( num_threads, false );
  saved_entry_point_.initialize( num_threads, false );
  current_positions_.resize( num_threads );
  saved_positions_.resize( num_threads );

#pragma omp parallel
  {
    const size_t tid = kernel().vp_manager.get_thread_id();
    sources_.at( tid ).resize( 0 );
    resize_sources();
  } // of omp parallel
}

void
nest::SourceTable::finalize()
{
  for ( size_t tid = 0; tid < static_cast< size_t >( sources_.size() ); ++tid )
  {
    if ( is_cleared_[ tid ].is_false() )
    {
      clear( tid );
    }
  }

  sources_.clear();
  current_positions_.clear();
  saved_positions_.clear();
}

bool
nest::SourceTable::is_cleared() const
{
  return is_cleared_.all_true();
}

std::vector< BlockVector< nest::Source > >&
nest::SourceTable::get_thread_local_sources( const size_t tid )
{
  return sources_[ tid ];
}

nest::SourceTablePosition
nest::SourceTable::find_maximal_position() const
{
  SourceTablePosition max_position( -1, -1, -1 );
  for ( size_t tid = 0; tid < kernel().vp_manager.get_num_threads(); ++tid )
  {
    if ( max_position < saved_positions_[ tid ] )
    {
      max_position = saved_positions_[ tid ];
    }
  }
  return max_position;
}

void
nest::SourceTable::clean( const size_t tid )
{
  // Find maximal position in source table among threads to make sure
  // unprocessed entries are not removed. Given this maximal position,
  // we can safely delete all larger entries since they will not be
  // touched any more.
  const SourceTablePosition max_position = find_maximal_position();

  // If this thread corresponds to max_position's thread, we can only
  // delete part of the sources table, with indices larger than those
  // in max_position; if this thread is larger than max_positions's
  // thread, we can delete all sources; otherwise we do nothing.
  if ( max_position.tid == static_cast< long >( tid ) )
  {
    for ( synindex syn_id = max_position.syn_id; syn_id < sources_[ tid ].size(); ++syn_id )
    {
      BlockVector< Source >& sources = sources_[ tid ][ syn_id ];
      if ( max_position.syn_id == syn_id )
      {
        // we need to add 2 to max_position.lcid since
        // max_position.lcid + 1 can contain a valid entry which we
        // do not want to delete.
        if ( max_position.lcid + 2 < static_cast< long >( sources.size() ) )
        {
          sources.erase( sources.begin() + max_position.lcid + 2, sources.end() );
        }
      }
      else
      {
        assert( max_position.syn_id < syn_id );
        sources.clear();
      }
    }
  }
  else if ( max_position.tid < static_cast< long >( tid ) )
  {
    sources_[ tid ].clear();
  }
  else
  {
    // do nothing
    assert( static_cast< long >( tid ) < max_position.tid );
  }
}

size_t
nest::SourceTable::get_node_id( const size_t tid, const synindex syn_id, const size_t lcid ) const
{
  if ( not kernel().connection_manager.get_keep_source_table() )
  {
    throw KernelException( "Cannot use SourceTable::get_node_id when get_keep_source_table is false" );
  }
  return sources_[ tid ][ syn_id ][ lcid ].get_node_id();
}

size_t
nest::SourceTable::remove_disabled_sources( const size_t tid, const synindex syn_id )
{
  if ( sources_[ tid ].size() <= syn_id )
  {
    return invalid_index; // no source table entry for this synapse model
  }

  BlockVector< Source >& mysources = sources_[ tid ][ syn_id ];
  const size_t max_size = mysources.size();
  if ( max_size == 0 )
  {
    return invalid_index; // no connections for this synapse model
  }

  // lcid needs to be signed, to allow lcid >= 0 check in while loop
  // to fail; afterwards we can be certain that it is non-negative and
  // we can static_cast it to index
  long lcid = max_size - 1;
  while ( lcid >= 0 and mysources[ lcid ].is_disabled() )
  {
    --lcid;
  }
  const size_t first_invalid_lcid = static_cast< size_t >( lcid + 1 ); // loop stopped on first valid entry or -1
  if ( first_invalid_lcid == max_size )
  {
    return invalid_index; // all lcids are valid, nothing to remove
  }

  mysources.erase( mysources.begin() + first_invalid_lcid, mysources.end() );
  return first_invalid_lcid;
}

void
nest::SourceTable::resize_sources()
{
  kernel().vp_manager.assert_thread_parallel();
  sources_.at( kernel().vp_manager.get_thread_id() ).resize( kernel().model_manager.get_num_connection_models() );
}

bool
nest::SourceTable::source_should_be_processed_( const size_t rank_start,
  const size_t rank_end,
  const Source& source ) const
{
  const size_t source_rank = kernel().mpi_manager.get_process_id_of_node_id( source.get_node_id() );

  return not( source.is_processed()
    or source.is_disabled()
    // is this thread responsible for this part of the MPI buffer?
    or source_rank < rank_start or rank_end <= source_rank );
}

bool
nest::SourceTable::next_entry_has_same_source_( const SourceTablePosition& current_position,
  const Source& current_source ) const
{
  assert( not current_position.is_invalid() );

  const auto& local_sources = sources_[ current_position.tid ][ current_position.syn_id ];
  const size_t next_lcid = current_position.lcid + 1;

  return (
    next_lcid < local_sources.size() and local_sources[ next_lcid ].get_node_id() == current_source.get_node_id() );
}

bool
nest::SourceTable::previous_entry_has_same_source_( const SourceTablePosition& current_position,
  const Source& current_source ) const
{
  assert( not current_position.is_invalid() );

  const auto& local_sources = sources_[ current_position.tid ][ current_position.syn_id ];
  const long previous_lcid = current_position.lcid - 1; // needs to be a signed type such that negative
                                                        // values can signal invalid indices

  return ( previous_lcid >= 0 and not local_sources[ previous_lcid ].is_processed()
    and local_sources[ previous_lcid ].get_node_id() == current_source.get_node_id() );
}

void
nest::SourceTable::populate_target_data_fields_( const SourceTablePosition& current_position,
  const Source& current_source,
  TargetData& next_target_data ) const
{
  const auto node_id = current_source.get_node_id();

  // set values of next_target_data
  next_target_data.set_source_lid( kernel().vp_manager.node_id_to_lid( node_id ) );
  next_target_data.set_source_tid( kernel().vp_manager.vp_to_thread( kernel().vp_manager.node_id_to_vp( node_id ) ) );
  next_target_data.reset_marker();
  next_target_data.set_syn_id( current_position.syn_id );
  next_target_data.set_target_lcid( current_position.lcid );
}

bool
nest::SourceTable::get_next_target_data( const size_t tid,
  const size_t rank_start,
  const size_t rank_end,
  size_t& source_rank,
  TargetData& next_target_data )
{
  SourceTablePosition& current_position = current_positions_[ tid ];

  if ( current_position.is_invalid() )
  {
    return false; // nothing to do here
  }

  // we stay in this loop either until we can return a valid
  // TargetData object or we have reached the end of the sources table
  while ( true )
  {
    current_position.seek_to_next_valid_index( sources_ );
    if ( current_position.is_invalid() )
    {
      return false; // reached the end of the sources table
    }

    // the current position contains an entry, so we retrieve it
    Source& current_source = sources_[ current_position.tid ][ current_position.syn_id ][ current_position.lcid ];

    if ( not source_should_be_processed_( rank_start, rank_end, current_source ) )
    {
      current_position.decrease();
      continue;
    }

    // we need to set a marker stating whether the entry following this
    // entry, if existent, has the same source
    kernel().connection_manager.set_source_has_more_targets( current_position.tid,
      current_position.syn_id,
      current_position.lcid,
      next_entry_has_same_source_( current_position, current_source ) );

    // no need to communicate this entry if the previous entry has the same source
    if ( previous_entry_has_same_source_( current_position, current_source ) )
    {
      current_source.set_processed( true ); // no need to look at this entry again
      current_position.decrease();
      continue;
    }

    // reaching this means we found an entry that should be
    // communicated via MPI, so we prepare to return the relevant data

    // set the source rank
    source_rank = kernel().mpi_manager.get_process_id_of_node_id( current_source.get_node_id() );

    populate_target_data_fields_( current_position, current_source, next_target_data );

    // we are about to return a valid entry, so mark it as processed
    current_source.set_processed( true );

    current_position.decrease();
    return true; // found a valid entry
  }
}

void
nest::SourceTable::dump_sources() const
{
  FULL_LOGGING_ONLY( for ( size_t tid = 0; tid < sources_.size(); ++tid ) {
    for ( size_t syn_id = 0; syn_id < sources_[ tid ].size(); ++syn_id )
    {
      for ( size_t lcid = 0; lcid < sources_[ tid ][ syn_id ].size(); ++lcid )
      {
        kernel().write_to_dump( String::compose( "src  : r%1 t%2 s%3 tg%4 l%5 tt%6",
          kernel().mpi_manager.get_rank(),
          kernel().vp_manager.get_thread_id(),
          sources_[ tid ][ syn_id ][ lcid ].get_node_id(),
          kernel().connection_manager.get_target_node_id( tid, syn_id, lcid ),
          lcid,
          tid ) );
      }
    }
  } )
}

