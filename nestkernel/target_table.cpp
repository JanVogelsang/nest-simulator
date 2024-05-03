/*
 *  target_table.cpp
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

// Includes from nestkernel:
#include "target_table.h"
#include "kernel_manager.h"

// Includes from libnestutil
#include "vector_util.h"

void
nest::TargetTable::initialize()
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  targets_.resize( num_threads );

#pragma omp parallel
  {
    const size_t tid = kernel().vp_manager.get_thread_id();
    // targets_[ tid ] = std::vector< std::vector< std::vector< Target > > >();
    targets_[ tid ] = std::vector< std::vector< Target > >();
  } // of omp parallel
}

void
nest::TargetTable::finalize()
{
  // std::vector< std::vector< std::vector< std::vector< Target > > > >().swap( targets_ );
  std::vector< std::vector< std::vector< Target > > >().swap( targets_ );
}

void
nest::TargetTable::prepare( const size_t tid )
{
  // add one to max_num_local_nodes to avoid possible overflow in case
  // of rounding errors
  const size_t num_local_nodes = kernel().node_manager.get_max_num_local_nodes() + 1;

  targets_[ tid ].resize( num_local_nodes );

  /*for ( size_t lid = 0; lid < num_local_nodes; ++lid )
  {
    // resize to maximal possible synapse-type index
    // TODO JV: This produces lots of empty vectors, compress this somehow
    targets_[ tid ][ lid ].resize( kernel().model_manager.get_num_connection_models() );
  }*/
}

void
nest::TargetTable::add_target( const size_t tid, const size_t target_rank, const TargetData& target_data )
{
  const size_t lid = target_data.get_source_lid();

  vector_util::grow( targets_[ tid ][ lid ] ); // TODO JV: Does this really make sense?

  // targets_[ tid ][ lid ][ target_data.get_syn_id() ].push_back( Target( target_rank, target_data.get_target_lcid() ) );
  targets_[ tid ][ lid ].push_back( Target( target_rank, target_data.get_target_lcid(), target_data.get_syn_id() ) );
}