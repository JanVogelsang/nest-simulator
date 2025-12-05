/*
 *  source.cpp
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

#include "source.h"
#include "kernel_manager.h"
#include "vp_manager_impl.h"

namespace nest
{

Source::Source( const size_t snode_id, const bool primary )
{
  const size_t source_vp = kernel().vp_manager.node_id_to_vp( snode_id );
  set_rank( kernel().mpi_manager.get_process_id_of_vp( source_vp ) );
  set_tid( kernel().vp_manager.vp_to_thread( source_vp ) );
  set_lid( kernel().vp_manager.node_id_to_lid( snode_id ) );
  set_primary( primary );
  set_processed( false );
}


size_t
Source::get_node_id() const
{
  const VPManager& vp_manager = kernel().vp_manager;
  return vp_manager.lid_to_node_id( get_lid(), vp_manager.thread_to_vp( get_tid(), get_rank() ) );
}

} // namespace nest
