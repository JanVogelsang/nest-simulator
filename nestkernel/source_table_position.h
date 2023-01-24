/*
 *  source_table_position.h
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

#ifndef SOURCE_TABLE_POSITION_H
#define SOURCE_TABLE_POSITION_H

// Includes from nestkernel:
#include "kernel_manager.h"

namespace nest
{

/**
 * Three-tuple to store position in 3d vector of sources.
 */
struct SourceTablePosition
{
  long tid;                        //!< thread index
  long syn_id;                     //!< synapse-type index
  long local_target_node_id;       //!< thread-local target node index
  long local_target_connection_id; //!< node-local target connection index

  SourceTablePosition();
  SourceTablePosition( const long tid,
    const long syn_id,
    const long local_target_node_id,
    const long local_target_connection_id );
  SourceTablePosition( const SourceTablePosition& rhs ) = default;
  SourceTablePosition& operator=( const SourceTablePosition& rhs ) = default;


  /**
   * Decreases the index.
   */
  void decrease();

  /**
   * Increases the index.
   */
  void increase();

  /**
   * Returns true if the indices point outside the SourceTable, e.g.,
   * to signal that the end was reached.
   */
  bool is_invalid() const;
};

inline SourceTablePosition::SourceTablePosition()
  : tid( -1 )
  , syn_id( -1 )
  , local_target_node_id( -1 )
  , local_target_connection_id( -1 )
{
}

inline SourceTablePosition::SourceTablePosition( const long tid,
  const long syn_id,
  const long local_target_node_id,
  const long local_target_connection_id )
  : tid( tid )
  , syn_id( syn_id )
  , local_target_node_id( local_target_node_id )
  , local_target_connection_id( local_target_connection_id )
{
}

inline void
SourceTablePosition::decrease()
{
  // first try finding a valid index by only decreasing target-local connection id
  --local_target_connection_id;
  if ( local_target_connection_id < 0 )
  {
    // then try finding a valid index by decreasing synapse index
    --syn_id;
    if ( syn_id < 0 )
    {
      // then try finding a valid index by decreasing target node id
      --local_target_node_id;
      if ( local_target_node_id < 0 )
      {
        // then try finding a valid index by decreasing thread index
        --tid;
        if ( tid < 0 )
        {
          return; // reached the end without finding a valid entry
        }

        local_target_node_id = kernel().node_manager.get_local_nodes( tid ).size() - 1;
      }

      syn_id = kernel().model_manager.get_num_connection_models() - 1;
    }

    local_target_connection_id = kernel()
                                   .node_manager.get_local_nodes( tid )
                                   .get_node_by_index( local_target_node_id )
                                   ->get_num_conn_type_sources( syn_id )
      - 1;
  }

  // if the found index is not valid, decrease further
  if ( local_target_connection_id
    == kernel()
         .node_manager.get_local_nodes( tid )
         .get_node_by_index( local_target_node_id )
         ->get_num_conn_type_sources( syn_id ) )
  {
    decrease();
  }
}

inline void
SourceTablePosition::increase()
{
  // first try finding a valid index by only increasing target-local connection id
  ++local_target_connection_id;
  if ( local_target_connection_id
    == kernel()
         .node_manager.get_local_nodes( tid )
         .get_node_by_index( local_target_node_id )
         ->get_num_conn_type_sources( syn_id ) )
  {
    // then try finding a valid index by increasing synapse index
    ++syn_id;
    if ( syn_id == kernel().model_manager.get_num_connection_models() )
    {
      // then try finding a valid index by increasing target node id
      ++local_target_node_id;
      if ( local_target_node_id == ernel().node_manager.get_local_nodes( tid ).size() )
      {
        // then try finding a valid index by increasing thread index
        ++tid;
        if ( tid == kernel().vp_manager.get_num_threads() )
        {
          return; // reached the end without finding a valid entry
        }

        local_target_node_id = 0;
      }

      syn_id = 0;
    }

    local_target_connection_id = 0;
  }

  // if the found index is still not valid, increase further
  if ( local_target_connection_id
    == kernel()
         .node_manager.get_local_nodes( tid )
         .get_node_by_index( local_target_node_id )
         ->get_num_conn_type_sources( syn_id ) )
  {
    increase();
  }
}

inline bool
SourceTablePosition::is_invalid() const
{
  return ( tid == -1 and syn_id == -1 and local_target_node_id == -1 and local_target_connection_id == -1 );
}

inline bool
operator==( const SourceTablePosition& lhs, const SourceTablePosition& rhs )
{
  return ( lhs.tid == rhs.tid and lhs.syn_id == rhs.syn_id and lhs.local_target_node_id == rhs.local_target_node_id
    and lhs.local_target_connection_id == rhs.local_target_connection_id );
}

inline bool
operator!=( const SourceTablePosition& lhs, const SourceTablePosition& rhs )
{
  return not operator==( lhs, rhs );
}

inline bool
operator<( const SourceTablePosition& lhs, const SourceTablePosition& rhs )
{
  if ( lhs.tid == rhs.tid )
  {
    if ( lhs.syn_id == rhs.syn_id )
    {
      if ( lhs.local_target_node_id == rhs.local_target_node_id )
      {
        return lhs.local_target_connection_id < rhs.local_target_connection_id;
      }
      else
      {
        return lhs.local_target_node_id < rhs.local_target_node_id;
      }
    }
    else
    {
      return lhs.syn_id < rhs.syn_id;
    }
  }
  else
  {
    return lhs.tid < rhs.tid;
  }
}

inline bool
operator>( const SourceTablePosition& lhs, const SourceTablePosition& rhs )
{
  return operator<( rhs, lhs );
}

inline bool
operator<=( const SourceTablePosition& lhs, const SourceTablePosition& rhs )
{
  return not operator>( lhs, rhs );
}

inline bool
operator>=( const SourceTablePosition& lhs, const SourceTablePosition& rhs )
{
  return not operator<( lhs, rhs );
}

} // namespace nest

#endif // SOURCE_TABLE_POSITION_H
