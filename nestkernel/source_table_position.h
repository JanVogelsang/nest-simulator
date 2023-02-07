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
// #include "kernel_manager.h"

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
   * Decreases the position.
   */
  void decrease();

  /**
   * Increases the position.
   */
  void increase();

  /**
   * Returns true if the indices point outside the SourceTable to signal that the end was reached.
   */
  bool reached_end() const;
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

inline bool
SourceTablePosition::reached_end() const
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
