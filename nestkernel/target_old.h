/*
 *  target.h
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

#ifndef TARGET_H
#define TARGET_H

// C++ includes:
#include <cassert>

// Includes from nestkernel:
#include "exceptions.h"
#include "nest_types.h"
#include "static_assert.h"

namespace nest
{

/**
 * This class implements a 64-bit target neuron identifier type. It uniquely identifies
 * a target neuron on a (remote) machine. Used in TargetTable for the presynaptic part
 * of the connection infrastructure.
 *
 * The bitwise layout of the neuron identifier for the "standard" CMAKE option:
 *  TODO JV
 *  +-------- processed flag
 *  |   +---- synapse-type id (syn_id)
 *  |   |
 *  ||----------||--thread--||---------rank----------||----local connection id (lcid)----|
 *  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000
 *  |       |  |       |  |       |  |       |  |       |  |       |  |       |  |       |
 *  63      56 55      48 47      40 39      32 31      24 23      16 15      8  7       0
 *
 * The bitwise layout of the neuron identifier for the "hpc" CMAKE option:
 *
 *  +-------- processed flag
 *  |   +---- synapse-type id (syn_id)
 *  |   |
 *  ||-----||---thread----||---------rank------------||----local connection id (lcid)----|
 *  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000  0000 0000
 *  |       |  |       |  |       |  |       |  |       |  |       |  |       |  |       |
 *  63      56 55      48 47      40 39      32 31      24 23      16 15      8  7       0
 *
 * Other custom layouts can be chosen by providing a list of 5 numbers, representing the bits required for rank, thread,
 * synapse id, local connection id and processed flag, respectively. The number of bits needs to sum to 64.
 * The processed flag must always use one bit.
 */

enum enum_status_target_id
{
  TARGET_ID_PROCESSED,
  TARGET_ID_UNPROCESSED
};

class Target
{
private:
  uint64_t remote_target_id_;

  // compile for both regular and adjacency list targets in any case, as devices always need the regular, non-adjacency-
  // list-based target version
  static constexpr uint8_t BITPOS_LOCAL_TARGET_NODE_ID = 0U;
  static constexpr uint8_t BITPOS_LOCAL_TARGET_CONNECTION_ID = NUM_BITS_LOCAL_NODE_ID;
  static constexpr uint8_t BITPOS_SYN_ID = BITPOS_LOCAL_TARGET_CONNECTION_ID + NUM_BITS_LOCAL_CONNECTION_ID;
  static constexpr uint8_t BITPOS_ADJACENCY_LIST_INDEX = 0U;
  static constexpr uint8_t BITPOS_RANK = NUM_BITS_LOCAL_NODE_ID + NUM_BITS_LOCAL_CONNECTION_ID + NUM_BITS_SYN_ID;
  static constexpr uint8_t BITPOS_TID = BITPOS_RANK + NUM_BITS_RANK;
  static constexpr uint8_t BITPOS_PROCESSED_FLAG = BITPOS_TID + NUM_BITS_TID;

  using bits_for_processed_flag = StaticAssert< NUM_BITS_PROCESSED_FLAG == 1U >::success;
  using position_of_processed_flag = StaticAssert< BITPOS_PROCESSED_FLAG == 63U >::success;

  // generate bit-masks used in bit-operations
  static constexpr uint64_t MASK_LOCAL_TARGET_NODE_ID =
    generate_bit_mask( NUM_BITS_LOCAL_NODE_ID, BITPOS_LOCAL_TARGET_NODE_ID );
  static constexpr uint64_t MASK_LOCAL_TARGET_CONNECTION_ID =
    generate_bit_mask( NUM_BITS_LOCAL_CONNECTION_ID, BITPOS_LOCAL_TARGET_CONNECTION_ID );
  static constexpr uint64_t MASK_SYN_ID = generate_bit_mask( NUM_BITS_SYN_ID, BITPOS_SYN_ID );
  static constexpr uint64_t MASK_ADJACENCY_LIST_INDEX =
    generate_bit_mask( NUM_BITS_ADJACENCY_LIST_INDEX, BITPOS_ADJACENCY_LIST_INDEX );
  static constexpr uint64_t MASK_RANK = generate_bit_mask( NUM_BITS_RANK, BITPOS_RANK );
  static constexpr uint64_t MASK_TID = generate_bit_mask( NUM_BITS_TID, BITPOS_TID );
  static constexpr uint64_t MASK_PROCESSED_FLAG = generate_bit_mask( NUM_BITS_PROCESSED_FLAG, BITPOS_PROCESSED_FLAG );

public:
  Target();
  Target( const Target& target );
  Target( const thread tid,
    const thread rank,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id );
  Target( const thread tid,
    const delay dendritic_delay,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id ); // only required for devices
  Target( const thread tid, const thread rank, const index adjacency_list_index );

  Target& operator=( const Target& );

  /**
   * Set thread-local target neuron ID. Only valid when not using adjacency lists.
   */
  void set_local_target_node_id( const index local_target_node_id );

  /**
   * Set node-local target connection ID. Only valid when not using adjacency lists.
   */
  void set_local_target_connection_id( const index local_target_connection_id );

  /**
   * Returns thread-local target neuron ID. Only valid when not using adjacency lists.
   */
  index get_local_target_node_id() const;

  /**
   * Returns node-local target connection ID. Only valid when not using adjacency lists.
   */
  index get_local_target_connection_id() const;
  /**
   * Set target adjacency list entry index. Only valid when using adjacency lists.
   */
  void set_adjacency_list_index( const index adjacency_list_index );

  /**
   * Returns target adjacency list entry index. Only valid when using adjacency lists.
   */
  index get_adjacency_list_index() const;

  /**
   * Set rank.
   */
  void set_rank( const thread rank );

  /**
   * Return rank.
   */
  thread get_rank() const;

  /**
   * Set dendritic delay (required for device nodes instead of rank).
   */
  void set_dendritic_delay( const delay dendritic_delay );

  /**
   * Get dendritic delay (required for device nodes instead of rank).
   */
  delay get_dendritic_delay() const;

  /**
   * Set thread id.
   */
  void set_tid( const thread tid );

  /**
   * Return thread id.
   */
  thread get_tid() const;

  /**
   * Set the synapse-type id.
   */
  void set_syn_id( const synindex syn_id );

  /**
   * Return synapse-type id.
   */
  synindex get_syn_id() const;

  /**
   * Set the status of the target identifier: processed or unprocessed.
   */
  void set_status( enum_status_target_id status );

  /**
   * Get the status of the target identifier: processed or unprocessed.
   */
  enum_status_target_id get_status() const;

  /**
   * Return the status od the target identifier: processed or unprocessed.
   */
  bool is_processed() const;

  /**
   * Return offset.
   */
  double get_offset() const;

  /**
   *  Set the status of the target identifier to processed
   */
  void mark_for_removal();
};

//!< check legal size
using success_target_size = StaticAssert< sizeof( Target ) == 8 >::success;

inline Target::Target()
  : remote_target_id_( 0 )
{
}

inline Target::Target( const Target& target )
  : remote_target_id_( target.remote_target_id_ )
{
  set_status( TARGET_ID_UNPROCESSED ); // initialize
}

inline Target&
Target::operator=( const Target& other )
{
  remote_target_id_ = other.remote_target_id_;
  set_status( TARGET_ID_UNPROCESSED );
  return *this;
}

inline Target::Target( const thread tid,
  const thread rank,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_target_connection_id )
  : remote_target_id_( 0 )
{
  assert( tid <= MAX_TID );
  assert( rank <= MAX_RANK );
  assert( syn_id <= MAX_SYN_ID );
  assert( local_target_node_id <= MAX_LOCAL_NODE_ID );
  assert( local_target_connection_id <= MAX_LOCAL_CONNECTION_ID );

  set_local_target_node_id( local_target_node_id );
  set_local_target_connection_id( local_target_connection_id );
  set_rank( rank );
  set_tid( tid );
  set_syn_id( syn_id );
  set_status( TARGET_ID_UNPROCESSED ); // initialize
}

inline Target::Target( const thread tid,
  const delay dendritic_delay,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_target_connection_id )
  : remote_target_id_( 0 )
{
  assert( tid <= MAX_TID );
  assert( syn_id <= MAX_SYN_ID );
  assert( local_target_node_id <= MAX_LOCAL_NODE_ID );
  assert( local_target_connection_id <= MAX_LOCAL_CONNECTION_ID );

  set_local_target_node_id( local_target_node_id );
  set_local_target_connection_id( local_target_connection_id );
  set_dendritic_delay( dendritic_delay );
  set_tid( tid );
  set_syn_id( syn_id );
  set_status( TARGET_ID_UNPROCESSED ); // initialize
}

inline Target::Target( const thread tid, const thread rank, const index adjacency_list_index )
  : remote_target_id_( 0 )
{
  assert( tid <= MAX_TID );
  assert( rank <= MAX_RANK );
  assert( adjacency_list_index <= MAX_ADJACENCY_LIST_INDEX );

  set_adjacency_list_index( adjacency_list_index );
  set_rank( rank );
  set_tid( tid );
  set_status( TARGET_ID_UNPROCESSED ); // initialize
}

inline void
Target::set_local_target_node_id( const index local_target_node_id )
{
  assert( local_target_node_id <= MAX_LOCAL_NODE_ID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_LOCAL_TARGET_NODE_ID ) )
    | ( static_cast< uint64_t >( local_target_node_id ) << BITPOS_LOCAL_TARGET_NODE_ID );
}

inline index
Target::get_local_target_node_id() const
{
  return ( ( remote_target_id_ & MASK_LOCAL_TARGET_NODE_ID ) >> BITPOS_LOCAL_TARGET_NODE_ID );
}

inline void
Target::set_local_target_connection_id( const index local_target_connection_id )
{
  assert( local_target_connection_id <= MAX_LOCAL_CONNECTION_ID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_LOCAL_TARGET_CONNECTION_ID ) )
    | ( static_cast< uint64_t >( local_target_connection_id ) << BITPOS_LOCAL_TARGET_CONNECTION_ID );
}

inline index
Target::get_local_target_connection_id() const
{
  return ( ( remote_target_id_ & MASK_LOCAL_TARGET_CONNECTION_ID ) >> BITPOS_LOCAL_TARGET_CONNECTION_ID );
}

inline void
Target::set_syn_id( const synindex syn_id )
{
  assert( syn_id <= MAX_SYN_ID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_SYN_ID ) ) | ( static_cast< uint64_t >( syn_id ) << BITPOS_SYN_ID );
}

inline synindex
Target::get_syn_id() const
{
  return ( ( remote_target_id_ & MASK_SYN_ID ) >> BITPOS_SYN_ID );
}

inline void
Target::set_adjacency_list_index( const index adjacency_list_index )
{
  assert( adjacency_list_index <= MAX_ADJACENCY_LIST_INDEX );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_ADJACENCY_LIST_INDEX ) )
    | ( static_cast< uint64_t >( adjacency_list_index ) << BITPOS_ADJACENCY_LIST_INDEX );
}

inline index
Target::get_adjacency_list_index() const
{
  return ( ( remote_target_id_ & MASK_ADJACENCY_LIST_INDEX ) >> BITPOS_ADJACENCY_LIST_INDEX );
}

inline void
Target::set_rank( const thread rank )
{
  assert( rank <= MAX_RANK );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_RANK ) ) | ( static_cast< uint64_t >( rank ) << BITPOS_RANK );
}

inline thread
Target::get_rank() const
{
  return ( ( remote_target_id_ & MASK_RANK ) >> BITPOS_RANK );
}

inline void
Target::set_dendritic_delay( const delay dendritic_delay )
{
  remote_target_id_ =
    ( remote_target_id_ & ( ~MASK_RANK ) ) | ( static_cast< uint64_t >( dendritic_delay ) << BITPOS_RANK );
}

inline delay
Target::get_dendritic_delay() const
{
  return ( ( remote_target_id_ & MASK_RANK ) >> BITPOS_RANK );
}

inline void
Target::set_tid( const thread tid )
{
  assert( tid <= MAX_TID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_TID ) ) | ( static_cast< uint64_t >( tid ) << BITPOS_TID );
}

inline thread
Target::get_tid() const
{
  return ( ( remote_target_id_ & MASK_TID ) >> BITPOS_TID );
}

inline void
Target::set_status( enum_status_target_id set_status_to )
{
  switch ( set_status_to )
  {
  case TARGET_ID_PROCESSED:
    remote_target_id_ = remote_target_id_ | MASK_PROCESSED_FLAG; // set single bit
    break;
  case TARGET_ID_UNPROCESSED:
    remote_target_id_ = remote_target_id_ & ~MASK_PROCESSED_FLAG; // clear single bit
    break;
  default:
    throw InternalError( "Invalid remote target id status." );
  }
}

inline enum_status_target_id
Target::get_status() const
{
  if ( ( remote_target_id_ & MASK_PROCESSED_FLAG ) >> BITPOS_PROCESSED_FLAG ) // test single bit
  {
    return ( TARGET_ID_PROCESSED );
  }
  return ( TARGET_ID_UNPROCESSED );
}

inline bool
Target::is_processed() const
{
  return ( get_status() == TARGET_ID_PROCESSED );
}

inline double
Target::get_offset() const
{
  return 0;
}

inline void
Target::mark_for_removal()
{
  set_status( TARGET_ID_PROCESSED );
}
class OffGridTarget : public Target
{
private:
  double offset_;

public:
  OffGridTarget();
  OffGridTarget( const Target& target, const double offset );
  double get_offset() const;
};

inline OffGridTarget::OffGridTarget()
  : Target()
  , offset_( 0 )
{
}

inline OffGridTarget::OffGridTarget( const Target& target, const double offset )
  : Target( target )
  , offset_( offset )
{
}

inline double
OffGridTarget::get_offset() const
{
  return offset_;
}

} // namespace nest

#endif // TARGET_H
