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

//enum enum_status_target_id
//{
//  TARGET_ID_PROCESSED,
//  TARGET_ID_UNPROCESSED
//};

class Target {
protected:
  unsigned int rank_ : NUM_BITS_RANK;         //!< target rank
  unsigned int tid_ : NUM_BITS_TID;           //!< thread index
  unsigned int syn_id_ : NUM_BITS_SYN_ID;     //!< synapse id

public:
  Target();
  Target( const Target& rhs );
  Target( const thread rank,
    const thread tid,
    const synindex syn_id );

  Target& operator=( const Target& rhs );

//  void set( const thread rank,
//    const thread tid,
//    const synindex syn_id );

  /**
   * Returns target rank.
   */
  thread get_rank() const;

  /**
   * Returns target thread index.
   */
  thread get_tid() const;

  /**
   * Returns synapse-type index.
   */
  synindex get_syn_id() const;

//  /**
//   * Set target rank.
//   */
//  void set_rank( const thread rank);
//
//  /**
//   * Set target thread index.
//   */
//  void set_tid( const thread tid );
//
//  /**
//   * Set synapse-type index.
//   */
//  void set_syn_id( const synindex syn_id );
};

//! check legal size
using success_spike_data_size = StaticAssert< sizeof( Target ) == 4 >::success;

inline Target::Target()
  : rank_( 0 )
  , tid_( 0 )
  , syn_id_( 0 )
{
}

inline Target::Target( const Target& rhs )
  : rank_( rhs.rank_ )
  , tid_( rhs.tid_ )
  , syn_id_( rhs.syn_id_ )
{
}

inline Target::Target( const thread rank,
  const thread tid,
  const synindex syn_id )
  : rank_( rank )
  , tid_( tid )
  , syn_id_( syn_id )
{

}

inline Target& Target::operator=( const Target& rhs ){
  rank_ = rhs.rank_;
  tid_ = rhs.tid_;
  syn_id_ = rhs.syn_id_;
  return *this;
}
//
//void Target::set( const thread rank,
//  const thread tid,
//  const synindex syn_id )
//  {
//  assert( 0 <= rank );
//  assert( rank <= MAX_RANK );
//  assert( 0 <= tid );
//  assert( tid <= MAX_TID );
//  assert( syn_id <= MAX_SYN_ID );
//
//  rank_ = rank;
//  tid_ = tid;
//  syn_id_ = syn_id;
//}

inline thread
Target::get_rank() const
{
  return rank_;
}

inline thread
Target::get_tid() const
{
  return tid_;
}

inline synindex
Target::get_syn_id() const
{
  return syn_id_;
}

//void Target::set_rank( const thread rank)
//{
//  rank_ = rank;
//}
//
//void Target::set_tid( const thread tid )
//{
//  tid_ = tid;
//}
//
//void Target::set_syn_id( const synindex syn_id )
//{
//  syn_id_ = syn_id;
//}

// TODO JV (pt): This has been reused for local devices, but can be written much cleaner
class LocalTarget
{
private:
  uint64_t remote_target_id_;

  static constexpr uint8_t BITPOS_LOCAL_TARGET_NODE_ID = 0U;
  static constexpr uint8_t BITPOS_LOCAL_TARGET_CONNECTION_ID = NUM_BITS_LOCAL_NODE_ID;
  static constexpr uint8_t BITPOS_RANK = BITPOS_LOCAL_TARGET_CONNECTION_ID + NUM_BITS_LOCAL_CONNECTION_ID;
  static constexpr uint8_t BITPOS_TID = BITPOS_RANK + NUM_BITS_RANK;
  static constexpr uint8_t BITPOS_SYN_ID = BITPOS_TID + NUM_BITS_TID;

  // generate bit-masks used in bit-operations
  static constexpr uint64_t MASK_LOCAL_TARGET_NODE_ID =
    generate_bit_mask( NUM_BITS_LOCAL_NODE_ID, BITPOS_LOCAL_TARGET_NODE_ID );
  static constexpr uint64_t MASK_LOCAL_TARGET_CONNECTION_ID =
    generate_bit_mask( NUM_BITS_LOCAL_CONNECTION_ID, BITPOS_LOCAL_TARGET_CONNECTION_ID );
  static constexpr uint64_t MASK_RANK = generate_bit_mask( NUM_BITS_RANK, BITPOS_RANK );
  static constexpr uint64_t MASK_TID = generate_bit_mask( NUM_BITS_TID, BITPOS_TID );
  static constexpr uint64_t MASK_SYN_ID = generate_bit_mask( NUM_BITS_SYN_ID, BITPOS_SYN_ID );

public:
  LocalTarget();
  LocalTarget( const LocalTarget& target );
  LocalTarget( const thread tid,
    const thread rank,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id );

  LocalTarget& operator=( const LocalTarget& );

  /**
   * Set thread-local target neuron ID.
   */
  void set_local_target_node_id( const index local_target_node_id );

  /**
   * Set node-local target connection ID.
   */
  void set_local_target_connection_id( const index local_target_connection_id );

  /**
   * Returns thread-local target neuron ID.
   */
  index get_local_target_node_id() const;

  /**
   * Returns node-local target connection ID.
   */
  index get_local_target_connection_id() const;

  /**
   * Set rank.
   */
  void set_rank( const thread rank );

  /**
   * Return rank.
   */
  thread get_rank() const;

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
   * Return offset.
   */
  double get_offset() const;

};

//!< check legal size
using success_target_size = StaticAssert< sizeof( LocalTarget ) == 8 >::success;

inline LocalTarget::LocalTarget()
  : remote_target_id_( 0 )
{
}

inline LocalTarget::LocalTarget( const LocalTarget& target )
  : remote_target_id_( target.remote_target_id_ )
{
}

inline LocalTarget&
LocalTarget::operator=( const LocalTarget& other )
{
  remote_target_id_ = other.remote_target_id_;
  return *this;
}

inline LocalTarget::LocalTarget( const thread tid,
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
}

inline void
LocalTarget::set_local_target_node_id( const index local_target_node_id )
{
  assert( local_target_node_id <= MAX_LOCAL_NODE_ID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_LOCAL_TARGET_NODE_ID ) )
    | ( static_cast< uint64_t >( local_target_node_id ) << BITPOS_LOCAL_TARGET_NODE_ID );
}

inline index
LocalTarget::get_local_target_node_id() const
{
  return ( ( remote_target_id_ & MASK_LOCAL_TARGET_NODE_ID ) >> BITPOS_LOCAL_TARGET_NODE_ID );
}

inline void
LocalTarget::set_local_target_connection_id( const index local_target_connection_id )
{
  assert( local_target_connection_id <= MAX_LOCAL_CONNECTION_ID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_LOCAL_TARGET_CONNECTION_ID ) )
    | ( static_cast< uint64_t >( local_target_connection_id ) << BITPOS_LOCAL_TARGET_CONNECTION_ID );
}

inline index
LocalTarget::get_local_target_connection_id() const
{
  return ( ( remote_target_id_ & MASK_LOCAL_TARGET_CONNECTION_ID ) >> BITPOS_LOCAL_TARGET_CONNECTION_ID );
}

inline void
LocalTarget::set_rank( const thread rank )
{
  assert( rank <= MAX_RANK );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_RANK ) ) | ( static_cast< uint64_t >( rank ) << BITPOS_RANK );
}

inline thread
LocalTarget::get_rank() const
{
  return ( ( remote_target_id_ & MASK_RANK ) >> BITPOS_RANK );
}

inline void
LocalTarget::set_tid( const thread tid )
{
  assert( tid <= MAX_TID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_TID ) ) | ( static_cast< uint64_t >( tid ) << BITPOS_TID );
}

inline thread
LocalTarget::get_tid() const
{
  return ( ( remote_target_id_ & MASK_TID ) >> BITPOS_TID );
}

inline void
LocalTarget::set_syn_id( const synindex syn_id )
{
  assert( syn_id <= MAX_SYN_ID );
  remote_target_id_ = ( remote_target_id_ & ( ~MASK_SYN_ID ) ) | ( static_cast< uint64_t >( syn_id ) << BITPOS_SYN_ID );
}

inline synindex
LocalTarget::get_syn_id() const
{
  return ( ( remote_target_id_ & MASK_SYN_ID ) >> BITPOS_SYN_ID );
}

inline double
LocalTarget::get_offset() const
{
  return 0;
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
