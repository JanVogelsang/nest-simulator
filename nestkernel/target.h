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
 * This class implements a 64-bit target neuron identifier type.
 *
 * It uniquely identifies a target neuron on a (remote) machine.
 * Used in TargetTable for the presynaptic part
 * of the connection infrastructure.
 *
 * The bitwise layout of the neuron identifier for the "standard" CMAKE option:
 *
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
 * Other custom layouts can be chosen by providing a list of 5
 * numbers, representing the bits required for rank, thread, synapse
 * id, local connection id and processed flag, respectively. The number
 * of bits needs to sum to 64. The processed flag must always use one
 * bit.
 */

enum enum_status_target_id
{
  TARGET_ID_PROCESSED,
  TARGET_ID_UNPROCESSED
};

class Target
{
private:
  size_t rank_ : NUM_BITS_RANK;
  size_t tid_ : NUM_BITS_TID;
  size_t syn_id_ : NUM_BITS_SYN_ID;
  size_t lcid_ : NUM_BITS_LCID;
  size_t processed_flag_ : 1;

public:
  Target() = delete;
  Target( const Target& target );
  Target( const size_t tid, const size_t rank, const synindex syn_id, const size_t lcid );

  Target& operator=( const Target& );

  /**
   * Set local connection id.
   */
  void set_lcid( const size_t lcid );

  /**
   * Return local connection id.
   */
  size_t get_lcid() const;

  /**
   * Set rank.
   */
  void set_rank( const size_t rank );

  /**
   * Return rank.
   */
  size_t get_rank() const;

  /**
   * Set thread id.
   */
  void set_tid( const size_t tid );

  /**
   * Return thread id.
   */
  size_t get_tid() const;

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

inline Target::Target( const Target& target )
  : rank_( target.rank_ )
  , tid_( target.tid_ )
  , syn_id_( target.syn_id_ )
  , lcid_( target.lcid_ )
{
  set_status( TARGET_ID_UNPROCESSED ); // initialize
}

inline Target&
Target::operator=( const Target& other )
{
  rank_ = other.rank_;
  tid_ = other.tid_;
  syn_id_ = other.syn_id_;
  lcid_ = other.lcid_;
  set_status( TARGET_ID_UNPROCESSED );
  return *this;
}

inline Target::Target( const size_t rank, const size_t tid, const synindex syn_id, const size_t lcid )
  : rank_( rank )
  , tid_( tid )
  , syn_id_( syn_id )
  , lcid_( lcid )
{
  set_status( TARGET_ID_UNPROCESSED ); // initialize
}

inline void
Target::set_lcid( const size_t lcid )
{
  assert( lcid < MAX_LCID );
  lcid_ = lcid;
}

inline size_t
Target::get_lcid() const
{
  return lcid_;
}

inline void
Target::set_rank( const size_t rank )
{
  assert( rank <= MAX_RANK ); // MAX_RANK is allowed since it is not used as invalid value
  rank_ = rank;
}

inline size_t
Target::get_rank() const
{
  return rank_;
}

inline void
Target::set_tid( const size_t tid )
{
  assert( tid <= MAX_TID ); // MAX_TID is allowed since it is not used as invalid value
  tid_ = tid;
}

inline size_t
Target::get_tid() const
{
  return tid_;
}

inline void
Target::set_syn_id( const synindex syn_id )
{
  assert( syn_id < MAX_SYN_ID );
  syn_id_ = syn_id;
}

inline synindex
Target::get_syn_id() const
{
  return syn_id_;
}

inline void
Target::set_status( enum_status_target_id set_status_to )
{
  switch ( set_status_to )
  {
  case TARGET_ID_PROCESSED:
    processed_flag_ = TARGET_ID_PROCESSED;
    break;
  case TARGET_ID_UNPROCESSED:
    processed_flag_ = TARGET_ID_UNPROCESSED;
    break;
  default:
    throw InternalError( "Invalid remote target id status." );
  }
}

inline enum_status_target_id
Target::get_status() const
{
  return static_cast< enum_status_target_id >( processed_flag_ );
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
  OffGridTarget() = delete;
  OffGridTarget( const Target& target, const double offset );
  double get_offset() const;
};

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

#endif /* #ifndef TARGET_H */
