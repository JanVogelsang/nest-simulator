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
 * This class implements a 64-bit target connection identifier type.
 *
 * It uniquely identifies a target neuron on a (remote) machine.
 * Used in TargetTable for the presynaptic part of the connection infrastructure.
 *
 */

enum enum_status_target_id : unsigned int
{
  TARGET_ID_PROCESSED,
  TARGET_ID_UNPROCESSED
};

class Target
{
private:
  unsigned int target_rank_ : NUM_BITS_RANK;
  enum_status_target_id status_ : 1;
  synindex syn_id_ : NUM_BITS_SYN_ID;
  unsigned int lcid_;

public:
  Target() = default;
  Target( const Target& target );
  Target( const size_t rank, const size_t lcid, const synindex syn_id );

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
   * Set synapse type id.
   */
  void set_syn_id( const synindex syn_id );

  /**
   * Return synapse type id.
   */
  size_t get_syn_id() const;

  /**
   * Set rank.
   */
  void set_rank( const size_t rank );

  /**
   * Return rank.
   */
  size_t get_rank() const;

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

inline Target::Target( const Target& other )
  : target_rank_( other.target_rank_ )
  , syn_id_( other.syn_id_ )
  , lcid_( other.lcid_ )
{
  set_status( TARGET_ID_UNPROCESSED ); // initialize
}

inline Target&
Target::operator=( const Target& other )
{
  target_rank_ = other.target_rank_;
  lcid_ = other.lcid_;
  syn_id_ = other.syn_id_;
  set_status( TARGET_ID_UNPROCESSED );
  return *this;
}

inline Target::Target( const size_t rank, const size_t lcid, const synindex syn_id )
{
  set_rank( rank );
  set_lcid( lcid );
  set_syn_id( syn_id );
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
Target::set_syn_id( const synindex syn_id )
{
  assert( syn_id < MAX_SYN_ID );
  syn_id_ = syn_id;
}

inline size_t
Target::get_syn_id() const
{
  return syn_id_;
}

inline void
Target::set_rank( const size_t rank )
{
  assert( rank <= MAX_RANK ); // MAX_RANK is allowed since it is not used as invalid value
  target_rank_ = rank;
}

inline size_t
Target::get_rank() const
{
  return target_rank_;
}

inline void
Target::set_status( enum_status_target_id set_status_to )
{
  status_ = set_status_to;
}

inline enum_status_target_id
Target::get_status() const
{
  return status_;
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

#endif /* #ifndef TARGET_H */
