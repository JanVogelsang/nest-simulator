/*
 *  target_data.h
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

#ifndef TARGET_DATA_H
#define TARGET_DATA_H

// C++ includes:
#include <limits>

// Includes from nestkernel:
#include "nest_types.h"
#include "static_assert.h"
#include "target.h"

namespace nest
{

enum enum_status_target_data_id
{
  TARGET_DATA_ID_DEFAULT,
  TARGET_DATA_ID_COMPLETE,
  TARGET_DATA_ID_END,
  TARGET_DATA_ID_INVALID
};

/**
 * Used to communicate part of the connection infrastructure from post- to presynaptic side. These are the elements of
 * the MPI buffers.
 * SeeAlso: SpikeData
 */
class TargetData
{
  // Members must be set explicitly -- no defaults
  // Done this way to create large vector without preconstruction and to handle variant fields

private:
  static constexpr uint8_t NUM_BITS_LID = 19U;
  static constexpr uint8_t NUM_BITS_MARKER = 2U;
  static constexpr uint8_t NUM_BITS_IS_PRIMARY = 1U;

  static constexpr int MAX_LID = generate_max_value( NUM_BITS_LID );

  unsigned int source_lid_ : NUM_BITS_LID; //!< local id of presynaptic neuron
  unsigned int source_tid_ : NUM_BITS_TID; //!< thread index of presynaptic neuron
  unsigned int marker_ : NUM_BITS_MARKER;
  //! TargetData has TargetDataFields else SecondaryTargetDataFields
  bool is_primary_ : NUM_BITS_IS_PRIMARY;
  unsigned int syn_id_ : NUM_BITS_SYN_ID; //!< synapse type
//  union
//  {
  thread target_tid_ : 25; //!< thread index of postsynaptic target neurons
//    unsigned int compressed_index_ : 25;
//    unsigned int secondary_recv_buffer_pos_ : 25;  // TODO JV (pt): Secondary events
//  };

public:

  void reset_marker();
  void set_complete_marker();
  void set_end_marker();
  void set_invalid_marker();
  bool is_complete_marker() const;
  bool is_end_marker() const;
  bool is_invalid_marker() const;
  void set_source_lid( const index source_lid );
  void set_source_tid( const thread source_tid );
  void set_target_tid( const thread target_tid );
  void set_syn_id( const synindex syn_id );
  // void set_secondary_recv_buffer_pos( const size_t secondary_recv_buffer_pos );
  index get_source_lid() const;
  thread get_source_tid() const;
  thread get_target_tid() const;
  synindex get_syn_id() const;
  // size_t get_secondary_recv_buffer_pos();
  void set_is_primary( const bool is_primary );
  bool is_primary() const;
};

//! check legal size
using success_target_data_size = StaticAssert< sizeof( TargetData ) == 8 >::success;

inline void
TargetData::reset_marker()
{
  marker_ = TARGET_DATA_ID_DEFAULT;
}

inline void
TargetData::set_complete_marker()
{
  marker_ = TARGET_DATA_ID_COMPLETE;
}

inline void
TargetData::set_end_marker()
{
  marker_ = TARGET_DATA_ID_END;
}

inline void
TargetData::set_invalid_marker()
{
  marker_ = TARGET_DATA_ID_INVALID;
}

inline bool
TargetData::is_complete_marker() const
{
  return marker_ == TARGET_DATA_ID_COMPLETE;
}

inline bool
TargetData::is_end_marker() const
{
  return marker_ == TARGET_DATA_ID_END;
}

inline bool
TargetData::is_invalid_marker() const
{
  return marker_ == TARGET_DATA_ID_INVALID;
}

inline void
TargetData::set_source_lid( const index source_lid )
{
  assert( source_lid < MAX_LID );
  source_lid_ = source_lid;
}

inline void
TargetData::set_source_tid( const thread source_tid )
{
  assert( source_tid < MAX_TID );
  source_tid_ = source_tid;
}

inline void
TargetData::set_target_tid( const thread target_tid )
{
  assert( target_tid < MAX_TID );
  target_tid_ = target_tid;
}

inline void
TargetData::set_syn_id( const synindex syn_id )
{
  assert( syn_id < MAX_SYN_ID );
  syn_id_ = syn_id;
}

inline index
TargetData::get_source_lid() const
{
  return source_lid_;
}

inline thread
TargetData::get_source_tid() const
{
  return source_tid_;
}

inline thread
TargetData::get_target_tid() const
{
  return target_tid_;
}

inline synindex
TargetData::get_syn_id() const
{
  return syn_id_;
}

inline void
TargetData::set_is_primary( const bool is_primary )
{
  is_primary_ = is_primary;
}

inline bool
TargetData::is_primary() const
{
  return is_primary_;
}
} // namespace nest

#endif // TARGET_DATA_H
