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
 * Used to communicate part of the connection infrastructure from
 * post- to presynaptic side. These are the elements of the MPI
 * buffers.
 *
 * SeeAlso: SpikeData
 */
class TargetData
{
  // Members must be set explicitly -- no defaults
  // Done this way to create large vector without pre-construction

private:
  static constexpr uint8_t NUM_BITS_LID = 17U;

  static constexpr int MAX_LID = generate_max_value( NUM_BITS_LID );

  size_t syn_id_ : NUM_BITS_SYN_ID;    //!< connection synapse type
  size_t source_lid_ : NUM_BITS_LID;   //!< local id of pre-synaptic neuron
  size_t source_tid_ : NUM_BITS_TID;   //!< thread index of pre-synaptic neuron
  size_t target_lcid_ : NUM_BITS_LCID; //!< first target connection local id
  size_t marker_ : 2;

public:
  void reset_marker();
  void set_complete_marker();
  void set_end_marker();
  void set_invalid_marker();
  bool is_complete_marker() const;
  bool is_end_marker() const;
  bool is_invalid_marker() const;
  void set_syn_id( const synindex syn_id );
  void set_source_lid( const size_t source_lid );
  void set_source_tid( const size_t source_tid );
  void set_target_lcid( const size_t lcid );
  synindex get_syn_id() const;
  size_t get_source_lid() const;
  size_t get_source_tid() const;
  size_t get_target_lcid() const;
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
TargetData::set_syn_id( const synindex syn_id )
{
  assert( syn_id < MAX_SYN_ID );
  syn_id_ = syn_id;
}

inline void
TargetData::set_source_lid( const size_t source_lid )
{
  assert( source_lid < MAX_LID );
  source_lid_ = source_lid;
}

inline void
TargetData::set_source_tid( const size_t source_tid )
{
  assert( source_tid < MAX_TID );
  source_tid_ = source_tid;
}

inline void
TargetData::set_target_lcid( const size_t lcid )
{
  assert( lcid < MAX_LCID );
  target_lcid_ = lcid;
}

inline synindex
TargetData::get_syn_id() const
{
  return syn_id_;
}

inline size_t
TargetData::get_source_lid() const
{
  return source_lid_;
}

inline size_t
TargetData::get_source_tid() const
{
  return source_tid_;
}

inline size_t
TargetData::get_target_lcid() const
{
  return target_lcid_;
}

} // namespace nest

#endif /* #ifndef TARGET_DATA_H */
