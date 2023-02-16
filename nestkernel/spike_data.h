/*
 *  spike_data.h
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

#ifndef SPIKE_DATA_H
#define SPIKE_DATA_H

// C++ includes
#include <cassert>

// Includes from nestkernel:
#include "nest_types.h"
#include "static_assert.h"

namespace nest
{

enum enum_status_spike_data_id
{
  SPIKE_DATA_ID_DEFAULT,
  SPIKE_DATA_ID_END,
  SPIKE_DATA_ID_COMPLETE,
  SPIKE_DATA_ID_INVALID,
};

/**
 * Used to communicate spikes. These are the elements of the MPI buffers.
 *
 * @see TargetData
 */
class SpikeData
{
  static constexpr int MAX_LAG = generate_max_value( NUM_BITS_LAG );

  unsigned int marker_ : NUM_BITS_MARKER_SPIKE_DATA;     //!< status flag
  unsigned int lag_ : NUM_BITS_LAG;                      //!< lag in this min-delay interval
  unsigned int target_tid_ : NUM_BITS_TID;                //!< target thread index
  synindex syn_id_ : NUM_BITS_SYN_ID;                    //!< synapse-type index
  unsigned int source_tid_ : NUM_BITS_TID;               //!< source thread index
  unsigned int source_lid_ : NUM_BITS_LOCAL_NODE_ID;     //!< source node local index
  // TODO JV: There still are bits to spare here, could we use them somehow?

public:
  SpikeData() = default;
  SpikeData( const thread target_tid,
    const synindex syn_id,
    const unsigned int lag,
    const thread source_tid,
    const index source_lid );
  SpikeData( const SpikeData& ) = default;
  SpikeData( SpikeData&& ) noexcept = default;
  SpikeData& operator=( const SpikeData& other ) = default;
  SpikeData& operator=( SpikeData&& other ) = default;

//  void set( const thread target_tid,
//    const synindex syn_id,
//    const unsigned int lag,
//    const thread source_tid,
//    const index source_lid);

  /**
   * Returns lag in min-delay interval.
   */
  unsigned int get_lag() const;

  /**
   * Returns target thread index.
   */
  thread get_target_tid() const;

  /**
   * Returns synapse-type index.
   */
  synindex get_syn_id() const;

  /**
   * Returns source thread index.
   */
  thread get_source_tid() const;

  /**
   * Returns source local node index.
   */
  thread get_source_lid() const;

  /**
   * Resets the status flag to default value.
   */
  void reset_marker();

  /**
   * Sets the status flag to complete marker.
   */
  void set_complete_marker();

  /**
   * Sets the status flag to end marker.
   */
  void set_end_marker();

  /**
   * Sets the status flag to invalid marker.
   */
  void set_invalid_marker();

  /**
   * Returns whether the marker is the complete marker.
   */
  bool is_complete_marker() const;

  /**
   * Returns whether the marker is the end marker.
   */
  bool is_end_marker() const;

  /**
   * Returns whether the marker is the invalid marker.
   */
  bool is_invalid_marker() const;

  /**
   * Returns offset.
   */
  double get_offset() const;
};

//! check legal size
using success_spike_data_size = StaticAssert< sizeof( SpikeData ) == 8 >::success;

//inline SpikeData::SpikeData()
//  : marker_( SPIKE_DATA_ID_DEFAULT )
//  , lag_( 0 )
//  , target_tid_( 0 )
//  , syn_id_( 0 )
//  , source_tid_( 0 )
//  , source_lid_( 0 )
//{
//}

inline SpikeData::SpikeData( const thread target_tid,
  const synindex syn_id,
  const unsigned int lag,
  const thread source_tid,
  const index source_lid )
  : marker_( SPIKE_DATA_ID_DEFAULT )
  , lag_( lag )
  , target_tid_( target_tid )
  , syn_id_( syn_id )
  , source_tid_( source_tid )
  , source_lid_( source_lid )
{
}

inline unsigned int
SpikeData::get_lag() const
{
  return lag_;
}

inline thread
SpikeData::get_target_tid() const
{
  return target_tid_;
}

inline synindex
SpikeData::get_syn_id() const
{
  return syn_id_;
}

inline thread
SpikeData::get_source_tid() const
{
  return source_tid_;
}

inline thread
SpikeData::get_source_lid() const
{
  return source_lid_;
}

//inline void
//SpikeData::set( const thread target_tid,
//  const synindex syn_id,
//  const unsigned int lag,
//  const thread source_tid,
//  const index source_lid)
//{
//  assert( 0 <= target_tid );
//  assert( target_tid <= MAX_TID );
//  assert( 0 <= source_tid );
//  assert( source_tid <= MAX_TID );
//  assert( syn_id <= MAX_SYN_ID );
//  assert( lag < MAX_LAG );
//
//  marker_ = SPIKE_DATA_ID_DEFAULT;
//  lag_ = lag;
//  target_tid_ = target_tid;
//  syn_id_ = syn_id;
//  source_tid_ = source_tid;
//  source_lid_ = source_lid;
//}

inline void
SpikeData::reset_marker()
{
  marker_ = SPIKE_DATA_ID_DEFAULT;
}

inline void
SpikeData::set_complete_marker()
{
  marker_ = SPIKE_DATA_ID_COMPLETE;
}

inline void
SpikeData::set_end_marker()
{
  marker_ = SPIKE_DATA_ID_END;
}

inline void
SpikeData::set_invalid_marker()
{
  marker_ = SPIKE_DATA_ID_INVALID;
}

inline bool
SpikeData::is_complete_marker() const
{
  return marker_ == SPIKE_DATA_ID_COMPLETE;
}

inline bool
SpikeData::is_end_marker() const
{
  return marker_ == SPIKE_DATA_ID_END;
}

inline bool
SpikeData::is_invalid_marker() const
{
  return marker_ == SPIKE_DATA_ID_INVALID;
}

inline double
SpikeData::get_offset() const
{
  return 0;
}

class OffGridSpikeData : public SpikeData
{
private:
  double offset_;

public:
  OffGridSpikeData() = default;
  OffGridSpikeData( const thread target_tid,
    const synindex syn_id,
    const unsigned int lag,
    const thread source_tid,
    const index source_lid,
    const double offset );
  OffGridSpikeData( const OffGridSpikeData& ) = default;
  OffGridSpikeData( OffGridSpikeData&& ) noexcept = default;
  OffGridSpikeData& operator=( const OffGridSpikeData& other ) = default;
  OffGridSpikeData& operator=( OffGridSpikeData&& other ) = default;

//  void set( const thread target_tid,
//    const synindex syn_id,
//    const unsigned int lag,
//    const thread source_tid,
//    const index source_lid,
//    const double offset );

  double get_offset() const;
};

//! check legal size
using success_offgrid_spike_data_size = StaticAssert< sizeof( OffGridSpikeData ) == 16 >::success;

//inline OffGridSpikeData::OffGridSpikeData()
//  : SpikeData()
//  , offset_( 0. )
//{
//}

inline OffGridSpikeData::OffGridSpikeData( const thread target_tid,
  const synindex syn_id,
  const unsigned int lag,
  const thread source_tid,
  const index source_lid,
  const double offset )
  : SpikeData( target_tid, syn_id, lag, source_tid, source_lid )
  , offset_( offset )
{
}


//inline void
//OffGridSpikeData::set( const thread _target_tid,
//  const synindex _syn_id,
//  const unsigned int _lag,
//  const thread _source_tid,
//  const index _source_lid,
//  const double _offset )
//{
//  assert( 0 <= _target_tid );
//  assert( _target_tid <= MAX_TID );
//  assert( 0 <= _source_tid );
//  assert( _source_tid <= MAX_TID );
//  assert( _syn_id <= MAX_SYN_ID );
//  assert( _lag < MAX_LAG );
//
//  marker = SPIKE_DATA_ID_DEFAULT;
//  lag = _lag;
//  target_tid = _target_tid;
//  syn_id = _syn_id;
//  offset_ = _offset;
//  source_tid = _source_tid;
//  source_lid = _source_lid;
//}


inline double
OffGridSpikeData::get_offset() const
{
  return offset_;
}

} // namespace nest

#endif /* SPIKE_DATA_H */
