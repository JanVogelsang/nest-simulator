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
#include "target.h"

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
protected:
  static constexpr int MAX_LAG = generate_max_value( NUM_BITS_LAG );

  index local_target_node_id_ : NUM_BITS_LOCAL_NODE_ID;             //!< thread-local neuron index
  index local_target_connection_id_ : NUM_BITS_LOCAL_CONNECTION_ID; //!< node-local connection index
  unsigned int marker_ : NUM_BITS_MARKER_SPIKE_DATA;                //!< status flag
  unsigned int lag_ : NUM_BITS_LAG;                                 //!< lag in this min-delay interval
  unsigned int tid_ : NUM_BITS_TID;                                 //!< thread index
  synindex syn_id_ : NUM_BITS_SYN_ID;                               //!< synapse-type index

public:
  SpikeData();
  SpikeData( const SpikeData& rhs );
  SpikeData( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id,
    const unsigned int lag );

  SpikeData& operator=( const SpikeData& rhs );

  void set( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id,
    const unsigned int lag,
    const double offset );

  template < class TargetT >
  void set( const TargetT& target, const unsigned int lag );

  /**
   * Returns thread-local target neuron ID.
   */
  index get_local_target_node_id() const;

  /**
   * Returns node-local target connection ID.
   */
  index get_local_target_connection_id() const;

  /**
   * Returns lag in min-delay interval.
   */
  unsigned int get_lag() const;

  /**
   * Returns thread index.
   */
  thread get_tid() const;

  /**
   * Returns synapse-type index.
   */
  synindex get_syn_id() const;

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

inline SpikeData::SpikeData()
  : local_target_node_id_( 0 )
  , local_target_connection_id_( 0 )
  , marker_( SPIKE_DATA_ID_DEFAULT )
  , lag_( 0 )
  , tid_( 0 )
  , syn_id_( 0 )
{
}

inline SpikeData::SpikeData( const SpikeData& rhs )
  : local_target_node_id_( rhs.local_target_node_id_ )
  , local_target_connection_id_( rhs.local_target_connection_id_ )
  , marker_( SPIKE_DATA_ID_DEFAULT )
  , lag_( rhs.lag_ )
  , tid_( rhs.tid_ )
  , syn_id_( rhs.syn_id_ )
{
}

inline SpikeData::SpikeData( const thread tid,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_target_connection_id,
  const unsigned int lag )
  : local_target_node_id_( local_target_node_id )
  , local_target_connection_id_( local_target_connection_id )
  , marker_( SPIKE_DATA_ID_DEFAULT )
  , lag_( lag )
  , tid_( tid )
  , syn_id_( syn_id )
{
}

inline SpikeData&
SpikeData::operator=( const SpikeData& rhs )
{
  local_target_node_id_ = rhs.local_target_node_id_;
  local_target_connection_id_ = rhs.local_target_connection_id_;
  marker_ = SPIKE_DATA_ID_DEFAULT;
  lag_ = rhs.lag_;
  tid_ = rhs.tid_;
  syn_id_ = rhs.syn_id_;
  return *this;
}

inline void
SpikeData::set( const thread tid,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_target_connection_id,
  const unsigned int lag,
  const double )
{
  assert( 0 <= tid );
  assert( tid <= MAX_TID );
  assert( syn_id <= MAX_SYN_ID );
  assert( local_target_node_id <= MAX_LOCAL_NODE_ID );
  assert( local_target_connection_id <= MAX_LOCAL_CONNECTION_ID );
  assert( lag < MAX_LAG );

  local_target_node_id_ = local_target_node_id;
  local_target_connection_id_ = local_target_connection_id;
  marker_ = SPIKE_DATA_ID_DEFAULT;
  lag_ = lag;
  tid_ = tid;
  syn_id_ = syn_id;
}


template < class TargetT >
inline void
SpikeData::set( const TargetT& target, const unsigned int lag )
{
  // the assertions in the above function are granted by the TargetT object!
  assert( lag < MAX_LAG );
  local_target_node_id_ = target.get_local_target_node_id();
  local_target_connection_id_ = target.get_local_target_connection_id();
  marker_ = SPIKE_DATA_ID_DEFAULT;
  lag_ = lag;
  tid_ = target.get_tid();
  syn_id_ = target.get_syn_id();
}

inline index
SpikeData::get_local_target_node_id() const
{
  return local_target_node_id_;
}

inline index
SpikeData::get_local_target_connection_id() const
{
  return local_target_connection_id_;
}

inline unsigned int
SpikeData::get_lag() const
{
  return lag_;
}

inline thread
SpikeData::get_tid() const
{
  return tid_;
}

inline synindex
SpikeData::get_syn_id() const
{
  return syn_id_;
}

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
  OffGridSpikeData();
  OffGridSpikeData( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id,
    const unsigned int lag,
    const double offset );
  void set( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id,
    const unsigned int lag,
    const double offset );

  template < class TargetT >
  void set( const TargetT& target, const unsigned int lag );
  double get_offset() const;
};

//! check legal size
using success_offgrid_spike_data_size = StaticAssert< sizeof( OffGridSpikeData ) == 16 >::success;

inline OffGridSpikeData::OffGridSpikeData()
  : SpikeData()
  , offset_( 0. )
{
}

inline OffGridSpikeData::OffGridSpikeData( const thread tid,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_target_connection_id,
  const unsigned int lag,
  const double offset )
  : SpikeData( tid, syn_id, local_target_node_id, local_target_connection_id, lag )
  , offset_( offset )
{
}


inline void
OffGridSpikeData::set( const thread tid,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_target_connection_id,
  const unsigned int lag,
  const double offset )
{
  assert( tid <= MAX_TID );
  assert( syn_id <= MAX_SYN_ID );
  assert( local_target_node_id <= MAX_LOCAL_NODE_ID );
  assert( local_target_connection_id <= MAX_LOCAL_CONNECTION_ID );
  assert( lag < MAX_LAG );

  local_target_node_id_ = local_target_node_id;
  local_target_connection_id_ = local_target_connection_id;
  marker_ = SPIKE_DATA_ID_DEFAULT;
  lag_ = lag;
  tid_ = tid;
  syn_id_ = syn_id;
  offset_ = offset;
}


template < class TargetT >
inline void
OffGridSpikeData::set( const TargetT& target, const unsigned int lag )
{
  SpikeData::set( target, lag );
  offset_ = target.get_offset();
}

inline double
OffGridSpikeData::get_offset() const
{
  return offset_;
}

} // namespace nest

#endif /* SPIKE_DATA_H */
