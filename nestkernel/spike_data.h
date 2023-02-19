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
//#include "kernel_manager.h"
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

struct SingleTargetSpikeData {
  index local_target_node_id : NUM_BITS_LOCAL_NODE_ID;               //!< thread-local neuron index
  index local_target_connection_id : NUM_BITS_LOCAL_CONNECTION_ID;   //!< node-local connection index
  unsigned int marker : NUM_BITS_MARKER_SPIKE_DATA;                  //!< status flag
  unsigned int lag : NUM_BITS_LAG;                                   //!< lag in this min-delay interval
  thread tid : NUM_BITS_TID;                                         //!< thread index
  synindex syn_id : NUM_BITS_SYN_ID;                                 //!< synapse-type index
};

struct CompressedSpikeData {
  index compressed_index : NUM_BITS_COMPRESSED_ID;                   //!< index of target adjacency list entry
  unsigned int marker : NUM_BITS_MARKER_SPIKE_DATA;                  //!< status flag
  unsigned int lag : NUM_BITS_LAG;                                   //!< lag in this min-delay interval
  thread tid : NUM_BITS_TID;                                         //!< thread index
  synindex syn_id : NUM_BITS_SYN_ID;                                 //!< synapse-type index
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

  union
  {
    SingleTargetSpikeData single_target_data_;
    CompressedSpikeData compressed_data_;
  };

public:
  SpikeData();
  SpikeData( const SpikeData& rhs );
  SpikeData( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id,
    const unsigned int lag );
  SpikeData( const thread tid,
    const synindex syn_id,
    const index compressed_id,
    const unsigned int lag );

  SpikeData& operator=( const SpikeData& rhs );

  void set( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id,
    const unsigned int lag,
    const double offset );
  void set( const thread tid,
    const synindex syn_id,
    const index compressed_id,
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
   * Returns index in compressed spike data structure.
   */
  index get_compressed_index() const;

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
  : compressed_data_{
    0,
    SPIKE_DATA_ID_DEFAULT,
    0,
    0,
    0 }
{
}

inline SpikeData::SpikeData( const SpikeData& rhs )
  : compressed_data_{ rhs.compressed_data_.compressed_index
    , SPIKE_DATA_ID_DEFAULT
    , rhs.compressed_data_.lag
    , rhs.compressed_data_.tid
    , rhs.compressed_data_.syn_id }
{
}

inline SpikeData::SpikeData( const thread tid,
  const synindex syn_id,
  const index local_target_node_id,
  const index local_target_connection_id,
  const unsigned int lag )
  : single_target_data_{
    local_target_node_id,
    local_target_connection_id,
    SPIKE_DATA_ID_DEFAULT,
    lag,
    tid,
    syn_id }
{
}

inline SpikeData::SpikeData( const thread tid,
  const synindex syn_id,
  const index compressed_index,
  const unsigned int lag )
  : compressed_data_{
    compressed_index,
    SPIKE_DATA_ID_DEFAULT,
    lag,
    tid,
    syn_id }
{
}

inline SpikeData&
SpikeData::operator=( const SpikeData& rhs )
{
  compressed_data_.compressed_index = rhs.compressed_data_.compressed_index;
  compressed_data_.marker = SPIKE_DATA_ID_DEFAULT;
  compressed_data_.lag = rhs.compressed_data_.lag;
  compressed_data_.tid = rhs.compressed_data_.tid;
  compressed_data_.syn_id = rhs.compressed_data_.syn_id;
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

  single_target_data_.local_target_node_id = local_target_node_id;
  single_target_data_.local_target_connection_id = local_target_connection_id;
  single_target_data_.marker = SPIKE_DATA_ID_DEFAULT;
  single_target_data_.lag = lag;
  single_target_data_.tid = tid;
  single_target_data_.syn_id = syn_id;
}

inline void
SpikeData::set( const thread tid,
  const synindex syn_id,
  const index compressed_index,
  const unsigned int lag,
  const double )
{
  assert( 0 <= tid );
  assert( tid <= MAX_TID );
  assert( syn_id <= MAX_SYN_ID );
  assert( compressed_index <= MAX_COMPRESSED_ID );
  assert( lag < MAX_LAG );

  compressed_data_.compressed_index = compressed_index;
  compressed_data_.marker = SPIKE_DATA_ID_DEFAULT;
  compressed_data_.lag = lag;
  compressed_data_.tid = tid;
  compressed_data_.syn_id = syn_id;
}

template < class TargetT >
inline void
SpikeData::set( const TargetT& target, const unsigned int lag )
{
  // the assertions in the above function are granted by the TargetT object!
  assert( lag < MAX_LAG );
  compressed_data_.compressed_index = target.get_compressed_index();
  compressed_data_.marker = SPIKE_DATA_ID_DEFAULT;
  compressed_data_.lag = lag;
  compressed_data_.tid = target.get_tid();
  compressed_data_.syn_id = target.get_syn_id();
}

inline index
SpikeData::get_local_target_node_id() const
{
//  assert( not kernel().connection_manager.use_compressed_spikes() );

  return single_target_data_.local_target_node_id;
}

inline index
SpikeData::get_local_target_connection_id() const
{
//  assert( not kernel().connection_manager.use_compressed_spikes() );

  return single_target_data_.local_target_connection_id;
}

inline index
SpikeData::get_compressed_index() const
{
//  assert( kernel().connection_manager.use_compressed_spikes() );

  return compressed_data_.compressed_index;
}

inline unsigned int
SpikeData::get_lag() const
{
  return compressed_data_.lag;
}

inline thread
SpikeData::get_tid() const
{
  return compressed_data_.tid;
}

inline synindex
SpikeData::get_syn_id() const
{
  return compressed_data_.syn_id;
}

inline void
SpikeData::reset_marker()
{
  compressed_data_.marker = SPIKE_DATA_ID_DEFAULT;
}

inline void
SpikeData::set_complete_marker()
{
  compressed_data_.marker = SPIKE_DATA_ID_COMPLETE;
}

inline void
SpikeData::set_end_marker()
{
  compressed_data_.marker = SPIKE_DATA_ID_END;
}

inline void
SpikeData::set_invalid_marker()
{
  compressed_data_.marker = SPIKE_DATA_ID_INVALID;
}

inline bool
SpikeData::is_complete_marker() const
{
  return compressed_data_.marker == SPIKE_DATA_ID_COMPLETE;
}

inline bool
SpikeData::is_end_marker() const
{
  return compressed_data_.marker == SPIKE_DATA_ID_END;
}

inline bool
SpikeData::is_invalid_marker() const
{
  return compressed_data_.marker == SPIKE_DATA_ID_INVALID;
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
  OffGridSpikeData( const thread tid,
    const synindex syn_id,
    const index compressed_index,
    const unsigned int lag,
    const double offset );
  void set( const thread tid,
    const synindex syn_id,
    const index local_target_node_id,
    const index local_target_connection_id,
    const unsigned int lag,
    const double offset );
  void set( const thread tid,
    const synindex syn_id,
    const index compressed_index,
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

inline OffGridSpikeData::OffGridSpikeData( const thread tid,
  const synindex syn_id,
  const index compressed_index,
  const unsigned int lag,
  const double offset )
  : SpikeData( tid, syn_id, compressed_index, lag )
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

  single_target_data_.local_target_node_id = local_target_node_id;
  single_target_data_.local_target_connection_id = local_target_connection_id;
  single_target_data_.marker = SPIKE_DATA_ID_DEFAULT;
  single_target_data_.lag = lag;
  single_target_data_.tid = tid;
  single_target_data_.syn_id = syn_id;
  offset_ = offset;
}

inline void
OffGridSpikeData::set( const thread tid,
  const synindex syn_id,
  const index compressed_index,
  const unsigned int lag,
  const double offset )
{
  assert( tid <= MAX_TID );
  assert( syn_id <= MAX_SYN_ID );
  assert( compressed_index <= MAX_COMPRESSED_ID );
  assert( lag < MAX_LAG );

  compressed_data_.compressed_index = compressed_index;
  compressed_data_.marker = SPIKE_DATA_ID_DEFAULT;
  compressed_data_.lag = lag;
  compressed_data_.tid = tid;
  compressed_data_.syn_id = syn_id;
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
