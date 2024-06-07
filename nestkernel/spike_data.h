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

//! Data for actual spikes in the communication buffers.
struct SpikeDataPayload
{
  size_t lcid : NUM_BITS_LCID;       //!< local connection index
  size_t lag : NUM_BITS_LAG;   //!< lag in this min-delay interval
  synindex syn_id : NUM_BITS_SYN_ID; //!< synapse-type index
  size_t dummy : 14;
};
//! check legal size
using success_spike_data_size = StaticAssert< sizeof( SpikeDataPayload ) == 8 >::success;

//! Meta-data sent for each source group to indicate the number of spikes sent per source group and their locations in
//! the recv-buffer.
struct SpikeDataMeta
{
  uint32_t offset;      //!< Position where source group payload starts
  uint32_t num_spikes;  //!< Number of spikes (payload) sent for this source group
};
//! check legal size
using success_spike_data_size = StaticAssert< sizeof( SpikeDataMeta ) == 8 >::success;

/**
 * Used to communicate spikes. These are the elements of the MPI buffers.
 *
 * @see TargetData
 */
class SpikeData
{
protected:
  union {
    SpikeDataPayload data;
    SpikeDataMeta meta;
  };

public:
  SpikeData();
  SpikeData( const SpikeData& rhs );
  SpikeData( const synindex syn_id, const size_t lcid, const unsigned int lag );
  SpikeData( const size_t offset, const size_t num_spikes );

  SpikeData& operator=( const SpikeData& rhs );

  // SpikeData( SpikeData&& rhs) = default;

  //! Required in connection with direct-send events.
  void set( const synindex syn_id, const size_t lcid, const unsigned int lag, const double offset );

  void set( const size_t offset, const size_t num_spikes );

  template < class TargetT >
  void set( const TargetT& target, const unsigned int lag );

  /**
   * Returns local connection ID.
   */
  size_t get_lcid() const;

  /**
   * Sets lcid value.
   *
   * @note Allows each rank to send the locally required buffer size per rank.
   */
  void set_lcid( size_t );

  /**
   * Returns lag in min-delay interval.
   */
  unsigned int get_lag() const;

  /**
   * Returns synapse-type index.
   */
  synindex get_syn_id() const;

  /**
   * Returns offset.
   */
  double get_offset() const;

    /**
   * Returns position where source group payload starts.
   */
  size_t get_source_group_offset() const;

  /**
   * Sets position where source group payload starts.
   */
  void set_source_group_offset( size_t );

    /**
   * Returns number of spikes (payload) sent for this source group.
   */
  size_t get_num_spikes() const;

  /**
   * Sets number of spikes (payload) sent for this source group.
   */
  void set_num_spikes( size_t );
};

//! check legal size
using success_spike_data_size = StaticAssert< sizeof( SpikeData ) == 8 >::success;

inline SpikeData::SpikeData()
  : data( 0, 0, 0, 0 )
{
}

inline SpikeData::SpikeData( const SpikeData& rhs )
  : data( rhs.get_lcid(), rhs.get_lag(), rhs.get_syn_id(), 0 )
{
}

inline SpikeData::SpikeData( const synindex syn_id, const size_t lcid, const unsigned int lag )
  : data( lcid, lag, syn_id, 0 )
{
}

inline SpikeData::SpikeData( const size_t offset, const size_t num_spikes )
  : meta( offset, num_spikes )
{
}

inline SpikeData&
SpikeData::operator=( const SpikeData& rhs )
{
  data.lcid = rhs.get_lcid();
  data.lag = rhs.get_lag();
  data.syn_id = rhs.get_syn_id();
  data.dummy = rhs.data.dummy;

  // meta.num_spikes = rhs.get_num_spikes();
  // meta.offset = rhs.get_source_group_offset();
  return *this;
}

inline void
SpikeData::set( const synindex syn_id, const size_t lcid, const unsigned int lag, const double )
{
  assert( syn_id < MAX_SYN_ID );
  assert( lcid < MAX_LCID );
  assert( lag < MAX_LAG );

  data.lcid = lcid;
  data.lag = lag;
  data.syn_id = syn_id;
}

inline void
SpikeData::set( const size_t offset, const size_t num_spikes )
{
  meta.offset = offset;
  meta.num_spikes = num_spikes;
}

template < class TargetT >
inline void
SpikeData::set( const TargetT& target, const unsigned int lag )
{
  // the assertions in the above function are granted by the TargetT object!
  assert( lag < MAX_LAG );
  data.lcid = target.get_lcid();
  data.lag = lag;
  data.syn_id = target.get_syn_id();
}

inline size_t
SpikeData::get_lcid() const
{
  return data.lcid;
}

inline void
SpikeData::set_lcid( size_t value )
{
  assert( value < MAX_LCID );
  data.lcid = value;
}

inline unsigned int
SpikeData::get_lag() const
{
  return data.lag;
}

inline synindex
SpikeData::get_syn_id() const
{
  return data.syn_id;
}

inline double
SpikeData::get_offset() const
{
  return 0;
}

inline size_t
SpikeData::get_source_group_offset() const
{
  return meta.offset;
}

inline void
SpikeData::set_source_group_offset( size_t offset )
{
  meta.offset = offset;
}

inline size_t
SpikeData::get_num_spikes() const
{
  return meta.num_spikes;
}

inline void
SpikeData::set_num_spikes( size_t num_spikes )
{
  meta.num_spikes = num_spikes;
}

class OffGridSpikeData : public SpikeData
{
private:
  double offset_;

public:
  OffGridSpikeData();
  OffGridSpikeData( const synindex syn_id, const size_t lcid, const unsigned int lag, const double offset );
  OffGridSpikeData( const OffGridSpikeData& rhs );
  OffGridSpikeData& operator=( const OffGridSpikeData& rhs );
  OffGridSpikeData& operator=( const SpikeData& rhs );
  void set( const synindex syn_id, const size_t lcid, const unsigned int lag, const double offset );

  template < class TargetT >
  void set( const TargetT& target, const unsigned int lag );
  double get_offset() const;
};

//! check legal size
using success_offgrid_spike_data_size = StaticAssert< sizeof( OffGridSpikeData ) == 16 >::success;

inline OffGridSpikeData::OffGridSpikeData()
  : SpikeData()
  , offset_( 0.0 )
{
}

inline OffGridSpikeData::OffGridSpikeData( const synindex syn_id,
  const size_t lcid,
  const unsigned int lag,
  const double offset )
  : SpikeData( syn_id, lcid, lag )
  , offset_( offset )
{
}

inline OffGridSpikeData::OffGridSpikeData( const OffGridSpikeData& rhs )
  : SpikeData( rhs )
  , offset_( rhs.offset_ )
{
}

inline OffGridSpikeData&
OffGridSpikeData::operator=( const OffGridSpikeData& rhs )
{
  data.lcid = rhs.get_lcid();
  data.lag = rhs.get_lag();
  data.syn_id = rhs.get_syn_id();
  offset_ = rhs.offset_;
  return *this;
}

inline OffGridSpikeData&
OffGridSpikeData::operator=( const SpikeData& rhs )
{
  // Need to use get_*() here, direct access to protected members of base-class instance is prohibited,
  // see example in https://en.cppreference.com/w/cpp/language/access.
  data.lcid = rhs.get_lcid();
  data.lag = rhs.get_lag();
  data.syn_id = rhs.get_syn_id();
  offset_ = 0;
  return *this;
}

inline void
OffGridSpikeData::set( const synindex syn_id, const size_t lcid, const unsigned int lag, const double offset )
{
  assert( syn_id < MAX_SYN_ID );
  assert( lcid < MAX_LCID );
  assert( lag < MAX_LAG );

  data.lcid = lcid;
  data.lag = lag;
  data.syn_id = syn_id;
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
