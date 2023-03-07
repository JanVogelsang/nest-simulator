/*
*  dynamic_spike_buffer.h
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

// C++ includes:
#include <algorithm>
#include <vector>

// Includes from nestkernel:
#include "nest_time.h"
#include "nest_types.h"


#ifndef NEST_DYNAMIC_SPIKE_BUFFER_H
#define NEST_DYNAMIC_SPIKE_BUFFER_H

namespace nest
{

struct SpikeBufferEntry
{
  Time t_stamp;  //!< Time when spike was emitted
  synindex syn_id : NUM_BITS_SYN_ID;  //!< Synapse type
  size_t local_connection_id : NUM_BITS_LOCAL_CONNECTION_ID;  //!< Neuron-local connection index
  delay t_syn_lag;  //!< The number of steps in the target slice until arriving at synapse
  // double offset;  //!< Precise spike offset  // TODO JV (pt): Templatize this to not always store offset

  SpikeBufferEntry( const Time t_stamp, const synindex syn_id, const size_t local_connection_id, const delay t_syn_lag )
    : t_stamp( t_stamp )
    , syn_id( syn_id )
    , local_connection_id( local_connection_id )
    , t_syn_lag( t_syn_lag )
  {}

  bool operator<( const SpikeBufferEntry& rhs ) const { return t_syn_lag < rhs.t_syn_lag; }
};

class DynamicSpikeBuffer
{
  std::vector< std::vector< SpikeBufferEntry > > spike_buffer_;

  size_t current_slice_;

  std::vector< SpikeBufferEntry >::const_iterator current_spike_;

public:
  DynamicSpikeBuffer() : spike_buffer_(), current_slice_( 0 ) {}

  void resize( const size_t max_delay_or_something );

  void push_back( const unsigned long slices_to_postpone, const Time t_stamp, const synindex syn_id, const size_t local_connection_id, const delay t_syn_lag );

  std::pair< std::vector< SpikeBufferEntry >::const_iterator &, const std::vector< SpikeBufferEntry >::const_iterator > get_next_spikes();

  void prepare_next_slice();
};

inline void
DynamicSpikeBuffer::resize( const size_t num_slices )
{
  spike_buffer_.resize( std::max( 1UL, num_slices ) );
}

inline void
DynamicSpikeBuffer::push_back( const unsigned long slices_to_postpone, const Time t_stamp, const synindex syn_id, const size_t local_connection_id, const delay t_syn_lag )
{
  spike_buffer_[ current_slice_ + slices_to_postpone ].emplace_back( t_stamp, syn_id, local_connection_id, t_syn_lag );
}

inline std::pair< std::vector< SpikeBufferEntry >::const_iterator &, const std::vector< SpikeBufferEntry >::const_iterator >
DynamicSpikeBuffer::get_next_spikes()
{
  return { current_spike_, spike_buffer_[ current_slice_ ].end() };
}

inline void
DynamicSpikeBuffer::prepare_next_slice()
{
  spike_buffer_[ current_slice_ ].clear();
  ++current_slice_;

  if ( current_slice_ == spike_buffer_.size() )  // TODO JV (pt): Profiling required
  {
    current_slice_ = 0;
  }

  std::sort(spike_buffer_[ current_slice_ ].begin(), spike_buffer_[ current_slice_ ].end() );
  current_spike_ = spike_buffer_[ current_slice_ ].begin();
}

}

#endif // NEST_DYNAMIC_SPIKE_BUFFER_H
