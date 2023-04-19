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

#ifndef DYNAMIC_SPIKE_BUFFER_H
#define DYNAMIC_SPIKE_BUFFER_H

// C++ includes:
#include <algorithm>
#include <vector>

// Includes from nestkernel:
#include "kernel_manager.h"
#include "nest_time.h"
#include "nest_types.h"
#include "static_assert.h"

namespace nest
{

struct SpikeBufferEntry
{
  Time t_stamp; //!< Time when spike was emitted
  delay axonal_delay : NUM_BITS_DELAY;
  synindex syn_id : NUM_BITS_SYN_ID;                         //!< Synapse type
  size_t local_connection_id : NUM_BITS_LOCAL_CONNECTION_ID; //!< Neuron-local connection index
  delay t_syn_lag : 9; //!< The number of steps in the target slice until arriving at synapse
  // double offset;  //!< Precise spike offset  // TODO JV (pt): Templatize this to not always store offset

  SpikeBufferEntry( const Time t_stamp,
    const delay axonal_delay,
    const synindex syn_id,
    const size_t local_connection_id,
    const delay t_syn_lag )
    : t_stamp( t_stamp )
    , axonal_delay( axonal_delay )
    , syn_id( syn_id )
    , local_connection_id( local_connection_id )
    , t_syn_lag( t_syn_lag )
  {
  }

  bool
  operator<( const SpikeBufferEntry& rhs ) const
  {
    return t_syn_lag < rhs.t_syn_lag;
  }
};

//! check legal size
using success_spike_buffer_entry_size = StaticAssert< sizeof( SpikeBufferEntry ) == 16 >::success;


class DynamicSpikeBuffer
{
  std::vector< std::vector< SpikeBufferEntry > > spike_buffer_;

  size_t current_slice_;

  std::vector< SpikeBufferEntry >::const_iterator current_spike_;

public:
  DynamicSpikeBuffer()
    : spike_buffer_( 1 )
    , current_slice_( 0 )
    , current_spike_( spike_buffer_[ 0 ].cbegin() )
  {
  }

  void resize( const delay max_axonal_delay );

  void push_back( const unsigned long slices_to_postpone,
    const Time t_stamp,
    const delay axonal_delay,
    const synindex syn_id,
    const size_t local_connection_id,
    const delay t_syn_lag );

  std::pair< std::vector< SpikeBufferEntry >::const_iterator&, const std::vector< SpikeBufferEntry >::const_iterator >
  get_next_spikes();

  void prepare_next_slice();
  void clean_slice();
  void increase_slice_index();
};

inline void
DynamicSpikeBuffer::resize( const delay max_axonal_delay )
{
  const delay min_delay = kernel().connection_manager.get_min_delay();
  const size_t num_slices = ( max_axonal_delay + min_delay - 1 ) / min_delay + 1;
  spike_buffer_.resize( num_slices );
  current_slice_ = 0;
  current_spike_ = spike_buffer_[ 0 ].cbegin();
}

inline void
DynamicSpikeBuffer::push_back( const unsigned long slices_to_postpone,
  const Time t_stamp,
  const delay axonal_delay,
  const synindex syn_id,
  const size_t local_connection_id,
  const delay t_syn_lag )
{
  const size_t spike_buffer_index = ( current_slice_ + slices_to_postpone ) % spike_buffer_.size();
  // insert into sorted vector at correct position to keep it sorted
  spike_buffer_[ spike_buffer_index ].emplace_back( t_stamp, axonal_delay, syn_id, local_connection_id, t_syn_lag );
}

inline std::pair< std::vector< SpikeBufferEntry >::const_iterator&,
  const std::vector< SpikeBufferEntry >::const_iterator >
DynamicSpikeBuffer::get_next_spikes()
{
  return { current_spike_, spike_buffer_[ current_slice_ ].cend() };
}

inline void
DynamicSpikeBuffer::prepare_next_slice()
{
  current_spike_ = spike_buffer_[ current_slice_ ].cbegin();
  std::sort( spike_buffer_[ current_slice_ ].begin(),
    spike_buffer_[ current_slice_ ].end(),
    []( const SpikeBufferEntry& lhs, const SpikeBufferEntry& rhs ) { return lhs.t_syn_lag < rhs.t_syn_lag; } );
}

inline void
DynamicSpikeBuffer::clean_slice()
{
  // TODO JV: It might make sense to free the memory here instead of just clearing
  // spike_buffer_[ current_slice_ ].clear();
  std::vector< SpikeBufferEntry >().swap( spike_buffer_[ current_slice_ ] );
  // TODO JV: Instead of clearing completely, it might make sense to just reduce the size by a certain factor to
  //  not allocate new space all the time again, but at the same time freeing memory after phases of high activity
  // spike_buffer_[ current_slice_ ].clear();
  // spike_buffer_[ current_slice_ ].resize( spike_buffer_[ current_slice_ ] / 2 );
}

inline void
DynamicSpikeBuffer::increase_slice_index()
{
  ++current_slice_;
  if ( current_slice_ == spike_buffer_.size() )
  {
    current_slice_ = 0;
  }
}

}

#endif // DYNAMIC_SPIKE_BUFFER_H
