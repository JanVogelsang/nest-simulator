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
#include <cstddef>
#include <vector>

// Includes from nestkernel:


#ifndef NEST_DYNAMIC_SPIKE_BUFFER_H
#define NEST_DYNAMIC_SPIKE_BUFFER_H

namespace nest
{

struct SpikeBufferEntry
{
  size_t local_connection_id;

  SpikeBufferEntry( const size_t local_connection_id ) : local_connection_id( local_connection_id ){}
};

class DynamicSpikeBuffer
{
  std::vector< std::vector< SpikeBufferEntry > > spike_buffer_;

public:
  void resize( const size_t max_delay_or_something );

  void push_back( const size_t some_index, const size_t local_connection_id );
};

inline void
DynamicSpikeBuffer::resize( const size_t max_delay_or_something )
{
  spike_buffer_.resize( max_delay_or_something );
}

inline void
DynamicSpikeBuffer::push_back( const size_t some_index, const size_t local_connection_id )
{
  spike_buffer_[ some_index ].emplace_back( local_connection_id );
}

}

#endif // NEST_DYNAMIC_SPIKE_BUFFER_H
