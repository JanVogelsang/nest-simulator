/*
 *  ring_buffer_impl.h
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

#ifndef RING_BUFFER_IMPL_H
#define RING_BUFFER_IMPL_H

#include "ring_buffer.h"

template < unsigned int num_channels >
nest::MultiChannelInputBuffer< num_channels >::MultiChannelInputBuffer()
  : buffer_( kernel().vp_manager.get_num_threads() )
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    buffer_[ tid ] = std::vector< std::array< double, num_channels > >(
      kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay() );
  }
}

template < unsigned int num_channels >
void
nest::MultiChannelInputBuffer< num_channels >::resize()
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  buffer_.resize( num_threads );
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    const size_t size = kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
    if ( buffer_[ tid ].size() != size )
    {
      buffer_[ tid ].resize( size );
    }
  }
}

template < unsigned int num_channels >
void
nest::MultiChannelInputBuffer< num_channels >::clear()
{
  resize(); // does nothing if size is fine
  // set all elements to 0.0
  for ( size_t slot = 0; slot < buffer_[ 0 ].size(); ++slot )
  {
    reset_values_all_channels( slot );
  }
}

#endif
