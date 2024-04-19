/*
 *  ring_buffer.cpp
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

#include "ring_buffer.h"

nest::RingBuffer::RingBuffer()
  : buffer_( kernel().vp_manager.get_num_threads() )
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    buffer_[ tid ] = std::vector< double >(
      kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay(), 0.0 );
  }
}

void
nest::RingBuffer::resize()
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  buffer_.resize( num_threads );
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    size_t size = kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
    if ( buffer_[ tid ].size() != size )
    {
      buffer_[ tid ].resize( size );
    }
  }
}

void
nest::RingBuffer::clear()
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  resize(); // does nothing if size is fine
            // clear all elements
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    buffer_[ tid ].assign( buffer_[ tid ].size(), 0.0 );
  }
}


nest::MultRBuffer::MultRBuffer()
  : buffer_( kernel().vp_manager.get_num_threads() )
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    buffer_[ tid ] = std::vector< double >(
      kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay(), 0.0 );
  }
}

void
nest::MultRBuffer::resize()
{

  const size_t num_threads = kernel().vp_manager.get_num_threads();
  buffer_.resize( num_threads );
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    size_t size = kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
    if ( buffer_[ tid ].size() != size )
    {
      buffer_[ tid ].resize( size );
    }
  }
}

void
nest::MultRBuffer::clear()
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  // clear all elements
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    buffer_[ tid ].assign( buffer_[ tid ].size(), 0.0 );
  }
}


nest::ListRingBuffer::ListRingBuffer()
  : buffer_( kernel().vp_manager.get_num_threads() )
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    buffer_[ tid ] = std::vector< std::list< double > >(
      kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay() );
  }
}

void
nest::ListRingBuffer::resize()
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  buffer_.resize( num_threads );
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    size_t size = kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
    if ( buffer_[ tid ].size() != size )
    {
      buffer_[ tid ].resize( size );
    }
  }
}

void
nest::ListRingBuffer::clear()
{
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  resize(); // does nothing if size is fine
            // clear all elements
  for ( size_t tid = 0; tid != num_threads; ++tid )
  {
    for ( unsigned int i = 0; i < buffer_[ tid ].size(); i++ )
    {
      buffer_[ tid ][ i ].clear();
    }
  }
}
