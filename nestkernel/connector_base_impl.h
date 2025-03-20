/*
 *  connector_base_impl.h
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

#include "connector_base.h"

// Includes from nestkernel:
#include "kernel_manager.h"

// Includes from models:
#include "weight_recorder.h"

#ifdef HAVE_SIONLIB
#include <sion.h>
#endif

#ifndef CONNECTOR_BASE_IMPL_H
#define CONNECTOR_BASE_IMPL_H

namespace nest
{

template < typename ConnectionT >
void
Connector< ConnectionT >::send_weight_event( const size_t tid,
  const unsigned int lcid,
  Event& e,
  const CommonSynapseProperties& cp )
{
  // If the pointer to the receiver node in the event is invalid,
  // the event was not sent, and a WeightRecorderEvent is therefore not created.
  if ( cp.get_weight_recorder() and e.receiver_is_valid() )
  {
    // Create new event to record the weight and copy relevant content.
    WeightRecorderEvent wr_e;
    wr_e.set_port( e.get_port() );
    wr_e.set_rport( e.get_rport() );
    wr_e.set_stamp( e.get_stamp() );
    wr_e.set_sender( e.get_sender() );
    wr_e.set_sender_node_id( kernel().connection_manager.get_source_node_id( tid, syn_id_, lcid ) );
    wr_e.set_weight( e.get_weight() );
    wr_e.set_delay_steps( e.get_delay_steps() );
    wr_e.set_receiver( *static_cast< Node* >( cp.get_weight_recorder() ) );
    // Set the node_id of the postsynaptic node as receiver node ID
    wr_e.set_receiver_node_id( e.get_receiver_node_id() );
    wr_e();
  }
}

#if defined( HAVE_SIONLIB ) && defined( HAVE_MPI )
template < typename ConnectionT >
void
Connector< ConnectionT >::dump_connections( const int sionlib_file_id,
  const size_t chunk_size,
  const size_t chunk_size_left,
  const size_t tid )
{
  size_t chars_left_in_chunk = chunk_size_left;
  size_t buf_start = chunk_size - chunk_size_left;
  char* buf = static_cast< char* >( malloc( chunk_size ) );
  unsigned char thread_id;
  unsigned int thread_lid;
  constexpr size_t num_chars_thread_id = ( sizeof( thread_id ) + sizeof( char ) / 2 ) / sizeof( char );
  constexpr size_t num_chars_thread_lid = ( sizeof( thread_lid ) + sizeof( char ) / 2 ) / sizeof( char );
  for ( ConnectionT& conn : C_ )
  {
    Node* target = conn.get_target( tid );
    thread_id = target->get_thread();
    thread_lid = target->get_thread_lid();

    assert( chunk_size - chars_left_in_chunk + 5 <= chunk_size );
    strncpy(
      buf + chunk_size - chars_left_in_chunk, reinterpret_cast< const char* >( &thread_id ), num_chars_thread_id );
    chars_left_in_chunk -= num_chars_thread_id;
    strncpy(
      buf + chunk_size - chars_left_in_chunk, reinterpret_cast< const char* >( &thread_lid ), num_chars_thread_lid );
    chars_left_in_chunk -= num_chars_thread_lid;
    if ( chars_left_in_chunk < num_chars_thread_id + num_chars_thread_lid )
    {
      chars_left_in_chunk -=
        ( sion_fwrite( buf + buf_start, sizeof( char ), chunk_size - chars_left_in_chunk - buf_start, sionlib_file_id )
          + sizeof( char ) / 2 )
        / sizeof( char );
      chars_left_in_chunk = chunk_size;
      buf_start = 0;
    }
  }
  if ( chars_left_in_chunk != chunk_size )
  {
    sion_fwrite( buf + buf_start, sizeof( char ), chunk_size - chars_left_in_chunk - buf_start, sionlib_file_id );
  }
}
#else
template < typename ConnectionT >
void
Connector< ConnectionT >::dump_connections( std::ostream& connections_out, const size_t tid )
{
  for ( ConnectionT& conn : C_ )
  {
    Node* target = conn.get_target( tid );
    const unsigned char thread_id = target->get_thread();
    const unsigned int thread_lid = target->get_thread_lid();
    connections_out.write( reinterpret_cast< const char* >( &thread_id ), sizeof( thread_id ) );
    connections_out.write( reinterpret_cast< const char* >( &thread_lid ), sizeof( thread_lid ) );

    double weight = conn.get_weight();
    double delay = conn.get_delay_ms();
  }
}
#endif

} // of namespace nest

#endif
