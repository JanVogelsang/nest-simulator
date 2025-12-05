/*
 *  source.h
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

#ifndef SOURCE_H
#define SOURCE_H

// Includes from nestkernel:
#include "nest_types.h"

#include <static_assert.h>

namespace nest
{

/**
 * Stores the node ID of a presynaptic neuron. Used in SourceTable.
 */
class Source
{
  static constexpr uint8_t NUM_BITS_LID = 16U;
  static constexpr uint8_t NUM_BITS_THREADS = 1;
  static constexpr uint8_t NUM_BITS_PROCESSES = 10;

private:
  static constexpr size_t PRIMARY_MASK = 1ULL;
  static constexpr size_t DISABLED_MASK = 1ULL << 1;
  static constexpr size_t PROCESSED_MASK = 1ULL << 2;
  static constexpr size_t LID_MASK = ( ( 1ULL << NUM_BITS_LID ) - 1 ) << 3;
  static constexpr size_t TID_MASK = ( ( 1ULL << NUM_BITS_THREADS ) - 1 ) << ( NUM_BITS_LID + 3 );
  static constexpr size_t RANK_MASK = ( ( 1ULL << NUM_BITS_PROCESSES ) - 1 ) << ( NUM_BITS_LID + NUM_BITS_THREADS + 3 );

  size_t node_id_ = 0;

public:
  Source();
  explicit Source( const size_t snode_id, const bool primary = true );
  Source( const Source source, const bool primary );
  Source( const size_t lid, const size_t tid, const size_t rank, const bool primary );

  void set_lid( const size_t lid );
  void set_tid( const size_t tid );
  void set_rank( const size_t rank );
  size_t get_lid() const;
  size_t get_tid() const;
  size_t get_rank() const;


  Source( const Source& ) = default;
  Source& operator=( const Source& ) = default;
  Source( Source&& ) = default;
  Source& operator=( Source&& ) = default;

  /**
   * Returns this Source's node ID.
   */
  size_t get_node_id() const;

  void set_processed( const bool processed );
  bool is_processed() const;

  /**
   * Sets whether Source is primary.
   */
  void set_primary( const bool primary );

  /**
   * Returns whether Source is primary.
   */
  bool is_primary() const;

  /**
   * Disables Source.
   */
  void disable();

  /**
   * Returns whether Source is disabled.
   */
  bool is_disabled() const;

  // size_t operator>>( const unsigned offset ) const { return node_id_ >> offset; }

  friend bool operator<( const Source& lhs, const Source& rhs );
  friend bool operator>( const Source& lhs, const Source& rhs );
  friend bool operator==( const Source& lhs, const Source& rhs );
};

//! check legal size
using success_source_size = StaticAssert< sizeof( Source ) == 8 >::success;

inline Source::Source()
  : node_id_( 0 )
{
  set_primary( true );
  set_processed( false );
}

inline Source::Source( const Source source, const bool primary )
  : node_id_( source.node_id_ )
{
  set_primary( primary );
}

inline Source::Source( const size_t lid, const size_t tid, const size_t rank, const bool primary )
{
  set_lid( lid );
  set_tid( tid );
  set_rank( rank );
  set_primary( primary );
}

inline void
Source::set_lid( const size_t lid )
{
  node_id_ = ( node_id_ & ~LID_MASK ) | ( ( lid << 3 ) & LID_MASK );
}

inline void
Source::set_tid( const size_t tid )
{
  node_id_ = ( node_id_ & ~TID_MASK ) | ( ( tid << ( NUM_BITS_LID + 3 ) ) & TID_MASK );
}

inline void
Source::set_rank( const size_t rank )
{
  node_id_ = ( node_id_ & ~RANK_MASK ) | ( ( rank << ( NUM_BITS_LID + NUM_BITS_THREADS + 3 ) ) & RANK_MASK );
}

inline size_t
Source::get_lid() const
{
  return ( node_id_ & LID_MASK ) >> 3;
}

inline size_t
Source::get_tid() const
{
  return ( node_id_ & TID_MASK ) >> ( NUM_BITS_LID + 3 );
}

inline size_t
Source::get_rank() const
{
  return ( node_id_ & RANK_MASK ) >> ( NUM_BITS_LID + NUM_BITS_THREADS + 3 );
}

inline void
Source::set_processed( const bool processed )
{
  if ( processed )
  {
    node_id_ |= PROCESSED_MASK;
  }
  else
  {
    node_id_ &= ~PROCESSED_MASK;
  }
}

inline bool
Source::is_processed() const
{
  return node_id_ & PROCESSED_MASK;
}

inline void
Source::disable()
{
  node_id_ |= PRIMARY_MASK;
}

inline bool
Source::is_disabled() const
{
  return node_id_ & DISABLED_MASK;
}

inline void
Source::set_primary( const bool primary )
{
  if ( primary )
  {
    node_id_ |= PRIMARY_MASK;
  }
  else
  {
    node_id_ &= ~PRIMARY_MASK;
  }
}

inline bool
Source::is_primary() const
{
  return node_id_ & PRIMARY_MASK;
}

inline bool
operator>( const Source& lhs, const Source& rhs )
{
  return operator<( rhs, lhs );
}

inline bool
operator==( const Source& lhs, const Source& rhs )
{
  return lhs.node_id_ == rhs.node_id_;
}

inline bool
operator<( const Source& lhs, const Source& rhs )
{
  return lhs.node_id_ < rhs.node_id_;
}

} // namespace nest

#endif /* #ifndef SOURCE_H */
