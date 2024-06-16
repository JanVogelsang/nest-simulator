/*
 *  target_identifier.h
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


#ifndef TARGET_IDENTIFIER_H
#define TARGET_IDENTIFIER_H

/**
 * @file Provide classes to be used as template arguments to Connection<T>.
 */

#include "compose.hpp"
#include "kernel_manager.h"

namespace nest
{

/**
 * Class providing classic target identified information with target pointer and
 * rport.
 *
 * This class represents a connection target using a pointer to the target
 * neuron and the rport. Connection classes with this class as template argument
 * provide "full" synapses.
 *
 * See Kunkel et al, Front Neuroinform 8:78 (2014), Sec 3.3.
 */
class TargetIdentifierPtr
{

public:
  TargetIdentifierPtr( const double delay )
    : target_( nullptr )
    , rport_( 0 )
    , more_targets_( false )
    , disabled_( false )
  {
    set_delay_ms( delay );
  }

  TargetIdentifierPtr( const TargetIdentifierPtr& t ) = default;
  TargetIdentifierPtr& operator=( const TargetIdentifierPtr& t ) = default;

  void
  get_status( DictionaryDatum& d ) const
  {
    // Do nothing if called on synapse prototype
    if ( target_ )
    {
      def< long >( d, names::rport, rport_ );
      def< long >( d, names::target, target_->get_node_id() );
    }
  }

  Node*
  get_target_ptr() const
  {
    return target_;
  }

  size_t
  get_rport() const
  {
    return rport_;
  }

  void
  set_target( Node* target )
  {
    target_ = target;
  }

  void set_target( const size_t target_thread, const size_t target_lid );

  void
  set_rport( size_t rprt )
  {
    rport_ = rprt;
  }

  /**
   * Return the delay of the connection in steps
   */
  double
  get_delay() const
  {
    return delay_;
  }

  /**
   * Set the delay of the connection specified in steps
   */
  void
  set_delay( const double d )
  {
    delay_ = d;
  }

  /**
   * Return the delay of the connection in ms
   */
  double
  get_delay_ms() const
  {
    return Time::delay_steps_to_ms( delay_ );
  }

  /**
   * Set the delay of the connection specified in ms
   */
  void
  set_delay_ms( const double d )
  {
    delay_ = Time::delay_ms_to_steps( d );
  }

  void
  set_source_has_more_targets( const bool more_targets )
  {
    more_targets_ = more_targets;
  }

  bool
  source_has_more_targets() const
  {
    return more_targets_;
  }

  /**
   * Disables the synapse.
   *
   * @see is_disabled
   */
  void
  disable()
  {
    disabled_ = true;
  }

  /**
   * Returns a flag denoting if the synapse is disabled.
   *
   * @see disable
   */
  bool
  is_disabled() const
  {
    return disabled_;
  }

private:
  Node* target_;       //!< Target node
  unsigned int rport_; //!< Receiver port at the target node
  unsigned int delay_ : NUM_BITS_DELAY;
  bool more_targets_ : 1;
  bool disabled_ : 1;
};

//! check legal size
using success_target_identifier_idx_size = StaticAssert< sizeof( TargetIdentifierPtr ) == 16 >::success;

/**
 * Class providing compact (hpc) target identified by index.
 *
 * This class represents a connection target using a thread-local index, instead of a pointer. Connection classes with
 * this class as template argument provide "hpc" synapses with minimal memory requirement.
 *
 * See Kunkel et al, Front Neuroinform 8:78 (2014), Sec 3.3.
 */
class TargetIdentifierIndex
{

public:
  TargetIdentifierIndex( const double delay )
    : target_lid_( invalid_targetindex )
    , target_thread_( invalid_thread )
    , more_targets_( false )
    , disabled_( false )
  {
    set_delay_ms( delay );
  }

  TargetIdentifierIndex( const TargetIdentifierIndex& t ) = default;
  TargetIdentifierIndex& operator=( const TargetIdentifierIndex& t ) = default;


  void
  get_status( DictionaryDatum& d ) const
  {
    // Do nothing if called on synapse prototype
    if ( target_lid_ != invalid_targetindex )
    {
      def< long >( d, names::rport, rport_ );
      def< long >( d, names::target, target_lid_ );
    }
  }

  Node*
  get_target_ptr() const
  {
    assert( target_lid_ != invalid_targetindex );
    return kernel().node_manager.thread_lid_to_node( target_thread_, target_lid_ );
  }

  size_t
  get_rport() const
  {
    return rport_;
  }

  void set_target( Node* target );
  void set_target( const size_t target_thread, const size_t target_lid );

  void
  set_rport( size_t rprt )
  {
    rport_ = rprt;
  }

  /**
   * Return the delay of the connection in steps
   */
  double
  get_delay() const
  {
    return delay_;
  }

  /**
   * Set the delay of the connection specified in steps
   */
  void
  set_delay( const double d )
  {
    delay_ = d;
  }

  /**
   * Return the delay of the connection in ms
   */
  double
  get_delay_ms() const
  {
    return Time::delay_steps_to_ms( delay_ );
  }

  /**
   * Set the delay of the connection specified in ms
   */
  void
  set_delay_ms( const double d )
  {
    delay_ = Time::delay_ms_to_steps( d );
  }

  void
  set_source_has_more_targets( const bool more_targets )
  {
    more_targets_ = more_targets;
  }

  bool
  source_has_more_targets() const
  {
    return more_targets_;
  }

  /**
   * Disables the synapse.
   *
   * @see is_disabled
   */
  void
  disable()
  {
    disabled_ = true;
  }

  /**
   * Returns a flag denoting if the synapse is disabled.
   *
   * @see disable
   */
  bool
  is_disabled() const
  {
    return disabled_;
  }

private:
  unsigned short rport_;
  targetindex target_lid_; //!< Target node local index
  unsigned int target_thread_ : NUM_BITS_TID;
  unsigned int delay_ : NUM_BITS_DELAY;
  bool more_targets_ : 1;
  bool disabled_ : 1;
};

//! check legal size
using success_target_identifier_idx_size = StaticAssert< sizeof( TargetIdentifierIndex ) == 8 >::success;

inline void
TargetIdentifierPtr::set_target( const size_t target_thread, const size_t target_lid )
{
  kernel().node_manager.ensure_valid_thread_local_ids();
  if ( target_lid > max_targetindex )
  {
    throw IllegalConnection(
      String::compose( "HPC synapses support at most %1 nodes per thread. "
                       "See Kunkel et al, Front Neuroinform 8:78 (2014), Sec 3.3.2.",
        max_targetindex ) );
  }
  target_ = kernel().node_manager.get_local_nodes( target_thread ).get_node_by_index( target_lid );
}

inline void
TargetIdentifierIndex::set_target( Node* target )
{
  kernel().node_manager.ensure_valid_thread_local_ids();

  size_t target_lid = target->get_thread_lid();
  if ( target_lid > max_targetindex )
  {
    throw IllegalConnection(
      String::compose( "HPC synapses support at most %1 nodes per thread. "
                       "See Kunkel et al, Front Neuroinform 8:78 (2014), Sec 3.3.2.",
        max_targetindex ) );
  }
  target_lid_ = target_lid;
  target_thread_ = target->get_thread();
}

inline void
TargetIdentifierIndex::set_target( const size_t target_thread, const size_t target_lid )
{
  if ( target_lid > max_targetindex )
  {
    throw IllegalConnection(
      String::compose( "HPC synapses support at most %1 nodes per thread. "
                       "See Kunkel et al, Front Neuroinform 8:78 (2014), Sec 3.3.2.",
        max_targetindex ) );
  }
  target_lid_ = target_lid;
  target_thread_ = target_thread;
}

} // namespace nest


#endif
