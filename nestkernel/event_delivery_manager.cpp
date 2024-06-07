/*
 *  event_delivery_manager.cpp
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

#include "event_delivery_manager.h"

// C++ includes:
#include <algorithm> // rotate
#include <numeric>   // accumulate

// Includes from nestkernel:
#include "connection_manager.h"
#include "connection_manager_impl.h"
#include "event_delivery_manager_impl.h"
#include "kernel_manager.h"
#include "mpi_manager_impl.h"
#include "send_buffer_position.h"
#include "source.h"
#include "vp_manager.h"
#include "vp_manager_impl.h"

// Includes from sli:
#include "dictutils.h"

#include "compose.hpp"

namespace nest
{


EventDeliveryManager::EventDeliveryManager()
  : off_grid_spiking_( false )
  , moduli_()
  , slice_moduli_()
  , emitted_spikes_register_()
  , off_grid_emitted_spikes_register_()
  , send_buffer_secondary_events_()
  , recv_buffer_secondary_events_()
  , local_spike_counter_()
  , send_buffer_spike_data_()
  , recv_buffer_spike_data_()
  , send_buffer_off_grid_spike_data_()
  , recv_buffer_off_grid_spike_data_()
  , send_buffer_target_data_()
  , recv_buffer_target_data_()
  , buffer_size_target_data_has_changed_( false )
  , global_max_spikes_per_rank_( 0 )
  , send_recv_buffer_shrink_limit_( 0.2 )
  , send_recv_buffer_shrink_spare_( 0.1 )
  , send_recv_buffer_grow_extra_( 0.5 )
  , send_buffer_segments_()
  , send_recv_buffer_resize_log_()
  , gather_completed_checker_()
{
}

EventDeliveryManager::~EventDeliveryManager()
{
}

void
EventDeliveryManager::initialize( const bool adjust_number_of_threads_or_rng_only )
{
  if ( not adjust_number_of_threads_or_rng_only )
  {
    init_moduli();
    reset_timers_for_preparation();
    reset_timers_for_dynamics();

    // Ensures that ResetKernel resets off_grid_spiking_
    off_grid_spiking_ = false;
    buffer_size_target_data_has_changed_ = false;
    send_recv_buffer_shrink_limit_ = 0.2;
    send_recv_buffer_shrink_spare_ = 0.1;
    send_recv_buffer_grow_extra_ = 0.5;
    send_recv_buffer_resize_log_.clear();
  }

  const size_t num_threads = kernel().vp_manager.get_num_threads();
  const size_t num_ranks = kernel().mpi_manager.get_num_processes();
  local_spike_counter_.resize( num_threads, 0 );
  reset_counters();

  send_buffer_segments_.reserve(2);
  emitted_spikes_register_.resize( num_threads );
  off_grid_emitted_spikes_register_.resize( num_threads );
  gather_completed_checker_.initialize( num_threads, false );

#pragma omp parallel
  {
    const size_t tid = kernel().vp_manager.get_thread_id();

    if ( not emitted_spikes_register_[ tid ] )
    {
      emitted_spikes_register_[ tid ] = new std::vector< std::vector< SpikeData > >( num_ranks );
    }

    if ( not off_grid_emitted_spikes_register_[ tid ] )
    {
      off_grid_emitted_spikes_register_[ tid ] = new std::vector< std::vector< OffGridSpikeData > >( num_ranks );
    }
  } // of omp parallel
}

void
EventDeliveryManager::finalize( const bool )
{
  // clear the spike buffers
  for ( auto& vec_spikedata_ptr : emitted_spikes_register_ )
  {
    delete vec_spikedata_ptr;
  }
  emitted_spikes_register_.clear(); // remove stale pointers

  for ( auto& vec_spikedata_ptr : off_grid_emitted_spikes_register_ )
  {
    delete vec_spikedata_ptr;
  }
  off_grid_emitted_spikes_register_.clear();

  send_buffer_secondary_events_.clear();
  recv_buffer_secondary_events_.clear();
  send_buffer_spike_data_.clear();
  recv_buffer_spike_data_.clear();
  send_buffer_off_grid_spike_data_.clear();
  recv_buffer_off_grid_spike_data_.clear();
  responsible_source_groups_.clear();
  thread_to_source_group_.clear();
  send_buffer_segments_.clear();
}

void
EventDeliveryManager::set_status( const DictionaryDatum& dict )
{
  updateValue< bool >( dict, names::off_grid_spiking, off_grid_spiking_ );

  double bsl = send_recv_buffer_shrink_limit_;
  if ( updateValue< double >( dict, names::spike_buffer_shrink_limit, bsl ) )
  {
    if ( bsl < 0 )
    {
      throw BadProperty( "buffer_shrink_limit >= 0 required." );
    }
    send_recv_buffer_shrink_limit_ = bsl;
  }

  double bss = send_recv_buffer_shrink_spare_;
  if ( updateValue< double >( dict, names::spike_buffer_shrink_spare, bss ) )
  {
    if ( bss < 0 or bss > 1 )
    {
      throw BadProperty( "0 <= buffer_shrink_spare <= 1 required." );
    }
    send_recv_buffer_shrink_spare_ = bss;
  }

  double bge = send_recv_buffer_grow_extra_;
  if ( updateValue< double >( dict, names::spike_buffer_grow_extra, bge ) )
  {
    if ( bge < 0 )
    {
      throw BadProperty( "buffer_grow_extra >= 0 required." );
    }
    send_recv_buffer_grow_extra_ = bge;
  }
}

void
EventDeliveryManager::get_status( DictionaryDatum& dict )
{
  def< bool >( dict, names::off_grid_spiking, off_grid_spiking_ );
  def< unsigned long >(
    dict, names::local_spike_counter, std::accumulate( local_spike_counter_.begin(), local_spike_counter_.end(), 0 ) );
  def< double >( dict, names::spike_buffer_shrink_limit, send_recv_buffer_shrink_limit_ );
  def< double >( dict, names::spike_buffer_shrink_spare, send_recv_buffer_shrink_spare_ );
  def< double >( dict, names::spike_buffer_grow_extra, send_recv_buffer_grow_extra_ );

  DictionaryDatum log_events = DictionaryDatum( new Dictionary );
  ( *dict )[ names::spike_buffer_resize_log ] = log_events;
  send_recv_buffer_resize_log_.to_dict( log_events );

#ifdef TIMER_DETAILED
  def< double >( dict, names::time_collocate_spike_data, sw_collocate_spike_data_.elapsed() );
  def< double >( dict, names::time_communicate_spike_data, sw_communicate_spike_data_.elapsed() );
  def< double >( dict, names::time_communicate_target_data, sw_communicate_target_data_.elapsed() );
#endif
}

void
EventDeliveryManager::prepare()
{
  const size_t num_source_groups_per_rank = kernel().mpi_manager.get_num_source_groups_per_rank();
  const size_t num_source_groups_per_thread = kernel().mpi_manager.get_num_source_groups_per_thread();
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  const size_t num_ranks = kernel().mpi_manager.get_num_processes();

  responsible_source_groups_.resize( num_threads );
  thread_to_source_group_.resize( num_threads );

  // Iterate over all source groups rank-wise and add the corresponding source group identifiers to the responsibility
  // list for the corresponding thread.
  size_t tid = 0;
  size_t source_group_counter = num_source_groups_per_thread;
  for ( size_t rank = 0; rank != num_ranks; ++rank )
  {
    for ( size_t s = 0; s != num_source_groups_per_rank; ++s )
    {
      assert(tid < num_threads);
      responsible_source_groups_[ tid ].emplace_back( rank, s );
      --source_group_counter;
      if ( source_group_counter == 0 )
      {
        ++tid;
        source_group_counter = num_source_groups_per_thread;
      }
    }
  }

  for ( tid = 0; tid != num_threads; ++tid )
  {
    // For each thread, set the corresponding rank-local source group index
    thread_to_source_group_[ tid ] = tid / ( num_threads / num_source_groups_per_rank );
  }
}

void
EventDeliveryManager::cleanup()
{
  responsible_source_groups_.clear();
  thread_to_source_group_.clear();
}

void
EventDeliveryManager::resize_send_recv_buffers_target_data()
{
  // compute send receive counts and allocate memory for buffers
  send_buffer_target_data_.resize( kernel().mpi_manager.get_buffer_size_target_data() );
  recv_buffer_target_data_.resize( kernel().mpi_manager.get_buffer_size_target_data() );
}

void
EventDeliveryManager::resize_send_recv_buffers_spike_data_()
{
  send_buffer_spike_data_.resize( kernel().mpi_manager.get_buffer_size_spike_data() );
  recv_buffer_spike_data_.resize( kernel().mpi_manager.get_buffer_size_spike_data() );
  send_buffer_off_grid_spike_data_.resize( kernel().mpi_manager.get_buffer_size_spike_data() );
  recv_buffer_off_grid_spike_data_.resize( kernel().mpi_manager.get_buffer_size_spike_data() );
}

void
EventDeliveryManager::configure_spike_data_buffers()
{
  assert( kernel().connection_manager.get_min_delay() != 0 );

  send_buffer_spike_data_.clear();
  send_buffer_off_grid_spike_data_.clear();

  resize_send_recv_buffers_spike_data_();
}

void
EventDeliveryManager::configure_secondary_buffers()
{
  send_buffer_secondary_events_.clear();
  send_buffer_secondary_events_.resize( kernel().mpi_manager.get_send_buffer_size_secondary_events_in_int() );
  recv_buffer_secondary_events_.clear();
  recv_buffer_secondary_events_.resize( kernel().mpi_manager.get_recv_buffer_size_secondary_events_in_int() );
}

void
EventDeliveryManager::init_moduli()
{
  long min_delay = kernel().connection_manager.get_min_delay();
  long max_delay = kernel().connection_manager.get_max_delay();
  assert( min_delay != 0 );
  assert( max_delay != 0 );

  // Ring buffers use modulos to determine where to store incoming events
  // with given time stamps, relative to the beginning of the slice in which
  // the spikes are delivered from the queue, ie, the slice after the one
  // in which they were generated. The pertaining offsets are 0..max_delay-1.
  moduli_.resize( min_delay + max_delay );

  for ( long d = 0; d < min_delay + max_delay; ++d )
  {
    moduli_[ d ] = ( kernel().simulation_manager.get_clock().get_steps() + d ) % ( min_delay + max_delay );
  }

  // Slice-based ring-buffers have one bin per min_delay steps,
  // up to max_delay.  Time is counted as for normal ring buffers.
  // The slice_moduli_ table maps time steps to these bins
  const size_t nbuff = static_cast< size_t >( std::ceil( static_cast< double >( min_delay + max_delay ) / min_delay ) );
  slice_moduli_.resize( min_delay + max_delay );
  for ( long d = 0; d < min_delay + max_delay; ++d )
  {
    slice_moduli_[ d ] = ( ( kernel().simulation_manager.get_clock().get_steps() + d ) / min_delay ) % nbuff;
  }
}

void
EventDeliveryManager::update_moduli()
{
  long min_delay = kernel().connection_manager.get_min_delay();
  long max_delay = kernel().connection_manager.get_max_delay();
  assert( min_delay != 0 );
  assert( max_delay != 0 );

  /*
   * Note that for updating the modulos, it is sufficient
   * to rotate the buffer to the left.
   */
  assert( moduli_.size() == ( size_t ) ( min_delay + max_delay ) );
  std::rotate( moduli_.begin(), moduli_.begin() + min_delay, moduli_.end() );

  // For the slice-based ring buffer, we cannot rotate the table, but
  // have to re-compute it, since max_delay_ may not be a multiple of
  // min_delay_.  Reference time is the time at the beginning of the slice.
  const size_t nbuff = static_cast< size_t >( std::ceil( static_cast< double >( min_delay + max_delay ) / min_delay ) );
  for ( long d = 0; d < min_delay + max_delay; ++d )
  {
    slice_moduli_[ d ] = ( ( kernel().simulation_manager.get_clock().get_steps() + d ) / min_delay ) % nbuff;
  }
}

void
EventDeliveryManager::reset_counters()
{
  for ( auto& spike_counter : local_spike_counter_ )
  {
    spike_counter = 0;
  }
}

void
EventDeliveryManager::reset_timers_for_preparation()
{
#ifdef TIMER_DETAILED
  sw_communicate_target_data_.reset();
#endif
}

void
EventDeliveryManager::reset_timers_for_dynamics()
{
#ifdef TIMER_DETAILED
  sw_collocate_spike_data_.reset();
  sw_communicate_spike_data_.reset();
#endif
}

void
EventDeliveryManager::write_done_marker_secondary_events_( const bool done )
{
  // write done marker at last position in every chunk
  for ( size_t rank = 0; rank < kernel().mpi_manager.get_num_processes(); ++rank )
  {
    send_buffer_secondary_events_[ kernel().mpi_manager.get_done_marker_position_in_secondary_events_send_buffer(
      rank ) ] = done;
  }
}

void
EventDeliveryManager::gather_secondary_events( const bool done )
{
  write_done_marker_secondary_events_( done );
  kernel().mpi_manager.communicate_secondary_events_Alltoallv(
    std::span( send_buffer_secondary_events_ ), std::span( recv_buffer_secondary_events_ ) );
}

bool
EventDeliveryManager::deliver_secondary_events( const size_t tid, const bool called_from_wfr_update )
{
  return kernel().connection_manager.deliver_secondary_events(
    tid, called_from_wfr_update, recv_buffer_secondary_events_ );
}

void
EventDeliveryManager::gather_spike_data()
{
  if ( off_grid_spiking_ )
  {
    gather_spike_data_(
      off_grid_emitted_spikes_register_, send_buffer_off_grid_spike_data_, recv_buffer_off_grid_spike_data_ );
  }
  else
  {
    gather_spike_data_( emitted_spikes_register_, send_buffer_spike_data_, recv_buffer_spike_data_ );
  }
}

template < typename SpikeDataT >
void
EventDeliveryManager::gather_spike_data_(
  std::vector< std::vector< std::vector< SpikeDataT > >* >& emitted_spikes_register,
  std::vector< SpikeDataT >& send_buffer,
  std::vector< SpikeDataT >& recv_buffer )
{
  const size_t num_ranks = kernel().mpi_manager.get_num_processes();
  size_t current_buff_size_per_rank = kernel().mpi_manager.get_send_recv_count_spike_data_per_rank();

  if ( global_max_spikes_per_rank_ < send_recv_buffer_shrink_limit_ * current_buff_size_per_rank )
  {
    const size_t new_buff_size_per_rank = std::max( ( 2 * kernel().mpi_manager.get_num_source_groups_per_rank() + 1 ),
      static_cast< size_t >( ( 1 + send_recv_buffer_shrink_spare_ ) * global_max_spikes_per_rank_ ) );
    kernel().mpi_manager.set_send_recv_count_spike_data_per_rank( new_buff_size_per_rank );
    resize_send_recv_buffers_spike_data_();
    send_recv_buffer_resize_log_.add_entry( global_max_spikes_per_rank_, new_buff_size_per_rank );
  }

  /* The following do-while loop is executed
   * - once if all spikes fit into current send buffers on all ranks
   * - twice if send buffer size needs to be increased to fit in all spikes
   */
  bool all_spikes_transmitted;
  std::span< SpikeDataT > send_buffer_span = send_buffer;
  std::span< SpikeDataT > recv_buffer_span = recv_buffer;
  SendBufferPosition spike_register_position;
  size_t last_buffer_size = 0;
  send_buffer_segments_.clear();
  do
  {
#ifdef TIMER_DETAILED
    {
      sw_collocate_spike_data_.start();
    }
#endif

    // Collocate spikes to send buffer
    collocate_spike_data_buffers_( emitted_spikes_register, send_buffer_span, spike_register_position );

#ifdef TIMER_DETAILED
    {
      sw_collocate_spike_data_.stop();
      sw_communicate_spike_data_.start();
    }
#endif

    kernel().mpi_manager.communicate_spike_data_Alltoall( send_buffer_span, recv_buffer_span );

#ifdef TIMER_DETAILED
    {
      sw_communicate_spike_data_.stop();
    }
#endif

    // Determine the actual maximum number of spikes any rank sent to any other rank.
    const auto [ max_spikes_sent, num_remaining_spikes ] = get_global_spike_counter_( recv_buffer_span );
    global_max_spikes_per_rank_ += max_spikes_sent;

    // The current size of the send/recv buffer including metadata
    current_buff_size_per_rank = kernel().mpi_manager.get_send_recv_count_spike_data_per_rank();
    all_spikes_transmitted = num_remaining_spikes == 0;
    send_buffer_segments_.push_back( ( current_buff_size_per_rank - last_buffer_size ) * num_ranks );
    last_buffer_size = current_buff_size_per_rank;

    if ( not all_spikes_transmitted )
    {
      // New size is equal to the size we would have needed to immediately communicate all spikes and an additional
      // growing factor to account for increased activity in the future.
      const size_t metadata_size_per_rank = kernel().mpi_manager.get_num_source_groups_per_rank() + 1;
      const size_t buffer_size_for_next_segment = ( num_remaining_spikes + metadata_size_per_rank );
      const size_t new_buff_size_per_rank = static_cast< size_t >( ( 1 + send_recv_buffer_grow_extra_ ) * ( buffer_size_for_next_segment + current_buff_size_per_rank ) );

      kernel().mpi_manager.set_send_recv_count_spike_data_per_rank( new_buff_size_per_rank );
      send_recv_buffer_resize_log_.add_entry( num_remaining_spikes + current_buff_size_per_rank, new_buff_size_per_rank );

      // Resizing the send/recv buffers for spike data requires setting the span for the next communication round
      // accordingly.
      resize_send_recv_buffers_spike_data_();
      // Adjust our spans to only use the parts of the communication buffers that are not already filled with spikes.
      send_buffer_span = std::span( send_buffer ).subspan( current_buff_size_per_rank * num_ranks );
      recv_buffer_span = std::span( recv_buffer ).subspan( current_buff_size_per_rank * num_ranks );
    }

  } while ( not all_spikes_transmitted );

  // We cannot shrink buffers here, because they first need to be read out by deliver events. Shrinking will happen at
  // beginning of next gather.
}

template < typename SpikeDataT >
void
EventDeliveryManager::collocate_spike_data_buffers_(
  std::vector< std::vector< std::vector< SpikeDataT > >* >& emitted_spikes_register,
  std::span< SpikeDataT > send_buffer,
  SendBufferPosition& spike_register_position )
{
  const size_t source_groups_per_rank = kernel().mpi_manager.get_num_source_groups_per_rank();
  const size_t num_ranks = kernel().mpi_manager.get_num_processes();
  const size_t num_threads = kernel().vp_manager.get_num_threads();
  const size_t send_recv_count_per_rank = send_buffer.size() / num_ranks;
  size_t end_idx_offset = send_recv_count_per_rank - source_groups_per_rank - 1;
  size_t max_num_remaining_spikes = 0;
  size_t max_num_spikes_sent = 0;

  // Initialize meta-data for each rank
  for ( size_t rank = 0; rank != num_ranks; ++rank )
  {
    for ( size_t s = 0; s != source_groups_per_rank; ++s )
    {
      send_buffer[ end_idx_offset + send_recv_count_per_rank * rank + s ].set_num_spikes( 0 );
      send_buffer[ end_idx_offset + send_recv_count_per_rank * rank + s ].set_source_group_offset( 0 );
    }
  }

  // TODO JV: Parallel for (reduction: max of max_num_remaining_spikes)
  // Loop over target ranks
  for ( size_t rank = 0; rank != num_ranks; ++rank )
  {
    size_t max_num_remaining_spikes_for_rank = 0;
    // Calculate send buffer starting position
    auto send_buffer_it = send_buffer.begin() + send_recv_count_per_rank * rank;
    size_t send_buffer_remaining_space = send_recv_count_per_rank - source_groups_per_rank - 1;

    while ( spike_register_position.tid_[ rank ] != num_threads and send_buffer_remaining_space != 0 )
    {
      auto& emitted_spikes_for_rank = ( *( emitted_spikes_register[ spike_register_position.tid_[ rank ] ] ) )[ rank ];
      const size_t source_group_idx = thread_to_source_group_[ spike_register_position.tid_[ rank ] ];
      const auto spikes_register_begin_it = emitted_spikes_for_rank.begin() + spike_register_position.idx_[ rank ];

      // This is the maximum number of spikes we could send immediately
      size_t spikes_to_send = emitted_spikes_for_rank.size() - spike_register_position.idx_[ rank ];
      // The actual number of spikes to send depends on the remaining space
      if ( spikes_to_send <= send_buffer_remaining_space )
      {
        spike_register_position.idx_[ rank ] = 0;
        ++spike_register_position.tid_[ rank ];
        send_buffer_remaining_space -= spikes_to_send;
      }
      else
      {
        max_num_remaining_spikes += spikes_to_send - send_buffer_remaining_space;
        spike_register_position.idx_[ rank ] += send_buffer_remaining_space;
        spikes_to_send = send_buffer_remaining_space;
        send_buffer_remaining_space = 0;
      }
      const auto spikes_register_end_it = spikes_register_begin_it + spikes_to_send;

      std::copy( spikes_register_begin_it, spikes_register_end_it, send_buffer_it );
      send_buffer_it += spikes_to_send;
      max_num_spikes_sent += spikes_to_send;

      // Set the start position of spikes for source group, which equals the number of spikes sent for the previous
      // source group plus the respective start position of that source group. Also set the number of spikes to expect.
      const size_t new_num_spikes = send_buffer[ end_idx_offset + send_recv_count_per_rank * rank + source_group_idx ].get_num_spikes() + spikes_to_send;
      send_buffer[ end_idx_offset + send_recv_count_per_rank * rank + source_group_idx ].set_num_spikes( new_num_spikes );
      send_buffer[ end_idx_offset + send_recv_count_per_rank * rank + source_group_idx + 1 ].set_source_group_offset(
        send_buffer[ end_idx_offset + send_recv_count_per_rank * rank + source_group_idx ].get_source_group_offset() + new_num_spikes );

      if ( spike_register_position.idx_[ rank ] == 0 )
      {
        emitted_spikes_for_rank.clear();
      }
    }

    for ( size_t tid = spike_register_position.tid_[ rank ] + 1; tid < num_threads; ++tid )
    {
      max_num_remaining_spikes_for_rank += ( *( emitted_spikes_register[ tid ] ) )[ rank ].size();
    }
    max_num_remaining_spikes = std::max( max_num_remaining_spikes_for_rank, max_num_remaining_spikes );
  }

  // In the very last entry, set the maximum number of spikes sent and remaining spikes to be sent for any the ranks.
  for ( size_t rank = 0; rank != num_ranks; ++rank )
  {
    send_buffer[ send_recv_count_per_rank * ( rank + 1 ) - 1 ].set_source_group_offset( max_num_spikes_sent );
    send_buffer[ send_recv_count_per_rank * ( rank + 1 ) - 1 ].set_num_spikes( max_num_remaining_spikes );
  }
}

template < typename SpikeDataT >
std::pair< size_t, size_t >
EventDeliveryManager::get_global_spike_counter_( std::span< SpikeDataT >& recv_buffer ) const
{
  size_t maximum_sent = 0;
  size_t maximum_remaining = 0;
  const size_t num_ranks = kernel().mpi_manager.get_num_processes();
  const size_t send_recv_count_per_rank = recv_buffer.size() / num_ranks;
  for ( size_t source_rank = 0; source_rank < num_ranks; ++source_rank )
  {
    maximum_sent = std::max(
      recv_buffer[ send_recv_count_per_rank * ( source_rank + 1 ) - 1 ].get_source_group_offset(), maximum_sent );
    maximum_remaining =
      std::max( recv_buffer[ send_recv_count_per_rank * ( source_rank + 1 ) - 1 ].get_num_spikes(), maximum_remaining );
  }
  return { maximum_sent, maximum_remaining };
}

void
EventDeliveryManager::deliver_events( const size_t tid )
{
  if ( off_grid_spiking_ )
  {
    deliver_events_( tid, recv_buffer_off_grid_spike_data_ );
  }
  else
  {
    deliver_events_( tid, recv_buffer_spike_data_ );
  }
}

template < typename SpikeDataT >
void
EventDeliveryManager::deliver_events_( const size_t tid, const std::vector< SpikeDataT >& recv_buffer )
{
  // deliver only at beginning of time slice
  if ( kernel().simulation_manager.get_from_step() > 0 )
  {
    return;
  }

  // TODO JV: Make this task-based instead of parallel for

  const size_t source_groups_per_rank = kernel().mpi_manager.get_num_source_groups_per_rank();
  const size_t num_ranks = kernel().mpi_manager.get_num_processes();
  const std::vector< ConnectorModel* >& cm = kernel().model_manager.get_connection_models( tid );

  // prepare Time objects for every possible time stamp within min_delay_
  std::vector< Time > prepared_timestamps( kernel().connection_manager.get_min_delay() );
  for ( size_t lag = 0; lag < static_cast< size_t >( kernel().connection_manager.get_min_delay() ); ++lag )
  {
    // Subtract min_delay because spikes were emitted in previous time slice and we use current clock.
    prepared_timestamps[ lag ] =
      kernel().simulation_manager.get_clock() + Time::step( lag + 1 - kernel().connection_manager.get_min_delay() );
  }

  assert(send_buffer_segments_.size() > 0);
  auto recv_buffer_start_it = recv_buffer.cbegin();
  for ( size_t send_buffer_segment_num : send_buffer_segments_ )
  {
    const size_t send_recv_count_per_rank = send_buffer_segment_num / num_ranks;
    for ( std::pair< size_t, size_t >& source_group : responsible_source_groups_[ tid ] )
    {
      const size_t rank = source_group.first;
      const size_t rank_source_group_idx = source_group.second;
      auto recv_buffer_it = recv_buffer_start_it + send_recv_count_per_rank * rank;
      const size_t recv_buffer_rank_source_group_begin =
        ( recv_buffer_it + send_recv_count_per_rank - 1 - source_groups_per_rank + rank_source_group_idx )
          ->get_source_group_offset();
      // Get number of spikes to deliver for the source group which we are responsible for
      const size_t num_spikes_received =
        ( recv_buffer_it + send_recv_count_per_rank - 1 - source_groups_per_rank + rank_source_group_idx )
          ->get_num_spikes();
      assert( static_cast<long>(num_spikes_received) <= static_cast<long>(send_recv_count_per_rank) - 1 - static_cast<long>(source_groups_per_rank) );
      assert( static_cast<long>(recv_buffer_rank_source_group_begin) <= static_cast<long>(send_recv_count_per_rank) - 1 - static_cast<long>(source_groups_per_rank) );
      recv_buffer_it += recv_buffer_rank_source_group_begin;

      // For each batch, extract data first from receive buffer into value-specific arrays, then deliver from these
      // arrays
      constexpr size_t SPIKES_PER_BATCH = 8;
      const size_t num_batches = num_spikes_received / SPIKES_PER_BATCH;
      const size_t num_remaining_entries = num_spikes_received - num_batches * SPIKES_PER_BATCH;

      SpikeEvent se_batch[ SPIKES_PER_BATCH ];
      size_t syn_id_batch[ SPIKES_PER_BATCH ];
      size_t lcid_batch[ SPIKES_PER_BATCH ];

      for ( size_t i = 0; i < num_batches; ++i )
      {
        for ( size_t j = 0; j < SPIKES_PER_BATCH; ++j )
        {
          assert(recv_buffer.size() > (recv_buffer_it + i * SPIKES_PER_BATCH + j) - recv_buffer.cbegin());
          const auto& spike_data = recv_buffer_it + i * SPIKES_PER_BATCH + j;
          syn_id_batch[ j ] = spike_data->get_syn_id();
          lcid_batch[ j ] = spike_data->get_lcid();
          se_batch[ j ].set_stamp( prepared_timestamps[ spike_data->get_lag() ] );
          se_batch[ j ].set_offset( spike_data->get_offset() );
          se_batch[ j ].set_sender_node_id_info( syn_id_batch[ j ], lcid_batch[ j ] );
        }
        for ( size_t j = 0; j < SPIKES_PER_BATCH; ++j )
        {
          kernel().connection_manager.send( tid, syn_id_batch[ j ], lcid_batch[ j ], cm, se_batch[ j ] );
        }
      }

      // Processed all regular-sized batches, now do remainder
      for ( size_t j = 0; j < num_remaining_entries; ++j )
      {
        assert(recv_buffer.size() > (recv_buffer_it + num_batches * SPIKES_PER_BATCH + j) - recv_buffer.cbegin());
        const auto& spike_data = recv_buffer_it + num_batches * SPIKES_PER_BATCH + j;
        syn_id_batch[ j ] = spike_data->get_syn_id();
        lcid_batch[ j ] = spike_data->get_lcid();
        se_batch[ j ].set_stamp( prepared_timestamps[ spike_data->get_lag() ] );
        se_batch[ j ].set_offset( spike_data->get_offset() );
        se_batch[ j ].set_sender_node_id_info( syn_id_batch[ j ], lcid_batch[ j ] );
      }
      for ( size_t j = 0; j < num_remaining_entries; ++j )
      {
        kernel().connection_manager.send( tid, syn_id_batch[ j ], lcid_batch[ j ], cm, se_batch[ j ] );
      }
    }
    recv_buffer_start_it += send_buffer_segment_num;
  }
}

void
EventDeliveryManager::gather_target_data( const size_t tid )
{
  assert( not kernel().connection_manager.is_source_table_cleared() );

  // assume all threads have some work to do
  gather_completed_checker_.set_false( tid );
  assert( gather_completed_checker_.all_false() );

  const AssignedRanks assigned_ranks = kernel().vp_manager.get_assigned_ranks( tid );

  kernel().connection_manager.prepare_target_table( tid );
  kernel().connection_manager.reset_source_table_entry_point( tid );

  while ( gather_completed_checker_.any_false() )
  {
    // assume this is the last gather round and change to false
    // otherwise
    gather_completed_checker_.set_true( tid );

#pragma omp master
    {
      if ( kernel().mpi_manager.adaptive_target_buffers() and buffer_size_target_data_has_changed_ )
      {
        resize_send_recv_buffers_target_data();
      }
    } // of omp master; (no barrier)
#pragma omp barrier

    kernel().connection_manager.restore_source_table_entry_point( tid );

    TargetSendBufferPosition send_buffer_position(
      assigned_ranks, kernel().mpi_manager.get_send_recv_count_target_data_per_rank() );

    const bool gather_completed = collocate_target_data_buffers_( tid, assigned_ranks, send_buffer_position );
    gather_completed_checker_.logical_and( tid, gather_completed );

    if ( gather_completed_checker_.all_true() )
    {
      set_complete_marker_target_data_( assigned_ranks, send_buffer_position );
    }
    kernel().connection_manager.save_source_table_entry_point( tid );
#pragma omp barrier
    kernel().connection_manager.clean_source_table( tid );

#pragma omp master
    {
#ifdef TIMER_DETAILED
      sw_communicate_target_data_.start();
#endif
      kernel().mpi_manager.communicate_target_data_Alltoall(
        std::span( send_buffer_target_data_ ), std::span( recv_buffer_target_data_ ) );
#ifdef TIMER_DETAILED
      sw_communicate_target_data_.stop();
#endif
    } // of omp master (no barriers!)
#pragma omp barrier

    const bool distribute_completed = distribute_target_data_buffers_( tid );
    gather_completed_checker_.logical_and( tid, distribute_completed );

    // resize mpi buffers, if necessary and allowed
    if ( gather_completed_checker_.any_false() and kernel().mpi_manager.adaptive_target_buffers() )
    {
#pragma omp master
      {
        buffer_size_target_data_has_changed_ = kernel().mpi_manager.increase_buffer_size_target_data();
      }
#pragma omp barrier
    }
  } // of while

  kernel().connection_manager.clear_source_table( tid );
}

bool
EventDeliveryManager::collocate_target_data_buffers_( const size_t tid,
  const AssignedRanks& assigned_ranks,
  TargetSendBufferPosition& send_buffer_position )
{
  size_t source_rank;
  TargetData next_target_data;
  bool valid_next_target_data;
  bool is_source_table_read = true;

  // no ranks to process for this thread
  if ( assigned_ranks.begin == assigned_ranks.end )
  {
    kernel().connection_manager.no_targets_to_process( tid );
    return is_source_table_read;
  }

  // reset markers
  for ( size_t rank = assigned_ranks.begin; rank < assigned_ranks.end; ++rank )
  {
    // reset last entry to avoid accidentally communicating done
    // marker
    send_buffer_target_data_[ send_buffer_position.end( rank ) - 1 ].reset_marker();
    // set first entry to invalid to avoid accidentally reading
    // uninitialized parts of the receive buffer
    send_buffer_target_data_[ send_buffer_position.begin( rank ) ].set_invalid_marker();
  }

  while ( true )
  {
    valid_next_target_data = kernel().connection_manager.get_next_target_data(
      tid, assigned_ranks.begin, assigned_ranks.end, source_rank, next_target_data );
    if ( valid_next_target_data ) // add valid entry to MPI buffer
    {
      if ( send_buffer_position.is_chunk_filled( source_rank ) )
      {
        // entry does not fit in this part of the MPI buffer any more,
        // so we need to reject it
        kernel().connection_manager.reject_last_target_data( tid );
        // after rejecting the last target, we need to save the
        // position to start at this point again next communication
        // round
        kernel().connection_manager.save_source_table_entry_point( tid );
        // we have just rejected an entry, so source table can not be
        // fully read
        is_source_table_read = false;
        if ( send_buffer_position.are_all_chunks_filled() ) // buffer is full
        {
          return is_source_table_read;
        }
        else
        {
          continue;
        }
      }
      else
      {
        send_buffer_target_data_[ send_buffer_position.idx( source_rank ) ] = next_target_data;
        send_buffer_position.increase( source_rank );
      }
    }
    else // all connections have been processed
    {
      // mark end of valid data for each rank
      for ( size_t rank = assigned_ranks.begin; rank < assigned_ranks.end; ++rank )
      {
        if ( send_buffer_position.idx( rank ) > send_buffer_position.begin( rank ) )
        {
          send_buffer_target_data_[ send_buffer_position.idx( rank ) - 1 ].set_end_marker();
        }
        else
        {
          send_buffer_target_data_[ send_buffer_position.begin( rank ) ].set_invalid_marker();
        }
      }
      return is_source_table_read;
    } // of else
  } // of while(true)
}

void
nest::EventDeliveryManager::set_complete_marker_target_data_( const AssignedRanks& assigned_ranks,
  const TargetSendBufferPosition& send_buffer_position )
{
  for ( size_t rank = assigned_ranks.begin; rank < assigned_ranks.end; ++rank )
  {
    const size_t idx = send_buffer_position.end( rank ) - 1;
    send_buffer_target_data_[ idx ].set_complete_marker();
  }
}

bool
nest::EventDeliveryManager::distribute_target_data_buffers_( const size_t tid )
{
  bool are_others_completed = true;
  const unsigned int send_recv_count_target_data_per_rank =
    kernel().mpi_manager.get_send_recv_count_target_data_per_rank();

  for ( size_t rank = 0; rank < kernel().mpi_manager.get_num_processes(); ++rank )
  {
    // Check last entry for completed marker
    if ( not recv_buffer_target_data_[ ( rank + 1 ) * send_recv_count_target_data_per_rank - 1 ].is_complete_marker() )
    {
      are_others_completed = false;
    }

    // Were any targets sent by this rank?
    if ( recv_buffer_target_data_[ rank * send_recv_count_target_data_per_rank ].is_invalid_marker() )
    {
      continue;
    }

    for ( unsigned int i = 0; i < send_recv_count_target_data_per_rank; ++i )
    {
      const TargetData& target_data = recv_buffer_target_data_[ rank * send_recv_count_target_data_per_rank + i ];
      if ( target_data.get_source_tid() == tid )
      {
        kernel().connection_manager.add_target( tid, rank, target_data );
      }

      // Is this the last target from this rank?
      if ( target_data.is_end_marker() )
      {
        break;
      }
    }
  }

  return are_others_completed;
}

} // of namespace nest
