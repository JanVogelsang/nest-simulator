/*
 *  archiving_node.cpp
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

#include "archiving_node.h"

// Includes from nestkernel:
#include "kernel_manager.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{

// member functions for ArchivingNode

ArchivingNode::ArchivingNode()
  : n_incoming_( 0 )
  , Kminus_( 0.0 )
  , Kminus_triplet_( 0.0 )
  , tau_minus_( 20.0 )
  , tau_minus_inv_( 1. / tau_minus_ )
  , tau_minus_triplet_( 110.0 )
  , tau_minus_triplet_inv_( 1. / tau_minus_triplet_ )
  , max_delay_( 0 )
  , trace_( 0.0 )
  , last_spike_( -1.0 )
  , has_stdp_ax_delay_( false )
{
  const size_t num_time_slots =
    kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
  correction_entries_stdp_ax_delay_.resize( num_time_slots );
}

ArchivingNode::ArchivingNode( const ArchivingNode& n )
  : StructuralPlasticityNode( n )
  , n_incoming_( n.n_incoming_ )
  , Kminus_( n.Kminus_ )
  , Kminus_triplet_( n.Kminus_triplet_ )
  , tau_minus_( n.tau_minus_ )
  , tau_minus_inv_( n.tau_minus_inv_ )
  , tau_minus_triplet_( n.tau_minus_triplet_ )
  , tau_minus_triplet_inv_( n.tau_minus_triplet_inv_ )
  , max_delay_( n.max_delay_ )
  , trace_( n.trace_ )
  , last_spike_( n.last_spike_ )
  , has_stdp_ax_delay_( false )
{
  const size_t num_time_slots =
    kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
  correction_entries_stdp_ax_delay_.resize( num_time_slots );
}

void
ArchivingNode::pre_run_hook_()
{
  const size_t num_time_slots =
    kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
  if ( correction_entries_stdp_ax_delay_.size() != num_time_slots )
  {
    correction_entries_stdp_ax_delay_.resize( num_time_slots );
  }
}

void
ArchivingNode::register_stdp_connection( double t_first_read, double delay )
{
  // Mark all entries in the deque, which we will not read in future as read by
  // this input, so that we safely increment the incoming number of
  // connections afterwards without leaving spikes in the history.
  // For details see bug #218. MH 08-04-22

  // JV: Remove this block
  for ( std::deque< histentry >::iterator runner = history_.begin();
        runner != history_.end() and ( t_first_read - runner->t_ > -1.0 * kernel().connection_manager.get_stdp_eps() );
        ++runner )
  {
    ( runner->access_counter_ )++;
  }

  n_incoming_++;

  max_delay_ = std::max( delay, max_delay_ );
}

double
ArchivingNode::get_K_value( double t )
{
  // case when the neuron has not yet spiked
  if ( history_.empty() )
  {
    trace_ = 0.;
    return trace_;
  }

  // search for the latest post spike in the history buffer that came strictly
  // before `t`
  int i = history_.size() - 1;
  while ( i >= 0 )
  {
    if ( t - history_[ i ].t_ > kernel().connection_manager.get_stdp_eps() )
    {
      trace_ = ( history_[ i ].Kminus_ * std::exp( ( history_[ i ].t_ - t ) * tau_minus_inv_ ) );
      return trace_;
    }
    --i;
  }

  // this case occurs when the trace was requested at a time precisely at or
  // before the first spike in the history
  trace_ = 0.;
  return trace_;
}

void
ArchivingNode::get_K_values( double t, double& K_value, double& nearest_neighbor_K_value, double& K_triplet_value )
{
  // case when the neuron has not yet spiked
  if ( history_.empty() )
  {
    K_triplet_value = Kminus_triplet_;
    nearest_neighbor_K_value = Kminus_;
    K_value = Kminus_;
    return;
  }

  // search for the latest post spike in the history buffer that came strictly
  // before `t`
  int i = history_.size() - 1;
  while ( i >= 0 )
  {
    if ( t - history_[ i ].t_ > kernel().connection_manager.get_stdp_eps() )
    {
      K_triplet_value =
        ( history_[ i ].Kminus_triplet_ * std::exp( ( history_[ i ].t_ - t ) * tau_minus_triplet_inv_ ) );
      K_value = ( history_[ i ].Kminus_ * std::exp( ( history_[ i ].t_ - t ) * tau_minus_inv_ ) );
      nearest_neighbor_K_value = std::exp( ( history_[ i ].t_ - t ) * tau_minus_inv_ );
      return;
    }
    --i;
  }

  // this case occurs when the trace was requested at a time precisely at or
  // before the first spike in the history
  K_triplet_value = 0.0;
  nearest_neighbor_K_value = 0.0;
  K_value = 0.0;
}

// JV: Remove this method
void
ArchivingNode::get_history( double t1,
  double t2,
  std::deque< histentry >::iterator* start,
  std::deque< histentry >::iterator* finish )
{
  *finish = history_.end();
  if ( history_.empty() )
  {
    *start = *finish;
    return;
  }
  std::deque< histentry >::reverse_iterator runner = history_.rbegin();
  const double t2_lim = t2 + kernel().connection_manager.get_stdp_eps();
  const double t1_lim = t1 + kernel().connection_manager.get_stdp_eps();
  while ( runner != history_.rend() and runner->t_ >= t2_lim )
  {
    ++runner;
  }
  *finish = runner.base();
  while ( runner != history_.rend() and runner->t_ >= t1_lim )
  {
    runner->access_counter_++;
    ++runner;
  }
  *start = runner.base();
}

void
ArchivingNode::deliver_event( const thread tid,
  const synindex syn_id,
  const index local_target_connection_id,
  const std::vector< ConnectorModel* >& cm,
  SpikeEvent& se )
{
  ConnectorBase* conn = connections_[ syn_id ];

  // STDP synapses need to make sure all post-synaptic spikes are known when delivering the spike to the synapse.
  // Spikes will therefore be stored in an intermediate spike buffer until no more post-synaptic spike could reach the
  // synapse before this spike will.

  // Only specific synapse types need to postpone the delivery
  /*if ( cm[ syn_id ]->requires_postponed_delivery() ) {
    if ( const double t_spike = se.get_stamp().get_ms(), delay dendritic_delay = 0, const Time& ori =
  kernel().simulation_manager.get_slice_origin(); t_spike + axonal - conn->get_connection_delay() < now ) {
      dynamic_spike_buffer_.push_back(se);
      return;
    }
  }*/

  // Send the event to the connection over which this event is transmitted to the node. The connection modifies the
  // event by adding a weight (TODO JV: only weight?).
  conn->send( tid, sources_[ syn_id ][ local_target_connection_id ].get_node_id(), local_target_connection_id, cm, se );

  handle( se );
}

void
ArchivingNode::set_spiketime( Time const& t_sp, double offset )
{
  StructuralPlasticityNode::set_spiketime( t_sp, offset );

  const double t_sp_ms = t_sp.get_ms() - offset;

  if ( n_incoming_ ) // JV: Remove this whole block, but still update the K value
  {
    // prune all spikes from history which are no longer needed
    // only remove a spike if:
    // - its access counter indicates it has been read out by all connected
    //   STDP synapses, and
    // - there is another, later spike, that is strictly more than
    //   (max_delay_ + eps) away from the new spike (at t_sp_ms)
    // [4.6, 5.4] -> 5.9, 6.3  --- 6.4+eps
    while ( history_.size() > 1 )
    {
      const double next_t_sp = history_[ 1 ].t_;
      if ( history_.front().access_counter_ >= n_incoming_
        and t_sp_ms - next_t_sp > 2 * max_delay_ + kernel().connection_manager.get_stdp_eps() )
      {
        history_.pop_front();
      }
      else
      {
        break;
      }
    }
    // update spiking history
    Kminus_ *= std::exp( ( last_spike_ - t_sp_ms ) * tau_minus_inv_ ) + 1.0;
    Kminus_triplet_ *= std::exp( ( last_spike_ - t_sp_ms ) * tau_minus_triplet_inv_ ) + 1.0;
    history_.push_back( histentry( t_sp_ms, Kminus_, Kminus_triplet_, 0 ) );
  }

  correct_synapses_stdp_ax_delay_( t_sp );
}

void
ArchivingNode::get_status( DictionaryDatum& d ) const
{
  def< double >( d, names::t_spike, get_spiketime_ms() );
  def< double >( d, names::tau_minus, tau_minus_ );
  def< double >( d, names::tau_minus_triplet, tau_minus_triplet_ );
  def< double >( d, names::post_trace, trace_ );
#ifdef DEBUG_ARCHIVER
  def< int >( d, names::archiver_length, history_.size() );
#endif

  // add status dict items from the parent class
  StructuralPlasticityNode::get_status( d );
}

void
ArchivingNode::set_status( const DictionaryDatum& d )
{
  // We need to preserve values in case invalid values are set
  double new_tau_minus = tau_minus_;
  double new_tau_minus_triplet = tau_minus_triplet_;
  updateValue< double >( d, names::tau_minus, new_tau_minus );
  updateValue< double >( d, names::tau_minus_triplet, new_tau_minus_triplet );

  if ( new_tau_minus <= 0.0 or new_tau_minus_triplet <= 0.0 )
  {
    throw BadProperty( "All time constants must be strictly positive." );
  }

  StructuralPlasticityNode::set_status( d );

  // do the actual update
  tau_minus_ = new_tau_minus;
  tau_minus_triplet_ = new_tau_minus_triplet;
  tau_minus_inv_ = 1. / tau_minus_;
  tau_minus_triplet_inv_ = 1. / tau_minus_triplet_;

  // check, if to clear spike history and K_minus
  bool clear = false;
  updateValue< bool >( d, names::clear, clear );
  if ( clear )
  {
    clear_history();
  }
}

void
ArchivingNode::clear_history()
{
  last_spike_ = -1.0;
  Kminus_ = 0.0;
  Kminus_triplet_ = 0.0;
  history_.clear();
}

void
ArchivingNode::add_correction_entry_stdp_ax_delay( SpikeEvent& spike_event,
  const double t_last_pre_spike,
  const double weight_revert,
  const double dendritic_delay )
{
  if ( not has_stdp_ax_delay_ )
  {
    has_stdp_ax_delay_ = true;

    const size_t num_time_slots =
      kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay();
    if ( correction_entries_stdp_ax_delay_.size() != num_time_slots )
    {
      correction_entries_stdp_ax_delay_.resize( num_time_slots );
    }
  }

  assert( correction_entries_stdp_ax_delay_.size()
    == static_cast< size_t >(
      kernel().connection_manager.get_min_delay() + kernel().connection_manager.get_max_delay() ) );

  const index idx = kernel().event_delivery_manager.get_modulo(
    spike_event.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() )
    - 2 * Time::delay_ms_to_steps( dendritic_delay ) );
  assert( static_cast< size_t >( idx ) < correction_entries_stdp_ax_delay_.size() );

  correction_entries_stdp_ax_delay_[ idx ].push_back(
    CorrectionEntrySTDPAxDelay( spike_event.get_sender_spike_data(), t_last_pre_spike, weight_revert ) );
}

void
ArchivingNode::reset_correction_entries_stdp_ax_delay_()
{
  if ( has_stdp_ax_delay_ )
  {
    const long mindelay_steps = kernel().connection_manager.get_min_delay();
    assert( correction_entries_stdp_ax_delay_.size()
      == static_cast< size_t >( mindelay_steps + kernel().connection_manager.get_max_delay() ) );

    for ( long lag = 0; lag < mindelay_steps; ++lag )
    {
      const long idx = kernel().event_delivery_manager.get_modulo( lag );
      assert( static_cast< size_t >( idx ) < correction_entries_stdp_ax_delay_.size() );

      correction_entries_stdp_ax_delay_[ idx ].clear();
    }
  }
}

void
ArchivingNode::correct_synapses_stdp_ax_delay_( const Time& t_spike )
{
  if ( has_stdp_ax_delay_ )
  {
    const Time& ori = kernel().simulation_manager.get_slice_origin();
    const Time& t_spike_rel = t_spike - ori;
    const long maxdelay_steps = kernel().connection_manager.get_max_delay();
    assert( correction_entries_stdp_ax_delay_.size()
      == static_cast< size_t >( kernel().connection_manager.get_min_delay() + maxdelay_steps ) );

    for ( long lag = t_spike_rel.get_steps() - 1; lag < maxdelay_steps + 1; ++lag )
    {
      const long idx = kernel().event_delivery_manager.get_modulo( lag );
      assert( static_cast< size_t >( idx ) < correction_entries_stdp_ax_delay_.size() );

      for ( auto it_corr_entry = correction_entries_stdp_ax_delay_[ idx ].begin();
            it_corr_entry < correction_entries_stdp_ax_delay_[ idx ].end();
            ++it_corr_entry )
      {
        connections_[ it_corr_entry->spike_data_.get_syn_id() ]
          ->correct_synapse_stdp_ax_delay( it_corr_entry->spike_data_,
            it_corr_entry->t_last_pre_spike_,
            &it_corr_entry->weight_revert_,
            t_spike.get_ms() );
      }
      // indicate that the new spike was processed by these STDP synapses
      history_.back().access_counter_ += correction_entries_stdp_ax_delay_[ idx ].size();
    }
  }
}

} // of namespace nest
