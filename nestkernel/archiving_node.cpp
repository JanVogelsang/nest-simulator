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
#include "connector_model.h"
#include "kernel_manager.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{

// member functions for ArchivingNode

ArchivingNode::ArchivingNode()
  : tau_minus_( 20.0 )
  , tau_minus_inv_( 1. / tau_minus_ )
  , tau_minus_triplet_( 110.0 )
  , tau_minus_triplet_inv_( 1. / tau_minus_triplet_ )
  , max_axonal_delay_( 0 )
  , max_dendritic_delay_( 0 )
{
}

ArchivingNode::ArchivingNode( const ArchivingNode& n )
  : StructuralPlasticityNode( n )
  , tau_minus_( n.tau_minus_ )
  , tau_minus_inv_( n.tau_minus_inv_ )
  , tau_minus_triplet_( n.tau_minus_triplet_ )
  , tau_minus_triplet_inv_( n.tau_minus_triplet_inv_ )
  , max_axonal_delay_( n.max_axonal_delay_ )
  , max_dendritic_delay_( n.max_dendritic_delay_ )
{
}

void ArchivingNode::init()
{
  Node::init();
  if ( max_axonal_delay_ )
  {
    intermediate_spike_buffer_.resize( max_axonal_delay_ );
  }
}

void
ArchivingNode::register_stdp_connection( const delay axonal_delay, const delay dendritic_delay, const synindex syn_id )
{
  max_axonal_delay_ = std::max( max_axonal_delay_, axonal_delay );
  max_dendritic_delay_ = std::max( max_dendritic_delay_, dendritic_delay );

  // check if synapse type has been registered already
  if( std::find(stdp_synapse_types_.begin(), stdp_synapse_types_.end(), syn_id ) == stdp_synapse_types_.end() ) {
    stdp_synapse_types_.push_back( syn_id );
  }
}

bool
ArchivingNode::has_stdp_connections() const
{
  return stdp_synapse_types_.size() > 0;
}

// double
// ArchivingNode::get_K_value( double t )
//{
//   // case when the neuron has not yet spiked
//   if ( history_.empty() )
//   {
//     trace_ = 0.;
//     return trace_;
//   }
//
//   // search for the latest post spike in the history buffer that came strictly before `t`
//   int i = history_.size() - 1;
//   while ( i >= 0 )
//   {
//     if ( t - history_[ i ].t > kernel().connection_manager.get_stdp_eps() )
//     {
//       trace_ = ( history_[ i ].Kminus * std::exp( ( history_[ i ].t - t ) * tau_minus_inv_ ) );
//       return trace_;
//     }
//     --i;
//   }
//
//   // this case occurs when the trace was requested at a time precisely at or
//   // before the first spike in the history
//   trace_ = 0.;
//   return trace_;
// }

// void
// ArchivingNode::get_K_values( double t, double& K_value, double& nearest_neighbor_K_value, double& K_triplet_value )
//{
//   // case when the neuron has not yet spiked
//   if ( history_.empty() )
//   {
//     K_triplet_value = Kminus_triplet_;
//     nearest_neighbor_K_value = Kminus_;
//     K_value = Kminus_;
//     return;
//   }
//
//   // search for the latest post spike in the history buffer that came strictly
//   // before `t`
//   int i = history_.size() - 1;
//   while ( i >= 0 )
//   {
//     if ( t - history_[ i ].t > kernel().connection_manager.get_stdp_eps() )
//     {
//       K_triplet_value =
//         ( history_[ i ].Kminus_triplet * std::exp( ( history_[ i ].t - t ) * tau_minus_triplet_inv_ ) );
//       K_value = ( history_[ i ].Kminus * std::exp( ( history_[ i ].t - t ) * tau_minus_inv_ ) );
//       nearest_neighbor_K_value = std::exp( ( history_[ i ].t - t ) * tau_minus_inv_ );
//       return;
//     }
//     --i;
//   }
//
//   // this case occurs when the trace was requested at a time precisely at or
//   // before the first spike in the history
//   K_triplet_value = 0.0;
//   nearest_neighbor_K_value = 0.0;
//   K_value = 0.0;
// }

// void
// ArchivingNode::get_history( double t1,
//   double t2,
//   std::deque< ArchivedSpikeTrace >::iterator* start,
//   std::deque< ArchivedSpikeTrace >::iterator* finish )
//{
//   *finish = history_.end();
//   if ( history_.empty() )
//   {
//     *start = *finish;
//     return;
//   }
//   std::deque< ArchivedSpikeTrace >::reverse_iterator runner = history_.rbegin();
//   const double t2_lim = t2 + kernel().connection_manager.get_stdp_eps();
//   const double t1_lim = t1 + kernel().connection_manager.get_stdp_eps();
//   while ( runner != history_.rend() and runner->t >= t2_lim )
//   {
//     ++runner;
//   }
//   *finish = runner.base();
//   while ( runner != history_.rend() and runner->t >= t1_lim )
//   {
//     runner->access_counter_++;
//     ++runner;
//   }
//   *start = runner.base();
// }

void
ArchivingNode::update_stdp_connections( const delay lag )
{
  if ( max_axonal_delay_ )
  {
    const std::vector< ConnectorModel* >& cm = kernel().model_manager.get_connection_models( get_thread() );
    auto [ current_pre_synaptic_spike, last_pre_synaptic_spike ] = intermediate_spike_buffer_.get_next_spikes();
    // Process all pending pre-synaptic spikes. Lag+1 is required, as lag marks the beginning of the current update step and
    // the update should be performed at the end instead.
    while ( current_pre_synaptic_spike != last_pre_synaptic_spike and current_pre_synaptic_spike->t_syn_lag == lag + 1 )
    {
      const double eps = kernel().connection_manager.get_stdp_eps();  // only get eps when there is a pre-synaptic spike
      const synindex syn_id = current_pre_synaptic_spike->syn_id;
      const index local_connection_id = current_pre_synaptic_spike->local_connection_id;
      const Time pre_spike_time = current_pre_synaptic_spike->t_stamp;

      const delay dendritic_delay = connections_[ syn_id ]->get_dendritic_delay( local_connection_id );
      const delay axonal_delay = current_pre_synaptic_spike->axonal_delay;
      const double last_pre_spike_time_ms = connections_[ syn_id ]->get_last_presynaptic_spike( local_connection_id );
      // If a pre-synaptic spike is about to be processed, make sure to process all post-synaptic spikes first, which
      // are due before or at the same time of the pre-synaptic spike.
      for ( std::deque< double >::const_iterator archived_spike_it = history_.cbegin(); archived_spike_it != history_.cend(); ++archived_spike_it )
      {
        const double post_spike_time_ms = *archived_spike_it + Time::delay_steps_to_ms( dendritic_delay );
        // Skip post-synaptic spikes which have been processed by the dendritic delay region already after the last
        // communication round.
        const double last_communication_time_ms = pre_spike_time.get_ms() + Time::delay_steps_to_ms( axonal_delay - lag - 1 );
        if ( last_communication_time_ms > post_spike_time_ms + eps or std::abs( last_communication_time_ms - post_spike_time_ms ) < eps )
        {
          continue;
        }

        // Only process spikes which reached the synapse after the last (already processed) pre-synaptic spike.
        if ( last_pre_spike_time_ms > post_spike_time_ms + eps or std::abs( last_pre_spike_time_ms - post_spike_time_ms ) < eps )
        {
          continue;
        }

        // TODO JV (pt): Check if this works for precise spikes
        // Don't process spikes which didn't reach the synapse yet.
        if ( post_spike_time_ms > pre_spike_time.get_ms() + Time::delay_steps_to_ms( axonal_delay ) + eps )
        {
          break;
        }

        // TODO JV (pt): Weight recorder event
        connections_[ syn_id ]->process_post_synaptic_spike( local_connection_id, post_spike_time_ms, cm[ syn_id ] );
      }

      Node::deliver_event( syn_id,
        local_connection_id,
        cm,
        pre_spike_time,
        axonal_delay,
        dendritic_delay,
        0 ); // TODO JV (pt): Precise spikes
      ++current_pre_synaptic_spike;
    }
  }
}

void ArchivingNode::prepare_update()
{
  intermediate_spike_buffer_.prepare_next_slice();

  const std::vector< ConnectorModel* >& cm = kernel().model_manager.get_connection_models( get_thread() );
  // Process all pre- and post-synaptic spikes in relative order. Processes all pre-synaptic spikes that were just
  // communicated and all post-synaptic spikes in the current archive.
  auto [ current_pre_synaptic_spike, last_pre_synaptic_spike ] = intermediate_spike_buffer_.get_next_spikes();
  const Time origin = kernel().simulation_manager.get_slice_origin();
  const delay min_delay = kernel().connection_manager.get_min_delay();
  for ( delay lag = 1; lag != min_delay + 1; ++lag)
  {
    const Time t_now = origin + Time::step( lag );
    std::deque< double >::const_iterator archived_spike_it = history_.cbegin();
    while ( archived_spike_it != history_.cend() )
    {
      // TODO JV (pt): This needs some more debugging, as the time module seems to do mysterious things with floats
      //  sometimes which leads to wrong values here (seems to be due to rounding).
      const tic_t time_since_post_spike_tics = ( t_now - Time::ms( *archived_spike_it ) ).get_tics();
      if ( time_since_post_spike_tics >= 0 )
      {
        // TODO JV (pt): Precise spike times will lead to a wrong delay here, as Time::ms does a round while a ceil would
        //  be required to "remove" the offset from the spike time.
        const delay dendritic_delay = Time( Time::tic( time_since_post_spike_tics ) ).get_steps();
        for ( const synindex stdp_syn_id : stdp_synapse_types_ )
        {
          // The trace for the previous timestep is updated here, as it is important to make sure the trace didn't get
          // updated before any pre-synaptic spikes were processed.
          connections_[ stdp_syn_id ]->update_trace( *archived_spike_it, dendritic_delay - 1, tau_minus_inv_ );
        }
        // as soon as dendritic delay > max dendritic delay here, the spike can be safely removed from history
        if ( dendritic_delay > max_dendritic_delay_ )
        {
          // This is guaranteed to always start erasing from the beginning of the history as entries are inserted into the
          // history in chronological order. Thus, this will never create any holes.
          archived_spike_it = history_.erase( archived_spike_it );  // it points to next element after erasing
          continue;
        }
        // process all pending post-synaptic spikes
        for ( const synindex stdp_syn_id : stdp_synapse_types_ )
        {
          connections_[ stdp_syn_id ]->update_stdp_connections( t_now.get_ms(), dendritic_delay, cm[ stdp_syn_id ] );
        }
      }
      ++archived_spike_it;
    }

    // find the next pre-synaptic spike for the current lag
    while ( current_pre_synaptic_spike != last_pre_synaptic_spike and current_pre_synaptic_spike->t_syn_lag == lag )
    {
      Node::deliver_event( current_pre_synaptic_spike->syn_id,
        current_pre_synaptic_spike->local_connection_id,
        cm,
        current_pre_synaptic_spike->t_stamp,
        current_pre_synaptic_spike->axonal_delay,
        0 ); // TODO JV (pt): Precise spikes
      ++current_pre_synaptic_spike;
    }
  }

  // prepare the next update cycle
  intermediate_spike_buffer_.clean_slice();
  intermediate_spike_buffer_.increase_slice_index();
  intermediate_spike_buffer_.prepare_next_slice();
}

void ArchivingNode::end_update()
{
  if ( max_axonal_delay_ > 0 )  // no need to clean the intermediate spike buffer if there is no axonal delay
  {
    intermediate_spike_buffer_.clean_slice();
  }
}

void
ArchivingNode::deliver_event( const synindex syn_id,
  const index local_target_connection_id,
  const std::vector< ConnectorModel* >& cm,
  const Time lag,
  const delay axonal_delay,
  const double offset )
{
  const delay t_syn = lag.get_steps() + axonal_delay;
  const delay min_delay = kernel().connection_manager.get_min_delay();
  const delay t_now = kernel().simulation_manager.get_slice_origin().get_steps() + min_delay;
  const delay time_until_reaching_synapse = t_syn - t_now;

  // STDP synapses need to make sure all post-synaptic spikes are known when delivering the spike to the synapse.
  // If we need to postpone the delivery, we need to check for how long. Furthermore, even if it was safe to deliver the
  // spike already, the weight would be incorrect if read out at the next possible time (i.e., after end of delivery).
  // The earliest possible time a spike can be delivered is right before the start of the next update cycle in which the
  // spike would reach the synapse.
  // However, dendritic delays smaller than min_delay add even more complexity here. If the dendritic delay is smaller
  // than the remaining time in the target slice (the "lag"), there might be a future post-synaptic slice that needs to
  // be processed before the pre-synaptic one, not known before the following update cycle. The pre-synaptic spike might
  // have to be delivered in that update cycle however, when the dendritic delay is small enough. Thus, these spikes
  // have to be handled during the actual update cycle instead of before.To handle all postponed spikes in the same
  // way, all spikes will therefore be processed during the update of the neuron at the correct times.
  // If a spike reached the synapse already at the current point in time, i.e. the end (inclusive) of the current slice,
  // the spike has to be processed instantly, but in the correct order relative to any possible post-synaptic spikes
  // which reached the synapse before the pre-synaptic spike. These spikes will be temporarily stored in the first slot
  // of the intermediate buffer which was previously occupied by the spikes of the last update cycle.


  // end of slice is inclusive, therefore subtract 1
  const unsigned long slices_to_postpone = static_cast< unsigned long >( ( time_until_reaching_synapse + min_delay - 1 ) / min_delay );
  const delay t_syn_lag = time_until_reaching_synapse + min_delay - slices_to_postpone * min_delay;
  intermediate_spike_buffer_.push_back( slices_to_postpone, lag, axonal_delay, syn_id, local_target_connection_id, t_syn_lag );
}

void
ArchivingNode::set_spiketime( Time const& t_sp, double offset )
{
  StructuralPlasticityNode::set_spiketime( t_sp, offset );

  const double t_sp_ms = t_sp.get_ms() - offset;

  if ( stdp_synapse_types_.size() > 0 ) // check if node has any incoming STDP connections
  {
    // update spiking history
    //    Kminus_ = Kminus_ * std::exp( ( last_spike_ - t_sp_ms ) * tau_minus_inv_ ) + 1.0;
    //    Kminus_triplet_ = Kminus_triplet_ * std::exp( ( last_spike_ - t_sp_ms ) * tau_minus_triplet_inv_ ) + 1.0;
    //    history_.push_back( ArchivedSpikeTrace( t_sp_ms, Kminus_, Kminus_triplet_ ) );
    history_.push_back( t_sp_ms );
  }
}

void
ArchivingNode::get_status( DictionaryDatum& d ) const
{
  // def< double >( d, names::t_spike, get_spiketime_ms() );
  def< double >( d, names::tau_minus, tau_minus_ );
  def< double >( d, names::tau_minus_triplet, tau_minus_triplet_ );
//  def< double >( d, names::post_trace, trace_ );
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
  history_.clear();
  for ( const synindex stdp_syn_id : stdp_synapse_types_ )
  {
    connections_[ stdp_syn_id ]->clear_history();
  }
}

} // of namespace nest
