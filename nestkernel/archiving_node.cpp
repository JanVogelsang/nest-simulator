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

void
ArchivingNode::init()
{
  Node::init();
  intermediate_spike_buffer_.resize( max_axonal_delay_ );
}

void
ArchivingNode::register_stdp_connection( const delay axonal_delay, const delay dendritic_delay, const synindex syn_id )
{
  max_axonal_delay_ = std::max( max_axonal_delay_, axonal_delay );
  max_dendritic_delay_ = std::max( max_dendritic_delay_, dendritic_delay );

  // check if synapse type has been registered already
  if ( std::find( stdp_synapse_types_.begin(), stdp_synapse_types_.end(), syn_id ) == stdp_synapse_types_.end() )
  {
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
  if ( max_axonal_delay_ > 0 )
  {
#ifdef TIMER_DETAILED
    if ( get_thread() == 0 )
      kernel().event_delivery_manager.sw_stdp_delivery_.start();
#endif
    const std::vector< ConnectorModel* >& cm = kernel().model_manager.get_connection_models( get_thread() );
    auto [ current_pre_synaptic_spike, last_pre_synaptic_spike ] = intermediate_spike_buffer_.get_next_spikes();
    // Process all pending pre-synaptic spikes. Lag+1 is required, as lag marks the beginning of the current update step
    // and the update should be performed at the end instead.
    while ( current_pre_synaptic_spike != last_pre_synaptic_spike and current_pre_synaptic_spike->t_syn_lag == lag + 1 )
    {
      const synindex syn_id = current_pre_synaptic_spike->syn_id;
      const index local_connection_id = current_pre_synaptic_spike->local_connection_id;
      const size_t dendritic_delay_id = current_pre_synaptic_spike->dendritic_delay_id;
      const Time pre_spike_time = current_pre_synaptic_spike->t_stamp;

      const delay axonal_delay = current_pre_synaptic_spike->axonal_delay;
      // Time of the last communication round, which can be easily calculated from the pre-synaptic spikes arrival time
      // at the synapse using the current lag (instead of querying the simulation manager).
      const double last_communication_time_steps = pre_spike_time.get_steps() + axonal_delay - lag - 1;
      // TODO JV (pt): Precise spikes
      process_spikes_until_pre_synaptic_spike( syn_id,
        axonal_delay,
        local_connection_id,
        dendritic_delay_id,
        pre_spike_time,
        last_communication_time_steps,
        0,
        cm[ syn_id ] );
      ++current_pre_synaptic_spike;
    }
#ifdef TIMER_DETAILED
    if ( get_thread() == 0 )
      kernel().event_delivery_manager.sw_stdp_delivery_.stop();
#endif
  }
}

void
ArchivingNode::deliver_event_with_trace( const synindex syn_id,
  const index local_target_connection_id,
  const size_t dendritic_delay_id,
  const ConnectorModel* cm,
  const Time lag,
  const delay axonal_delay,
  const double offset )
{
  SpikeEvent se;
  se.set_stamp( lag + Time::step( axonal_delay ) );
  se.set_offset( offset ); // TODO JV (help): Why can't offset be incorporated into lag?
  se.set_sender_node_id_info( get_thread(), syn_id, get_node_id(), local_target_connection_id );

  // Send the event to the connection over which this event is transmitted to the node. The connection modifies the
  // event by adding a weight and optionally updates its internal state as well.
  connections_[ syn_id ]->send(
    get_thread(), get_node_id(), local_target_connection_id, dendritic_delay_id, tau_minus_inv_, history_, cm, se );

  // TODO JV (pt): Optionally, the rport can be set here (somehow). For example by just handing it as a parameter to
  //  handle, or just handing the entire local connection id to the handle function (and storing an array of rports
  //  which can be indexed by the local connection id).

  handle( se );
}

void
ArchivingNode::set_spiketime( const Time& t_sp, double offset )
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
