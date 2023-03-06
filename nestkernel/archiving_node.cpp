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
  , max_dendritic_delay_( 0 )
{
}

ArchivingNode::ArchivingNode( const ArchivingNode& n )
  : StructuralPlasticityNode( n )
  , tau_minus_( n.tau_minus_ )
  , tau_minus_inv_( n.tau_minus_inv_ )
  , tau_minus_triplet_( n.tau_minus_triplet_ )
  , tau_minus_triplet_inv_( n.tau_minus_triplet_inv_ )
  , max_dendritic_delay_( n.max_dendritic_delay_ )
{
}

void
ArchivingNode::register_stdp_connection( const double dendritic_delay, const synindex syn_id )
{
  max_dendritic_delay_ = std::max( max_dendritic_delay_, dendritic_delay );
  stdp_synapse_types_.insert( syn_id );
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
ArchivingNode::update_stdp_connections( Time const& origin, const long from, const long to )
{
  const std::vector< ConnectorModel* >& cm = kernel().model_manager.get_connection_models( get_thread() );
  for ( long lag = from + 1; lag <= to; ++lag )
  {
    // let STDP connections process previous spikes of this neuron
    for ( auto archived_spike_it = history_.cbegin(); archived_spike_it < history_.cend(); ++archived_spike_it )
    {
      const double dendritic_delay = Time::delay_steps_to_ms( origin.get_steps() + lag ) - *archived_spike_it;
      if ( dendritic_delay < 0 )
      {
        continue;
      }

      for ( const synindex stdp_syn_id : stdp_synapse_types_ )
      {
        connections_[ stdp_syn_id ]->update_stdp_connections( *archived_spike_it, dendritic_delay, tau_minus_inv_, cm );
      }
      if ( dendritic_delay > max_dendritic_delay_
        or std::abs( dendritic_delay - max_dendritic_delay_ ) < kernel().connection_manager.get_stdp_eps() )
      {
        // This is guaranteed to always start erasing from the beginning of the history as entries are inserted into the
        // history in chronological order. Thus, this will never create any holes.
        history_.erase( archived_spike_it );
      }
    }
  }
}

void
ArchivingNode::deliver_event( const thread tid,
  const synindex syn_id,
  const index local_target_connection_id,
  const std::vector< ConnectorModel* >& cm,
  SpikeEvent& se )
{
  ConnectorBase* conn = connections_[ syn_id ];

  const double axonal_delay = conn->get_axonal_delay( local_target_connection_id );
  const double t_now = kernel().simulation_manager.get_slice_origin().get_ms();

  // STDP synapses need to make sure all post-synaptic spikes are known when delivering the spike to the synapse.
  // Spikes will therefore be stored in an intermediate spike buffer until no more post-synaptic spike could reach the
  // synapse before this spike will. 
  if ( axonal_delay > kernel().connection_manager.get_min_delay() )
  {

    intermediate_spike_buffer_.push_back(  );
  }
  else
  {
    process_event( tid, syn_id, local_target_connection_id, cm, se );
  }
}

void
ArchivingNode::process_event( const thread tid,
  const synindex syn_id,
  const index local_target_connection_id,
  const std::vector< ConnectorModel* >& cm,
  SpikeEvent& se )
{
  // TODO JV (pt): Think about removing access to connections_ in derived classes of node
  ConnectorBase* conn = connections_[ syn_id ];

  // Send the event to the connection over which this event is transmitted to the node. The connection modifies the
  // event by adding a weight.
  conn->send( tid, local_target_connection_id, cm, se, this );

  handle( se );
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
