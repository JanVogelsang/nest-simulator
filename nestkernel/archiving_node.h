/*
 *  archiving_node.h
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

#ifndef ARCHIVING_NODE_H
#define ARCHIVING_NODE_H

// C++ includes:
#include <algorithm>
#include <deque>

// Includes from nestkernel:
#include "archived_spike.h"
#include "connector_base_impl.h"
#include "nest_time.h"
#include "nest_types.h"
#include "node.h"
#include "structural_plasticity_node.h"

// Includes from sli:
#include "dictdatum.h"

#ifdef HAVE_BOOST
// TODO JV (pt): This is a fixed-size container now, but this will not work with offgrid spikes, as they can't be
//  combined into a single slot for a simulation step and in theory multiple spikes could be emitted in a single
//  simulation step. There is "circular_buffer_space_optimized" which is resizable, but is not iterable anymore.
//  Maybe circular buffer is not the correct choice in general here.
#include <boost/circular_buffer.hpp>
template < typename T >
using SpikeArchive = boost::circular_buffer< T, std::allocator< T > >; // TODO JV (pt): Rethink allocator choice
#else
assert( false ); // TODO JV (pt): Make RingBuffer more general, so it could also be used here
#endif

#define DEBUG_ARCHIVER 1

namespace nest
{

class ConnectorModel;

/**
 * A node which archives spike history for the purposes of spike-timing
 * dependent plasticity (STDP)
 */
class ArchivingNode : public StructuralPlasticityNode
{
public:
  /**
   * \fn ArchivingNode()
   * Constructor.
   */
  ArchivingNode();

  /**
   * \fn ArchivingNode()
   * Copy Constructor.
   */
  ArchivingNode( const ArchivingNode& );

  /**
   * \fn double get_K_value(long t)
   * return the Kminus (synaptic trace) value at t (in ms). When the trace is
   * requested at the exact same time that the neuron emits a spike, the trace
   * value as it was just before the spike is returned.
   */
  // double get_K_value( double t ) override;

  std::pair< double, std::vector< double > > get_stdp_history( const double last_pre_spike_time,
    const double pre_spike_time,
    const double dendritic_delay,
    const synindex syn_id ) override;

  /**
   * \fn void get_K_values( double t,
   *   double& Kminus,
   *   double& nearest_neighbor_Kminus,
   *   double& Kminus_triplet )
   * write the Kminus (eligibility trace for STDP),
   * nearest_neighbour_Kminus (eligibility trace for nearest-neighbour STDP:
   *   like Kminus, but increased to 1, rather than by 1, on a spike
   *   occurrence),
   * and Kminus_triplet
   * values at t (in ms) to the provided locations.
   * @throws UnexpectedEvent
   */
  // void get_K_values( double t, double& Kminus, double& nearest_neighbor_Kminus, double& Kminus_triplet ) override;

  /**
   * \fn void get_K_values( double t,
   *   double& Kminus,
   *   double& Kminus_triplet )
   * The legacy version of the function, kept for compatibility
   * after changing the function signature in PR #865.
   * @throws UnexpectedEvent
   */
  //  void
  //  get_K_values( double t, double& Kminus, double& Kminus_triplet )
  //  {
  //    double nearest_neighbor_Kminus_to_discard;
  //    get_K_values( t, Kminus, nearest_neighbor_Kminus_to_discard, Kminus_triplet );
  //  }

  /**
   * \fn double get_K_triplet_value(std::deque<ArchivedSpikeTrace>::iterator &iter)
   * return the triplet Kminus value for the associated iterator.
   */
  // double get_K_triplet_value( const std::deque< ArchivedSpikeTrace >::iterator& iter );

  /**
   * \fn void get_history(long t1, long t2,
   * std::deque<Archiver::ArchivedSpikeTrace>::iterator* start,
   * std::deque<Archiver::ArchivedSpikeTrace>::iterator* finish)
   * return the spike times (in steps) of spikes which occurred in the range
   * (t1,t2].
   */
  // TODO JV (pt): This function could still exist, as there is no reason why the history should not be accessed. This
  //  function will be commented out for now to make sure no currently implemented model depends on it anymore.
  //  void get_history( double t1,
  //    double t2,
  //    std::deque< ArchivedSpikeTrace >::iterator* start,
  //    std::deque< ArchivedSpikeTrace >::iterator* finish ) override;

  /**
   * Register a new incoming STDP connection.
   *
   * t_first_read: The newly registered synapse will read the history entries
   * with t > t_first_read.
   */
  void register_stdp_connection( const double dendritic_delay, const synindex syn_id ) override;

  /**
   * Postponed delivery is required for STDP synapses with predominantly axonal delay. The archiving node supports this
   * feature by utilizing an intermediate spike buffer.
   */
  inline virtual bool
  supports_postponed_delivery() const
  {
    return false;
  }

  /**
   * When receiving an incoming event, forward it to the corresponding connection and handle the event updated by the
   * connection.
   */
  void deliver_event( const thread tid,
    const synindex syn_id,
    const index local_target_connection_id,
    const std::vector< ConnectorModel* >& cm,
    SpikeEvent& se ) override;

  void get_status( DictionaryDatum& d ) const override;
  void set_status( const DictionaryDatum& d ) override;

  /**
   * Framework for STDP with predominantly axonal delays:
   * Buffer a correction entry for a short time window.
   */
  void add_correction_entry_stdp_ax_delay( SpikeEvent& spike_event,
    const double t_last_pre_spike,
    const double weight_revert,
    const double dendritic_delay );


protected:
  void pre_run_hook_();

  /**
   * \fn void set_spiketime(Time const & t_sp, double offset)
   * record spike history
   */
  void set_spiketime( Time const& t_sp, double offset = 0.0 );

  /**
   * \fn double get_spiketime()
   * return most recent spike time in ms
   */
  double get_spiketime_ms() const;

  /**
   * Inform all incoming STDP connections of a post-synaptic spike to update the synaptic weight.
   */
  void update_stdp_connections( Time const& origin, const long from, const long to ) override;

  /**
   * clear spike history
   */
  void clear_history();

  void reset_correction_entries_stdp_ax_delay_();

  /**
   * Return if the node has any incoming stdp connections.
   */
  virtual bool has_stdp_connections() const;

private:
  // sum exp(-(t-ti)/tau_minus)
  // double Kminus_;

  // sum exp(-(t-ti)/tau_minus_triplet)
  // double Kminus_triplet_;

  double tau_minus_;
  double tau_minus_inv_;

  // time constant for triplet low pass filtering of "post" spike train
  double tau_minus_triplet_;
  double tau_minus_triplet_inv_;

  // double trace_;

  // double last_spike_;

  /**
   * Maximum dendritic delay of all incoming connections.
   */
  double max_dendritic_delay_;

  /**
   * Saves all synapse types registered to this node that need to get informed of post-synaptic spikes and support
   * axonal delay.
   * TODO JV (pt): There has to be a better (more object-oriented) way of handling specific connections differently.
   */
  std::set< synindex > stdp_synapse_types_;

  // spiking history needed by stdp synapses
  // TODO JV (pt): This has to be more generic somehow to support any synapse type, but without increasing memory usage
  //  for nodes without these special synapse type which need to store other values than just Kminus. In add_connection
  //  one will see which synapse type the neuron has to support and the synapses should actually manage the values in
  //  the history themselves somehow as this value depends on the exact synapse type. Even for "regular" STDP synapses
  //  the Kminus will be different for different pairing schemes.
  // SpikeArchive< ArchivedSpikeTrace > history_;
  std::deque< double > history_; // TODO JV (pt): Implement custom data structure that uses less memory than deque

  /**
   * Data structure to store incoming spikes sent over STDP synapses with predominantly axonal delay. If at time of
   * delivery there might occur any post-synaptic spike that reaches the synapse before the pre-synaptic spike
   * (i.e., after the dendritic delay), the pre-synaptic spike has to be buffered until no more critical post-synaptic
   * spike can occur.
   * Arranged in a 2d structure: X|X
   */
  // std::vector <  > intermediate_spike_buffer_;

  /**
   * Framework for STDP with predominantly axonal delays:
   * Due to the long axonal delays, relevant spikes of this neuron might not yet be available at the time when incoming
   * synapses are updated (spike delivery). Therefore, for each spike received through an STDP synapse with
   * predominantly axonal delay, information is stored for a short period of time allowing for retrospective correction
   * of the synapse and the already delivered spike.
   */
  struct CorrectionEntrySTDPAxDelay
  {
    CorrectionEntrySTDPAxDelay( const synindex syn_id,
      const index local_connection_id,
      const double t_last_pre_spike,
      const double weight_revert )
      : syn_id_( syn_id )
      , local_connection_id_( local_connection_id )
      , t_last_pre_spike_( t_last_pre_spike )
      , weight_revert_( weight_revert )
    {
    }

    synindex syn_id_;           //!< index of synapse type
    index local_connection_id_; //!< index of connection in node
    double t_last_pre_spike_;   //!< time of the last pre-synaptic spike before this spike
    double weight_revert_;      //!< synaptic weight to revert to (STDP depression needs to be undone)
  };

  /**
   * Framework for STDP with predominantly axonal delays:
   * Buffer of correction entries sorted by t_spike_pre + delay
   * (i.e., the actual arrival time at this neuron).
   */
  std::vector< std::vector< CorrectionEntrySTDPAxDelay > > correction_entries_stdp_ax_delay_;

  /**
   * Framework for STDP with predominantly axonal delays:
   * Triggered when this neuron spikes, to correct all relevant incoming STDP synapses with predominantly axonal delays
   * and corresponding received spikes.
   */
  void correct_synapses_stdp_ax_delay_( const Time& t_spike );
};

inline double
ArchivingNode::get_spiketime_ms() const
{
  assert( false ); // TODO JV (pt)
  // return last_spike_;
}

inline std::pair< double, std::vector< double > >
ArchivingNode::get_stdp_history( const double last_pre_spike_time,
  const double pre_spike_time,
  const double dendritic_delay,
  const synindex syn_id )
{
  return connections_[ syn_id ]->get_stdp_history(
    last_pre_spike_time, pre_spike_time, dendritic_delay, tau_minus_inv_, history_.cbegin(), history_.cend() );
}

} // of namespace
#endif
