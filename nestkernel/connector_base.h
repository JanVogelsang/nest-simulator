/*
 *  connector_base.h
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

#ifndef CONNECTOR_BASE_H
#define CONNECTOR_BASE_H

// C++ includes:
#include <vector>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection_label.h"
#include "nest_datums.h"
#include "nest_names.h"
#include "spikecounter.h"

// Includes from sli:
#include "archived_spike.h"
#include "arraydatum.h"
#include "dictutils.h"

namespace nest
{

class ConnectorModel;
class Event;

/**
 * Base class to allow storing Connectors for different synapse types
 * in vectors. We define the interface here to avoid casting.
 */
class ConnectorBase
{

public:
  // Destructor needs to be declared virtual to avoid undefined
  // behavior, avoid possible memory leak and needs to be defined to
  // avoid linker error, see, e.g., Meyers, S. (2005) p40ff
  virtual ~ConnectorBase() {};

  /**
   * Return syn_id_ of the synapse type of this Connector (index in list of synapse prototypes).
   */
  virtual synindex get_syn_id() const = 0;

  /**
   * Return the number of connections/sources in this Connector.
   */
  virtual size_t size() const = 0;

  /**
   * Write status of the connection at position lcid to the dictionary
   * dict.
   */
  virtual void get_synapse_status( const thread tid, const index lcid, DictionaryDatum& dict ) const = 0;

  /**
   * Set status of the connection at position lcid according to the
   * dictionary dict.
   */
  virtual void set_synapse_status( const index lcid, const DictionaryDatum& dict, ConnectorModel& cm ) = 0;

  /**
   * Get the proportion of the transmission delay attributed to the dendrite of a connection.
   */
  virtual delay get_dendritic_delay( const index lcid ) const = 0;

  /**
   * Get the time of the last pre-synaptic spike that was emitted over the specified connection.
   */
  virtual double get_last_presynaptic_spike( const index lcid ) const = 0;

#ifndef USE_ADJACENCY_LIST
  /**
   * Get the proportion of the transmission delay attributed to the axon of a connection.
   */
  virtual double get_axonal_delay( const index lcid ) const = 0;
#endif

  /**
   * Get the source node id of a specific connection.
   */
  virtual index get_source( const index lcid ) = 0;

  /**
   * Get all source node ids.
   */
  virtual std::vector< index > get_sources() = 0;

  /**
   * Get the indices of all connections corresponding to a specific source node id with specific label.
   */
  virtual std::vector< index > get_connection_indices( const index source_node_id,
    const long connection_label = UNLABELED_CONNECTION ) const = 0;

  /**
   * Get the indices of all connections with specific label.
   */
  virtual std::vector< index > get_connection_indices( const long connection_label = UNLABELED_CONNECTION ) const = 0;

  /**
   * Remove source information of all connections in this container.
   */
  virtual void clear_sources() = 0;

  virtual long get_connection_label( const index lcid ) const = 0;

  virtual bool is_connection_disabled( const index lcid ) const = 0;

  /**
   * Send the event e to the connection at position lcid. Return bool
   * indicating whether the following connection belongs to the same
   * source.
   */
  virtual void send( const thread tid,
    const index lcid,
    const delay axonal_delay,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) = 0;

  virtual void send( const thread tid,
    const index lcid,
    const delay axonal_delay,
    const delay dendritic_delay,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) = 0;

  virtual void send_weight_event( const thread tid,
    const index lcid,
    Event& e,
    const CommonSynapseProperties& cp,
    Node* target ) = 0;

  virtual void update_stdp_connections( const double post_spike_time_syn,
    const delay dendritic_delay,
    const ConnectorModel* cm ) = 0;

  virtual void
  process_post_synaptic_spike( const index lcid, const double post_spike_time_syn, const ConnectorModel* cm ) = 0;

  virtual void
  update_trace( const double post_spike_time, const delay dendritic_delay, const double tau_minus_inv ) = 0;

  virtual double get_trace( const double pre_spike_time,
    const double dendritic_delay,
    const double tau_minus_inv,
    const std::deque< double >::const_iterator history_begin,
    const std::deque< double >::const_iterator history_end ) = 0;

  virtual std::pair< double, double > get_Kminus( const double dendritic_delay ) = 0;

  virtual void clear_history() = 0;

  /**
   * Update weights of dopamine modulated STDP connections.
   */
  virtual void trigger_update_weight( const long vt_node_id,
    const thread tid,
    const std::vector< spikecounter >& dopa_spikes,
    const double t_trig,
    const std::vector< ConnectorModel* >& cm ) = 0;

  /**
   * Sort connections according to dendritic delay.
   */
  virtual void prepare_connections( const thread tid, const index target_lid ) = 0;

  /**
   * Disable the transfer of events through the connection at position
   * lcid.
   */
  virtual void disable_connection( const index lcid ) = 0;

  /**
   * Remove disabled connections from the connector.
   */
  virtual void remove_disabled_connections( const index first_disabled_index ) = 0;
};

/**
 * Homogeneous connector, contains synapses of one particular type (syn_id_).
 */
template < typename ConnectionT >
class Connector : public ConnectorBase
{
  //! A region of connections in the Connector where all connections share the same dendritic delay
  struct DelayRegion
  {
    index start; //!< first connection of delay region
    index end;   //!< last connection of relay region
    //! the post-synaptic trace in synapse time, i.e. the trace of the neuron "dendritic delay"-milliseconds ago
    double Kminus;
    double last_post_spike; //!< time of the last update in synapse time
  };

private:
  const synindex syn_id_;
  std::map< delay, DelayRegion > dendritic_delay_regions_;
  std::map< delay, std::vector< index > > connection_indices_by_delay_;
  std::vector< delay > axonal_delays_;
  std::vector< ConnectionT > C_;

  /**
   * This data structure stores the node IDs of presynaptic neurons connected to this neuron. If structural plasticity
   * is disabled, it is only relevant during postsynaptic connection creation, before the connection information has
   * been transferred to the presynaptic side.
   * Arranged in a 2d array:
   * 1st dimension: synapse types
   * 2nd dimension: source node IDs
   * After all connections have been created, the information stored in this structure is transferred to the presynaptic
   * side and the sources vector can be cleared, unless further required for structural plasticity.
   */
  std::vector< index > sources_;

public:
  explicit Connector( const synindex syn_id )
    : syn_id_( syn_id )
  {
    // TODO JV: Benchmark this both with and without reserve
    axonal_delays_.reserve( 11250 );
    C_.reserve( 11250 );
    sources_.reserve( 11250 );
  }

  ~Connector() override
  {
    std::vector< ConnectionT >().swap( C_ );
    auto s = connection_indices_by_delay_.begin()->second;
    std::vector< index >().swap( sources_ );
    std::map< delay, DelayRegion >().swap( dendritic_delay_regions_ );
    std::map< delay, std::vector< index > >().swap( connection_indices_by_delay_ );
#ifndef USE_ADJACNECY_LIST
    std::vector< delay >().swap( axonal_delays_ );
#endif
  }

  synindex
  get_syn_id() const override
  {
    return syn_id_;
  }

  size_t
  size() const override
  {
    return C_.size();
  }

  void
  get_synapse_status( const thread tid, const index lcid, DictionaryDatum& dict ) const override
  {
    assert( lcid < C_.size() );

    C_[ lcid ].get_status( dict );

    // get target node ID here, where tid is available
    // necessary for hpc synapses using TargetIdentifierIndex
    // def< long >( dict, names::target, C_[ lcid ].get_target( tid )->get_node_id() );
  }

  void
  set_synapse_status( const index lcid, const DictionaryDatum& dict, ConnectorModel& cm ) override
  {
    assert( lcid < C_.size() );

    C_[ lcid ].set_status( dict, static_cast< GenericConnectorModel< ConnectionT >& >( cm ) );
  }

  delay
  get_dendritic_delay( const index lcid ) const override
  {
    assert( lcid < C_.size() );

    // TODO JV: Profiling - lower bound vs find vs simply iterating
    // TODO JV (pt): Further optimize this as much as possible (after additional profiling)
    auto v = std::lower_bound( dendritic_delay_regions_.begin(),
      dendritic_delay_regions_.end(),
      lcid,
      []( const std::pair< const delay, DelayRegion >& r, const index idx ) -> const bool
      { return r.second.end < idx; } );
    return v->first;
  }

  double
  get_last_presynaptic_spike( const index lcid ) const override
  {
    return C_[ lcid ].get_last_presynaptic_spike();
  }

#ifndef USE_ADJACENCY_LIST
  double
  get_axonal_delay( const index lcid ) const override
  {
    assert( lcid < axonal_delays_.size() );

    return axonal_delays_[ lcid ];
  }
#endif

  index
  get_source( const index lcid ) override
  {
    return sources_[ lcid ];
  }

  std::vector< index >
  get_sources() override
  {
    return sources_;
  }

  const index
  add_device_connection( ConnectionT& c, const index source_node_id )
  {
    C_.push_back( c );
    sources_.push_back( source_node_id );
    if ( C_.size() > MAX_LOCAL_CONNECTION_ID )
    {
      throw KernelException(
        String::compose( "Too many connections: at most %1 connections supported per virtual "
                         "process and synapse model to a specific target neuron.",
          MAX_LOCAL_CONNECTION_ID ) );
    }
    return C_.size() - 1;
  }

  void
  add_connection( const ConnectionT& c,
    const index source_node_id,
    const delay axonal_delay,
    const delay dendritic_delay )
  {
    connection_indices_by_delay_[ dendritic_delay ].push_back(
      C_.size() ); // TODO JV: Save continuous indices as start and end instead of each individual
    axonal_delays_.push_back( axonal_delay );
    C_.push_back( c );
    sources_.push_back( source_node_id );
    if ( C_.size() > MAX_LOCAL_CONNECTION_ID )
    {
      throw KernelException(
        String::compose( "Too many connections: at most %1 connections supported per virtual "
                         "process and synapse model to a specific target neuron.",
          MAX_LOCAL_CONNECTION_ID ) );
    }
  }

  std::vector< index >
  get_connection_indices( const index source_node_id,
    const long connection_label = UNLABELED_CONNECTION ) const override
  {
    // binary search in sorted sources
    std::vector< index >::const_iterator it = sources_.cbegin();
    const std::vector< index >::const_iterator end = sources_.cend();

    std::vector< index > indices;
    while ( ( it = std::find( it, end, source_node_id ) ) != sources_.end() )
    {
      const index lcid = std::distance( sources_.cbegin(), it );
      // Connection is disabled
      if ( not C_[ lcid ].is_disabled() )
      {
        if ( connection_label == UNLABELED_CONNECTION or C_[ lcid ].get_label() == connection_label )
        {
          indices.push_back( lcid );
        }
      }
      ++it;
    }

    return indices;
  }

  std::vector< index >
  get_connection_indices( const long connection_label = UNLABELED_CONNECTION ) const override
  {
    std::vector< index > indices;
    for ( auto conn_it = C_.cbegin(); conn_it != C_.cend(); ++conn_it )
    {
      const index lcid = std::distance( C_.cbegin(), conn_it );
      // Connection is disabled
      if ( not C_[ lcid ].is_disabled() )
      {
        if ( connection_label == UNLABELED_CONNECTION or C_[ lcid ].get_label() == connection_label )
        {
          indices.push_back( lcid );
        }
      }
    }
    return indices;
  }

  long
  get_connection_label( const index lcid ) const override
  {
    return C_[ lcid ].get_label();
  }

  bool
  is_connection_disabled( const index lcid ) const override
  {
    return C_[ lcid ].is_disabled();
  }

  void
  clear_sources() override
  {
    std::vector< index >().swap( sources_ );
  }

  void
  send( const thread tid,
    const index lcid,
    const delay axonal_delay,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) override
  {
    send( tid, lcid, axonal_delay, get_dendritic_delay( lcid ), cm, e, target );
  }

  void
  send( const thread tid,
    const index lcid,
    const delay axonal_delay,
    const delay dendritic_delay,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) override
  {
    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties();

    e.set_port( lcid );
    if ( not C_[ lcid ].is_disabled() )
    {
      C_[ lcid ].send( e, tid, axonal_delay, dendritic_delay, cp, target );
      send_weight_event( tid, lcid, e, cp, target );
    }
  }

  std::pair< double, double >
  get_Kminus( const double dendritic_delay ) override
  {
    // get the post-synaptic trace in synapse time, i.e. the trace of the neuron "dendritic delay"-milliseconds ago
    const delay dendritic_delay_steps = Time::delay_ms_to_steps( dendritic_delay );
    const DelayRegion& delay_region = dendritic_delay_regions_.find( dendritic_delay_steps )->second;
    return { delay_region.last_post_spike, delay_region.Kminus };
  }

  void send_weight_event( const thread tid,
    const index lcid,
    Event& e,
    const CommonSynapseProperties& cp,
    Node* target ) override;

  void
  trigger_update_weight( const long vt_node_id,
    const thread tid,
    const std::vector< spikecounter >& dopa_spikes,
    const double t_trig,
    const std::vector< ConnectorModel* >& cm ) override
  {
    for ( size_t i = 0; i < C_.size(); ++i )
    {
      if ( static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )
             ->get_common_properties()
             .get_vt_node_id()
        == vt_node_id )
      {
        C_[ i ].trigger_update_weight( tid,
          dopa_spikes,
          t_trig,
          get_dendritic_delay( i ),
          static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties() );
      }
    }
  }

  void update_stdp_connections( const double post_spike_time_syn,
    const delay dendritic_delay,
    const ConnectorModel* cm ) override;

  void
  process_post_synaptic_spike( const index lcid, const double post_spike_time_syn, const ConnectorModel* cm ) override
  {
    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< const GenericConnectorModel< ConnectionT >* >( cm )->get_common_properties();
    C_[ lcid ].process_post_synaptic_spike( post_spike_time_syn, cp );
  }

  void
  update_trace( const double post_spike_time, const delay dendritic_delay, const double tau_minus_inv ) override
  {
    auto group_it = dendritic_delay_regions_.find( dendritic_delay );

    if ( group_it != dendritic_delay_regions_.end() )
    {
      // update post-synaptic trace
      group_it->second.Kminus = group_it->second.Kminus
          * std::exp(
            ( group_it->second.last_post_spike - ( post_spike_time + Time::delay_steps_to_ms( dendritic_delay ) ) )
            * tau_minus_inv )
        + 1;
      group_it->second.last_post_spike = post_spike_time + Time::delay_steps_to_ms( dendritic_delay );
    }
  }

  double get_trace( const double pre_spike_time,
    const double dendritic_delay,
    const double tau_minus_inv,
    const std::deque< double >::const_iterator history_begin,
    const std::deque< double >::const_iterator history_end ) override;

  void
  clear_history() override
  {
    for ( auto& group : dendritic_delay_regions_ )
    {
      group.second.Kminus = 0;
      group.second.last_post_spike = -1;
    }
  }

  void prepare_connections( const thread tid, const index target_lid ) override;

  void
  disable_connection( const index lcid ) override
  {
    assert( false ); // TODO JV (pt): Structural plasticity

    /*assert( not C_[ lcid ].is_disabled() );
    C_[ lcid ].disable();
    sources_[ lcid ] = DISABLED_NODE_ID;*/
  }

  void
  remove_disabled_connections( const index first_disabled_index ) override
  {
    assert( false ); // TODO JV (pt): Structural plasticity

    /*assert( C_[ first_disabled_index ].is_disabled() );
    C_.erase( C_.begin() + first_disabled_index, C_.end() );*/
  }
};

} // of namespace nest

#endif
