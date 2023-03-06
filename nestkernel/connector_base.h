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
#include "nest_datums.h"
#include "nest_names.h"
#include "source.h"
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
  virtual double get_dendritic_delay( const index lcid ) = 0;

  /**
   * Get the proportion of the transmission delay attributed to the axon of a connection.
   */
  virtual double get_axonal_delay( const index lcid ) = 0;

  /**
   * Get information about the source node of a specific connection.
   */
  virtual Source& get_source( const index local_target_connection_id ) = 0;

  /**
   * Get the indices of all connection corresponding to a specific source node.
   */
  virtual std::vector< index > get_connection_indices( const index source_node_id ) const = 0;

  /**
   * Remove source information of all connections in this container.
   */
  virtual void clear_sources() = 0;

  /**
   * Reset all processed flags of all sources in this container.
   */
  virtual void reset_sources_processed_flags() = 0;

  /**
   * Add ConnectionID with given source_node_id and lcid to conns. If
   * target_neuron_node_ids is given, only add connection if
   * target_neuron_node_ids contains the node ID of the target of the connection.
   */
  virtual void get_connection_with_specified_targets( const index source_node_id,
    const std::vector< size_t >& target_neuron_node_ids,
    const thread tid,
    const index lcid,
    const long synapse_label,
    std::deque< ConnectionID >& conns ) const = 0;

  /**
   * Send the event e to the connection at position lcid. Return bool
   * indicating whether the following connection belongs to the same
   * source.
   */
  virtual void send( const thread tid,
    const index local_target_connection_id,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) = 0;

  virtual void send( const thread tid,
    const index local_target_connection_id,
    const delay dendritic_delay,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) = 0;

  virtual void send_weight_event( const thread tid,
    const index local_target_connection_id,
    Event& e,
    const CommonSynapseProperties& cp,
    Node* target ) = 0;

  virtual void update_stdp_connections( const double post_spike_time,
    const double dendritic_delay,
    const double tau_minus_inv,
    const std::vector< ConnectorModel* >& cm ) = 0;

  virtual std::pair< double, std::vector< double > > get_stdp_history( const double last_pre_spike_time,
    const double pre_spike_time,
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
  virtual void prepare_connections() = 0;

  /**
   * Disable the transfer of events through the connection at position
   * lcid.
   */
  virtual void disable_connection( const index lcid ) = 0;

  /**
   * Remove disabled connections from the connector.
   */
  virtual void remove_disabled_connections( const index first_disabled_index ) = 0;


protected:
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
  std::vector< Source > sources_; // TODO JV: Move to Connector
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
    index start; //<! first connection of delay region
    index end;   //<! last connection of relay region
    //! the post-synaptic trace in synapse time, i.e. the trace of the neuron "dendritic delay"-milliseconds ago
    double Kminus;
    double last_post_spike; //<! time of the last update in synapse time
  };

private:
  const synindex syn_id_;
  std::map< delay, DelayRegion > dendritic_delay_regions_;
  std::map< delay, std::vector< index > > connection_indices_by_delay_;
  std::vector< ConnectionT > C_;

public:
  explicit Connector( const synindex syn_id )
    : syn_id_( syn_id )
  {
  }

  ~Connector() override
  {
    C_.clear();
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

  double
  get_dendritic_delay( const index lcid ) override
  {
    assert( lcid < C_.size() );

    return std::lower_bound( dendritic_delay_regions_.begin(),
      dendritic_delay_regions_.end(),
      lcid,
      []( const std::pair< const delay, DelayRegion >& r, const index idx ) -> bool { return r.second.start < idx; } )
      ->first;
  }

  double
  get_axonal_delay( const index lcid ) override
  {
    assert( lcid < C_.size() );

    return C_[ lcid ].get_axonal_delay();
  }

  const index
  add_connection( const ConnectionT& c, const Source src, const delay dendritic_delay )
  {
    const size_t next_index = C_.size();
    connection_indices_by_delay_[ dendritic_delay ].push_back( next_index );
    dendritic_delay_regions_[ dendritic_delay ];
    C_.push_back( c );
    sources_.push_back( src );
    // Return index of added item
    return next_index;
  }

  Source&
  get_source( const index local_target_connection_id ) override
  {
    return sources_[ local_target_connection_id ];
  }

  std::vector< index >
  get_connection_indices( const index source_node_id ) const override
  {
    std::vector< index > indices;

    std::vector< Source >::const_iterator it = sources_.cbegin();
    while ( ( it = std::find_if( it,
                sources_.cend(),
                [ source_node_id ]( const Source src ) { return src.get_node_id() == source_node_id; } ) )
      != sources_.end() )
    {
      // Do something with iter
      indices.push_back( std::distance( sources_.begin(), it ) );
      ++it;
    }

    return indices;
  }

  void
  get_connection_with_specified_targets( const index source_node_id,
    const std::vector< size_t >& target_neuron_node_ids,
    const thread tid,
    const index lcid,
    const long synapse_label,
    std::deque< ConnectionID >& conns ) const override
  {
    assert( false ); // TODO JV (pt): Structural plasticity

    /*if ( not C_[ lcid ].is_disabled() )
    {
      if ( synapse_label == UNLABELED_CONNECTION or C_[ lcid ].get_label() == synapse_label )
      {
        const index current_target_node_id = C_[ lcid ].get_target( tid )->get_node_id();
        if ( std::find( target_neuron_node_ids.begin(), target_neuron_node_ids.end(), current_target_node_id )
          != target_neuron_node_ids.end() )
        {
          conns.push_back(
            ConnectionDatum( ConnectionID( source_node_id, current_target_node_id, tid, syn_id_, lcid ) ) );
        }
      }
    }*/
  }

  void
  clear_sources() override
  {
    // TODO JV: Is this actually safe? Or is swap the better way of clearing 2D vectors?
    sources_.clear();
  }

  void
  reset_sources_processed_flags() override
  {
    for ( std::vector< Source >::iterator source = sources_.begin(); source != sources_.end(); ++source )
    {
      source->set_processed( false );
    }
  }

  void
  send( const thread tid,
    const index local_target_connection_id,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) override
  {
    send( tid, local_target_connection_id, get_dendritic_delay( local_target_connection_id ), cm, e, target );
  }

  void
  send( const thread tid,
    const index local_target_connection_id,
    const delay dendritic_delay,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) override
  {
    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties();

    e.set_port( local_target_connection_id );
    if ( not C_[ local_target_connection_id ].is_disabled() )
    {
      C_[ local_target_connection_id ].send( e, tid, dendritic_delay, cp, target );
      send_weight_event( tid, local_target_connection_id, e, cp, target );
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
    const index local_target_connection_id,
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

  void
  update_stdp_connections( const double post_spike_time,
    const double dendritic_delay,
    const double tau_minus_inv,
    const std::vector< ConnectorModel* >& cm ) override
  {
    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties();

    // TODO JV (pt): Verify for precise spike times
    // TODO JV: Make sure double->long conversion yields the correct delay
    auto group_it = dendritic_delay_regions_.find( Time::delay_ms_to_steps( dendritic_delay ) );
    if ( group_it != dendritic_delay_regions_.end() ) // check if there are connections with given dendritic delay
    {
      const auto end = C_.begin() + group_it->second.end;
      for ( auto it = C_.begin() + group_it->second.start; it != end; ++it )
      {
        it->process_post_synaptic_spike( post_spike_time + dendritic_delay, cp );
      }

      // update post-synaptic trace
      group_it->second.Kminus = group_it->second.Kminus
          * std::exp( ( group_it->second.last_post_spike - ( post_spike_time + dendritic_delay ) ) * tau_minus_inv )
        + 1;
      group_it->second.last_post_spike = post_spike_time + dendritic_delay;
    }
  }

  std::pair< double, std::vector< double > > get_stdp_history( const double last_pre_spike_time,
    const double pre_spike_time,
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

  void
  prepare_connections() override
  {
    std::vector< ConnectionT > temp_connections;
    std::vector< Source > temp_sources;
    temp_connections.reserve( C_.size() );
    temp_sources.reserve( C_.size() );
    for ( const auto& region : connection_indices_by_delay_ )
    {
      const auto region_start = temp_sources.cend();
      for ( auto idx : region.second )
      {
        temp_connections.push_back( std::move( C_[ idx ] ) );
        temp_sources.push_back( std::move( sources_[ idx ] ) );
      }
      dendritic_delay_regions_[ region.first ] = { static_cast< index >(
                                                     std::distance( temp_sources.cbegin(), region_start ) ),
        static_cast< index >( std::distance( temp_sources.cbegin(), temp_sources.cend() ) ),
        0,
        -1 };
    }
    C_ = std::move( temp_connections );
    sources_ = std::move( temp_sources );
  }

  void
  disable_connection( const index lcid ) override
  {
    assert( false ); // TODO JV (pt): Structural plasticity

    /*assert( not C_[ lcid ].is_disabled() );
    C_[ lcid ].disable();
    sources_[ lcid ].disable();*/
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
