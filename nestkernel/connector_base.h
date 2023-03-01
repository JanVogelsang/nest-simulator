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

// Includes from libnestutil:
#include "sort.h"

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

  virtual void correct_synapse_stdp_ax_delay( const index local_target_connection_id,
    const double t_last_pre_spike,
    double* weight_revert,
    const double t_post_spike,
    Node* target ) = 0;

  virtual void send_weight_event( const thread tid,
    const index local_target_connection_id,
    Event& e,
    const CommonSynapseProperties& cp,
    Node* target ) = 0;

  virtual void update_stdp_connections( const ArchivedSpikeTrace &archived_spike, const double dendritic_delay ) = 0;

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
  virtual void sort_connections_by_dendritic_delay() = 0;

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
  // TODO JV: This should be converted from type Source to index once the simulation starts
  std::vector< Source > sources_;
};

/**
 * Homogeneous connector, contains synapses of one particular type (syn_id_).
 */
template < typename ConnectionT >
class Connector : public ConnectorBase
{
private:
  const synindex syn_id_;
  std::map< delay, std::pair< typename std::vector< ConnectionT >::iterator, typename std::vector< ConnectionT >::iterator > > dendritic_delay_regions_;
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

    return C_[ lcid ].get_dendritic_delay();
  }

  double
  get_axonal_delay( const index lcid ) override
  {
    assert( lcid < C_.size() );

    return C_[ lcid ].get_axonal_delay();
  }

  const index
  add_connection( const ConnectionT& c, const Source src )
  {
    C_.push_back( c );
    sources_.push_back( src );
    // Return index of added item
    return C_.size() - 1;
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
    while ( ( it = std::find_if( it, sources_.cend(), [ source_node_id ]( const Source src ) { return src.get_node_id() == source_node_id; } ) ) != sources_.end() )
    {
      // Do something with iter
      indices.push_back( std::distance( sources_.begin(), it ));
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
    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties();

    e.set_port( local_target_connection_id );
    if ( not C_[ local_target_connection_id ].is_disabled() )
    {
      C_[ local_target_connection_id ].send( e, tid, cp, target );
      send_weight_event( tid, local_target_connection_id, e, cp, target );
    }
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
          static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties() );
      }
    }
  }

  void update_stdp_connections( const ArchivedSpikeTrace &archived_spike, const double dendritic_delay ) override
  {
    // TODO JV (pt): Verify for precise spike times
    // TODO JV: Make sure double->long conversion yields the correct delay
    auto group_it = dendritic_delay_regions_.find( Time::delay_ms_to_steps(dendritic_delay) );
    if ( group_it != dendritic_delay_regions_.end() )  // check if there are connections with given dendritic delay
    {
      const auto end = group_it->second.second;
      for ( auto it = group_it->second.first; it != end; ++it )
      {
        // TODO JV (help): How to only implement this function on some synapse types and still make this compile?
        // it->process_post_synaptic_spike( archived_spike.t + dendritic_delay, cp );
      }
    }
  }

  void correct_synapse_stdp_ax_delay( const index local_target_connection_id,
    const double t_last_pre_spike,
    double* weight_revert,
    const double t_post_spike,
    Node* target ) override;

  void
  sort_connections_by_dendritic_delay() override
  {
    std::sort( C_.begin(), C_.end(), []( const ConnectionT& lhs, const ConnectionT& rhs )
      { return lhs.get_dendritic_delay() < rhs.get_dendritic_delay(); } );

    auto start_it = C_.begin();
    delay last_delay = start_it->get_dendritic_delay();
    for ( auto it = C_.begin() + 1; it != C_.end(); ++it )
    {
      if ( it->get_dendritic_delay() != last_delay )  // detect switch to new dendritic delay region
      {
        dendritic_delay_regions_[ last_delay ] = { start_it, it };
        start_it = it;
        last_delay = it->get_dendritic_delay();
      }
    }
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
