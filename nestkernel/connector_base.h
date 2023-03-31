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
  virtual double get_dendritic_delay( const index lcid ) const = 0;

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
  virtual std::vector< index > get_connection_indices( const index source_node_id, const long connection_label = UNLABELED_CONNECTION ) const = 0;

    /**
     * Get the indices of all connections with specific label.
     */
  virtual std::vector< index > get_connection_indices( const long connection_label = UNLABELED_CONNECTION ) const = 0;

  /**
   * Remove source information of all connections in this container.
   */
  virtual void clear_sources() = 0;

  virtual long
  get_connection_label( const index lcid ) const = 0;

  virtual bool
  is_connection_disabled( const index lcid ) const = 0;

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

  virtual void correct_synapse_stdp_ax_delay( const index local_target_connection_id,
    const double t_last_pre_spike,
    const double axonal_delay,
    double* weight_revert,
    const double t_post_spike,
    Node* target ) = 0;

  virtual void send_weight_event( const thread tid,
    const index lcid,
    Event& e,
    const CommonSynapseProperties& cp,
    Node* target ) = 0;

  /**
   * Update weights of dopamine modulated STDP connections.
   */
  virtual void trigger_update_weight( const long vt_node_id,
    const thread tid,
    const std::vector< spikecounter >& dopa_spikes,
    const double t_trig,
    const std::vector< ConnectorModel* >& cm ) = 0;

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
private:
  const synindex syn_id_;
  std::vector< ConnectionT > C_;
  std::vector< delay > axonal_delays_;

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

  ~Connector() override = default;
  //  ~Connector() override
  //  {
  //    C_.clear();
  //  }

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
  get_dendritic_delay( const index lcid ) const override
  {
    assert( lcid < C_.size() );

    return C_[ lcid ].get_dendritic_delay();
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
  add_connection( const ConnectionT& c, const index source_node_id, const delay axonal_delay )
  {
    C_.push_back( c );
    sources_.push_back( source_node_id );
#ifndef USE_ADJACENCY_LIST
    axonal_delays_.push_back( axonal_delay );
#endif
    if ( C_.size() > MAX_LOCAL_CONNECTION_ID )
    {
      throw KernelException(
        String::compose( "Too many connections: at most %1 connections supported per virtual "
                         "process and synapse model to a specific target neuron.",
          MAX_LOCAL_CONNECTION_ID ) );
    }
    // Return index of added item
    return C_.size() - 1;
  }

  std::vector< index >
  get_connection_indices( const index source_node_id, const long connection_label = UNLABELED_CONNECTION ) const override
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
    sources_.clear();
  }

  void
  send( const thread tid,
    const index lcid,
    const delay axonal_delay,
    const std::vector< ConnectorModel* >& cm,
    Event& e,
    Node* target ) override
  {
    assert( C_.size() > lcid );

    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties();

    e.set_port( lcid );
    if ( not C_[ lcid ].is_disabled() )
    {
      C_[ lcid ].send( e, tid, axonal_delay, cp, target );
      send_weight_event( tid, lcid, e, cp, target );
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

  void correct_synapse_stdp_ax_delay( const index local_target_connection_id,
    const double t_last_pre_spike,
    const double axonal_delay,
    double* weight_revert,
    const double t_post_spike,
    Node* target ) override;

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
