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
  ConnectorBase() : last_visited_connection_( 0 ) {}

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
   * Get the proportion of the transmission delay attributed to the dendrite.
   */
  virtual double get_connection_delay( const index lcid, ConnectorModel& cm ) = 0;

  /**
   * Get the source node id of a specific connection.
   */
  virtual index& get_source( const index local_target_connection_id ) = 0;

  /**
   * Get the first and last index of all connections corresponding to a specific source node.
   */
  virtual std::pair< index, index > get_connection_indices( const index source_node_id ) const = 0;

  /**
   * Get the first index of all connections corresponding to a specific source node.
   */
  virtual index get_first_connection_index( const index source_node_id ) const = 0;

  /**
   * Sets the index of the connection that last received a spike.
   */
  void set_last_visited_connection( const size_t last_visited_connection )
  {
    last_visited_connection_ = last_visited_connection;
  }

  /**
   * Returns the index of the connection that last received a spike.
   */
  size_t get_last_visited_connection()
  {
    return last_visited_connection_;
  }

  void reset_last_visited_connection()
  {
    last_visited_connection_ = 0;
  }

  /**
   * Remove source information of all connections in this container.
   */
  virtual void clear_sources() = 0;

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
    const ConnectorModel* cm,
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

  /**
   * Update weights of dopamine modulated STDP connections.
   */
  virtual void trigger_update_weight( const long vt_node_id,
    const thread tid,
    const std::vector< spikecounter >& dopa_spikes,
    const double t_trig,
    const std::vector< ConnectorModel* >& cm ) = 0;

  /**
   * Sort connections according to source node IDs.
   */
  virtual void sort_connections_and_sources() = 0;

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
  // TODO JV (help): Would it work to precompute the start-indices per source rank after sorting to save time when searching?
  //  What else could we do here to speed up the lookup? Most searching algorithms are quite fast when it comes to
  //  finding an element in a sorted sequence, but is it also equally fast to guarantee an element is not in the list as
  //  this will be the case which will occur much more often here. Some sort of pre-filtering might make sense here as
  //  well, by somehow intelligently analyzing all sources and thus enabling faster lookup (but how?).
  std::vector< index > sources_;

  size_t last_visited_connection_;
};

/**
 * Homogeneous connector, contains synapses of one particular type (syn_id_).
 */
template < typename ConnectionT >
class Connector : public ConnectorBase
{
private:
  std::vector< ConnectionT > C_;
  const synindex syn_id_;

public:
  explicit Connector( const synindex syn_id )
    : ConnectorBase()
    , syn_id_( syn_id )
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
  get_connection_delay( const index lcid, ConnectorModel& cm ) override
  {
    assert( lcid < C_.size() );

    return C_[ lcid ].get_delay();
  }

  const index
  add_connection( const ConnectionT& c, const index source_node_id )
  {
    C_.push_back( c );
    sources_.push_back( source_node_id );
    // Return index of added item
    return C_.size() - 1;
  }

  index&
  get_source( const index local_target_connection_id ) override
  {
    return sources_[ local_target_connection_id ];
  }

  std::pair< index, index >
  get_connection_indices( const index source_node_id ) const override
  {
    assert( source_node_id >= sources_[ last_visited_connection_ ] );

    // binary search in sorted sources
    const std::vector< index >::const_iterator begin = sources_.begin() + last_visited_connection_;
    const std::vector< index >::const_iterator end = sources_.end();
    auto first_last_source = std::equal_range( begin, end, source_node_id );

    return { first_last_source.first - begin, first_last_source.second - begin };
  }

  index
  get_first_connection_index( const index source_node_id ) const override
  {
    // binary search in sorted sources
    const std::vector< index >::const_iterator begin = sources_.begin() + last_visited_connection_;
    const std::vector< index >::const_iterator end = sources_.end();

    auto first_source = std::lower_bound( begin, end, source_node_id );

    return first_source - begin;
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
    std::vector< index >().swap( sources_ );
  }

  void
  send( const thread tid,
    const ConnectorModel* cm,
    Event& e,
    Node* target ) override
  {
    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< const GenericConnectorModel< ConnectionT >* >( cm )->get_common_properties();

    const index local_target_connection_id = e.get_local_connection_id();
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

  void correct_synapse_stdp_ax_delay( const index local_target_connection_id,
    const double t_last_pre_spike,
    double* weight_revert,
    const double t_post_spike,
    Node* target ) override;

  void
  sort_connections_and_sources() override
  {
    nest::sort( sources_, C_ );
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
