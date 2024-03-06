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
  virtual void
  get_synapse_status( const index lcid, DictionaryDatum& dict, const size_t dendritic_delay_id = 0 ) const = 0;

  /**
   * Set status of the connection at position lcid according to the
   * dictionary dict.
   */
  virtual void set_synapse_status( const index lcid,
    const DictionaryDatum& dict,
    const ConnectorModel& cm,
    const size_t dendritic_delay_id = 0 ) = 0;

  /**
   * Get the proportion of the transmission delay attributed to the dendrite of a connection.
   */
  virtual delay get_dendritic_delay( const size_t dendritic_delay_id ) const = 0;

  /**
   * Get the time of the last pre-synaptic spike that was emitted over the specified connection.
   */
  virtual double get_last_presynaptic_spike( const index lcid, const size_t dendritic_delay_id ) const = 0;

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

  virtual long get_connection_label( const index lcid, const size_t dendritic_delay_id = 0 ) const = 0;

  virtual bool is_connection_disabled( const index lcid, const size_t dendritic_delay_id = 0 ) const = 0;

  /**
   * Send the event e to the connection at position lcid. Return bool
   * indicating whether the following connection belongs to the same
   * source.
   */
  virtual void send( const thread tid,
    const index target_node_id,
    const index lcid,
    const delay total_delay,
    const ConnectorModel* cm,
    Event& e ) = 0;
  virtual void send( const thread tid,
    const index target_node_id,
    const index lcid,
    const size_t dendritic_delay_id,
    const double tau_minus_inv,
    const std::deque< double >& history,
    const ConnectorModel* cm,
    Event& e ) = 0;

  virtual void
  send_weight_event( const thread tid, Event& e, const CommonSynapseProperties& cp, const index target_node_id ) = 0;

  virtual void update_stdp_connections( const index node_id,
    const thread tid,
    const double post_spike_time_syn,
    const delay dendritic_delay,
    const ConnectorModel* cm ) = 0;

  virtual void process_post_synaptic_spike( const index lcid,
    const double post_spike_time_syn,
    const ConnectorModel* cm,
    const size_t dendritic_delay_id = 0 ) = 0;

  virtual void
  update_trace( const double post_spike_time, const delay dendritic_delay, const double tau_minus_inv ) = 0;

  virtual double get_trace( const double pre_spike_time,
    const size_t dendritic_delay_id,
    const double tau_minus_inv,
    const std::deque< double >& history ) = 0;

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
  // TODO JV (pt): Structural plasticity
  // std::vector< index > sources_;

public:
  explicit Connector( const synindex syn_id, const size_t initial_size )
    : syn_id_( syn_id )
  {
    C_.reserve( initial_size );
  }

  ~Connector() override
  {
    std::vector< ConnectionT >().swap( C_ );
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
  get_synapse_status( const index lcid, DictionaryDatum& dict, const size_t ) const override
  {
    assert( lcid < C_.size() );

    C_[ lcid ].get_status( dict );

    // get target node ID here, where tid is available
    // necessary for hpc synapses using TargetIdentifierIndex
    // def< long >( dict, names::target, C_[ lcid ].get_target( tid )->get_node_id() );
  }

  void
  set_synapse_status( const index lcid, const DictionaryDatum& dict, const ConnectorModel& cm, const size_t ) override
  {
    assert( lcid < C_.size() );

    C_[ lcid ].set_status( dict, static_cast< const GenericConnectorModel< ConnectionT >& >( cm ) );
  }

  delay
  get_dendritic_delay( const size_t ) const override
  {
    throw UnexpectedEvent( "This connector type does not support retrieving dendritic delays." );
  }

  double
  get_last_presynaptic_spike( const index lcid, const size_t ) const override
  {
    return C_[ lcid ].get_last_presynaptic_spike();
  }

  index
  get_source( const index lcid ) override
  {
    assert( false ); // TODO JV (pt)
    // return sources_[ lcid ];
  }

  std::vector< index >
  get_sources() override
  {
    assert( false ); // TODO JV (pt)
    // return sources_;
  }

  index
  add_connection( const ConnectionT& c )
  {
    C_.push_back( c );
    if ( C_.size() > MAX_LOCAL_CONNECTION_ID )
    {
      throw KernelException(
        String::compose( "Too many connections: at most %1 connections supported per virtual "
                         "process and synapse model to a specific target neuron.",
          MAX_LOCAL_CONNECTION_ID ) );
    }
    return C_.size() - 1;
  }

  std::vector< index >
  get_connection_indices( const index source_node_id,
    const long connection_label = UNLABELED_CONNECTION ) const override
  {
    // TODO JV (pt): Structural plasticity
    //    // binary search in sorted sources
    //    std::vector< index >::const_iterator it = sources_.cbegin();
    //    const std::vector< index >::const_iterator end = sources_.cend();
    //
    std::vector< index > indices;
    //    while ( ( it = std::find( it, end, source_node_id ) ) != sources_.end() )
    //    {
    //      const index lcid = std::distance( sources_.cbegin(), it );
    //      // Connection is disabled
    //      if ( not C_[ lcid ].is_disabled() )
    //      {
    //        if ( connection_label == UNLABELED_CONNECTION or C_[ lcid ].get_label() == connection_label )
    //        {
    //          indices.push_back( lcid );
    //        }
    //      }
    //      ++it;
    //    }
    //
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
  get_connection_label( const index lcid, const size_t ) const override
  {
    return C_[ lcid ].get_label();
  }

  bool
  is_connection_disabled( const index lcid, const size_t ) const override
  {
    return C_[ lcid ].is_disabled();
  }

  void send( const thread tid,
    const index target_node_id,
    const index lcid,
    const delay total_delay,
    const ConnectorModel* cm,
    Event& e ) override;

  void
  send( const thread tid,
    const index target_node_id,
    const index lcid,
    const size_t dendritic_delay_id,
    const double tau_minus_inv,
    const std::deque< double >& history,
    const ConnectorModel* cm,
    Event& e ) override
  {
    throw UnexpectedEvent( "" );
  }

  void send_weight_event( const thread tid,
    Event& e,
    const CommonSynapseProperties& cp,
    const index target_node_id ) override;

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
  update_stdp_connections( const index node_id,
    const thread,
    const double post_spike_time_syn,
    const delay dendritic_delay,
    const ConnectorModel* cm ) override
  {
    throw UnexpectedEvent( "This connector type does not support STDP connections." );
  }

  void
  process_post_synaptic_spike( const index lcid,
    const double post_spike_time_syn,
    const ConnectorModel* cm,
    const size_t ) override
  {
    throw UnexpectedEvent( "This connector type does not support STDP connections." );
  }

  void
  update_trace( const double post_spike_time, const delay dendritic_delay, const double tau_minus_inv ) override
  {
    throw UnexpectedEvent( "This connector type does not support STDP connections." );
  }

  double
  get_trace( const double pre_spike_time,
    const size_t dendritic_delay_id,
    const double tau_minus_inv,
    const std::deque< double >& history ) override
  {
    throw UnexpectedEvent( "This connector type does not support STDP connections." );
  }

  void
  clear_history() override
  {
    throw UnexpectedEvent( "This connector type does not support STDP connections." );
  }

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

/**
 * Homogeneous connector, contains synapses of one particular type (syn_id_), additionally storing dendritic delays for
 * each connection.
 */
template < typename ConnectionT >
class DendriticDelayConnector : public ConnectorBase
{
  //! A region of connections in the Connector where all connections share the same dendritic delay
  struct DelayRegion
  {
    delay dendritic_delay; //!< Dendritic delay of this region
    std::vector< ConnectionT > connections;
    //! the post-synaptic trace in synapse time, i.e. the trace of the neuron "dendritic delay"-milliseconds ago
    double Kminus;
    double last_post_spike; //!< time of the last update in synapse time

    DelayRegion( delay dendritic_delay )
      : dendritic_delay( dendritic_delay )
      , Kminus( 0 )
      , last_post_spike( -1 )
    {
    }
  };

private:
  const synindex syn_id_;
  const size_t initial_size_;

  // TODO JV (pt): There will usually be very few elements in the map and no insertions at runtime, so using an ordered
  //  vector can speed up performance in this case. However, this might actually yield worse performance when having
  //  many different delays in the network. Requires benchmarking.
  std::vector< DelayRegion > dendritic_delay_regions_;
  std::map< delay, size_t > dendritic_delay_region_indices_;
  // std::map< delay, DelayRegion > dendritic_delay_regions_;

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
  // std::vector< index > sources_;

public:
  explicit DendriticDelayConnector( const synindex syn_id, const size_t initial_size )
    : syn_id_( syn_id )
    , initial_size_( initial_size )
  {
  }

  ~DendriticDelayConnector() override
  {
    std::vector< DelayRegion >().swap( dendritic_delay_regions_ );
  }

  synindex
  get_syn_id() const override
  {
    return syn_id_;
  }

  size_t
  size() const override
  {
    return std::accumulate( dendritic_delay_regions_.begin(),
      dendritic_delay_regions_.end(),
      0,
      []( const size_t sum, const DelayRegion& r ) -> size_t { return sum + r.connections.size(); } );
  }

  void
  get_synapse_status( const index lcid, DictionaryDatum& dict, const size_t dendritic_delay_id ) const override
  {
    assert( dendritic_delay_id < dendritic_delay_regions_.size() );
    assert( lcid < dendritic_delay_regions_[ dendritic_delay_id ].connections.size() );

    dendritic_delay_regions_[ dendritic_delay_id ].connections[ lcid ].get_status( dict );
  }

  void
  set_synapse_status( const index lcid,
    const DictionaryDatum& dict,
    const ConnectorModel& cm,
    const size_t dendritic_delay_id ) override
  {
    assert( dendritic_delay_id < dendritic_delay_regions_.size() );
    assert( lcid < dendritic_delay_regions_[ dendritic_delay_id ].connections.size() );

    dendritic_delay_regions_[ dendritic_delay_id ].connections[ lcid ].set_status(
      dict, static_cast< const GenericConnectorModel< ConnectionT >& >( cm ) );
  }

  delay
  get_dendritic_delay( const size_t dendritic_delay_id ) const override
  {
    assert( dendritic_delay_id < dendritic_delay_regions_.size() );

    return dendritic_delay_regions_[ dendritic_delay_id ].dendritic_delay;
  }

  double
  get_last_presynaptic_spike( const index lcid, const size_t dendritic_delay_id ) const override
  {
    return dendritic_delay_regions_[ dendritic_delay_id ].connections[ lcid ].get_last_presynaptic_spike();
  }

  index
  get_source( const index lcid ) override
  {
    assert( false ); // TODO JV (pt)
    // return sources_[ lcid ];
  }

  std::vector< index >
  get_sources() override
  {
    assert( false ); // TODO JV (pt)
    // return sources_;
  }

  std::pair< index, size_t >
  add_connection( const ConnectionT& c, const delay dendritic_delay )
  {
    size_t dendritic_delay_id = 0;
    bool delay_region_exists = false;
    for ( DelayRegion& region : dendritic_delay_regions_ )
    {
      if ( dendritic_delay == region.dendritic_delay )
      {
        region.connections.push_back( c );
        delay_region_exists = true;
        break;
      }
      ++dendritic_delay_id;
    }
    if ( not delay_region_exists )
    {
      dendritic_delay_region_indices_[ dendritic_delay ] = dendritic_delay_regions_.size();
      dendritic_delay_regions_.emplace_back( dendritic_delay );
      dendritic_delay_regions_[ dendritic_delay_id ].connections.reserve( initial_size_ );
      dendritic_delay_regions_[ dendritic_delay_id ].connections.push_back( c );
    }

    if ( dendritic_delay_regions_[ dendritic_delay_id ].connections.size() > MAX_LOCAL_CONNECTION_ID )
    {
      throw KernelException(
        String::compose( "Too many connections: at most %1 connections supported per virtual process and synapse model"
                         " to a specific target neuron for specific dendritic delay.",
          MAX_LOCAL_CONNECTION_ID ) );
    }

    return { dendritic_delay_regions_[ dendritic_delay_id ].connections.size() - 1, dendritic_delay_id };
  }

  std::vector< index >
  get_connection_indices( const index source_node_id,
    const long connection_label = UNLABELED_CONNECTION ) const override
  {
    // TODO JV (pt): Structural plasticity
    //    // binary search in sorted sources
    //    std::vector< index >::const_iterator it = sources_.cbegin();
    //    const std::vector< index >::const_iterator end = sources_.cend();
    //
    std::vector< index > indices;
    //    while ( ( it = std::find( it, end, source_node_id ) ) != sources_.end() )
    //    {
    //      const index lcid = std::distance( sources_.cbegin(), it );
    //      // Connection is disabled
    //      if ( not C_[ lcid ].is_disabled() )
    //      {
    //        if ( connection_label == UNLABELED_CONNECTION or C_[ lcid ].get_label() == connection_label )
    //        {
    //          indices.push_back( lcid );
    //        }
    //      }
    //      ++it;
    //    }
    //
    return indices;
  }

  std::vector< index >
  get_connection_indices( const long connection_label = UNLABELED_CONNECTION ) const override
  {
    // TODO JV (pt): Structural plasticity
    std::vector< index > indices;
    //    for ( auto conn_it = C_.cbegin(); conn_it != C_.cend(); ++conn_it )
    //    {
    //      const index lcid = std::distance( C_.cbegin(), conn_it );
    //      // Connection is disabled
    //      if ( not C_[ lcid ].is_disabled() )
    //      {
    //        if ( connection_label == UNLABELED_CONNECTION or C_[ lcid ].get_label() == connection_label )
    //        {
    //          indices.push_back( lcid );
    //        }
    //      }
    //    }
    return indices;
  }

  long
  get_connection_label( const index lcid, const size_t dendritic_delay_id ) const override
  {
    return dendritic_delay_regions_[ dendritic_delay_id ].connections[ lcid ].get_label();
  }

  bool
  is_connection_disabled( const index lcid, const size_t dendritic_delay_id ) const override
  {
    return dendritic_delay_regions_[ dendritic_delay_id ].connections[ lcid ].is_disabled();
  }
  void
  send( const thread tid,
    const index target_node_id,
    const index lcid,
    const delay total_delay,
    const ConnectorModel* cm,
    Event& e ) override
  {
    throw UnexpectedEvent( "" );
  }

  void send( const thread tid,
    const index target_node_id,
    const index lcid,
    const size_t dendritic_delay_id,
    const double tau_minus_inv,
    const std::deque< double >& history,
    const ConnectorModel* cm,
    Event& e ) override;

  void send_weight_event( const thread tid,
    Event& e,
    const CommonSynapseProperties& cp,
    const index target_node_id ) override;

  void
  trigger_update_weight( const long vt_node_id,
    const thread tid,
    const std::vector< spikecounter >& dopa_spikes,
    const double t_trig,
    const std::vector< ConnectorModel* >& cm ) override
  {
    // TODO JV (pt): Dopa synapses
    //    for ( size_t i = 0; i < C_.size(); ++i )
    //    {
    //      if ( static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )
    //             ->get_common_properties()
    //             .get_vt_node_id()
    //        == vt_node_id )
    //      {
    //        C_[ i ].trigger_update_weight( tid,
    //          dopa_spikes,
    //          t_trig,
    //          get_dendritic_delay( i ),
    //          static_cast< GenericConnectorModel< ConnectionT >* >( cm[ syn_id_ ] )->get_common_properties() );
    //      }
    //    }
  }

  void update_stdp_connections( const index node_id,
    const thread tid,
    const double post_spike_time_syn,
    const delay dendritic_delay,
    const ConnectorModel* cm ) override;

  void
  process_post_synaptic_spike( const index lcid,
    const double post_spike_time_syn,
    const ConnectorModel* cm,
    const size_t dendritic_delay_id ) override
  {
    typename ConnectionT::CommonPropertiesType const& cp =
      static_cast< const GenericConnectorModel< ConnectionT >* >( cm )->get_common_properties();
    dendritic_delay_regions_[ dendritic_delay_id ].connections[ lcid ].process_post_synaptic_spike(
      post_spike_time_syn, cp );
  }

  void
  update_trace( const double post_spike_time, const delay dendritic_delay, const double tau_minus_inv ) override
  {
    const auto delay_region_it = dendritic_delay_region_indices_.find( dendritic_delay );

    if ( delay_region_it != dendritic_delay_region_indices_.end() )
    {
      DelayRegion& delay_region = dendritic_delay_regions_[ delay_region_it->second ];
      // update post-synaptic trace
      delay_region.Kminus = delay_region.Kminus
          * std::exp(
            ( delay_region.last_post_spike - ( post_spike_time + Time::delay_steps_to_ms( dendritic_delay ) ) )
            * tau_minus_inv )
        + 1;
      delay_region.last_post_spike = post_spike_time + Time::delay_steps_to_ms( dendritic_delay );
    }
  }

  double get_trace( const double pre_spike_time,
    const size_t dendritic_delay_id,
    const double tau_minus_inv,
    const std::deque< double >& history ) override;

  void
  clear_history() override
  {
    for ( DelayRegion& delay_region : dendritic_delay_regions_ )
    {
      delay_region.Kminus = 0;
      delay_region.last_post_spike = -1;
    }
  }

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