/*
 *  node.h
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

#ifndef NODE_H
#define NODE_H

// C++ includes:
#include <deque>
#include <string>
#include <vector>

// Includes from nestkernel:
#include "archived_spike.h"
#include "connection_type_enum.h"
#include "connector_base.h"
#include "deprecation_warning.h"
#include "event.h"
#include "nest_names.h"
#include "nest_time.h"
#include "nest_types.h"
#include "node_collection.h"
#include "secondary_event.h"

// Includes from sli:
#include "dictdatum.h"

/** @file node.h
 * Declarations for base class Node
 */

namespace nest
{
class ConnectorModel;
class Model;
class TimeConverter;


/**
 * @defgroup user_interface Model developer interface.
 * Functions and classes important for implementing new Node and
 * Model classes.
 */

/**
 * Base class for all NEST network objects.
 *
 * Class Node is the top of the simulation object hierarchy. It
 * defines the most general interface to a network element.
 *
 * Class Node provide the interface for
 * - updating the dynamic state of an object
 * - connecting nodes, using particular Events
 * - accepting connection requests
 * - handling incoming events
 * A new type of Node must be derived from this base class and
 * implement its interface.
 * In order to keep the inheritance hierarchy flat, it is encouraged
 * to directly subclass from base class Node.
 *
 * @see class Event
 * @ingroup user_interface
 */

/** @BeginDocumentation

   Name: Node - General properties of all nodes.

   Parameters:
   frozen     booltype    - Whether the node is updated during simulation
   global_id  integertype - The node ID of the node (cf. local_id)
   local      booltype    - Whether the node is available on the local process
   model      literaltype - The model type the node was created from
   state      integertype - The state of the node (see the help on elementstates
                            for details)
   thread     integertype - The id of the thread the node is assigned to (valid
                            locally)
   vp         integertype - The id of the virtual process the node is assigned
                            to (valid globally)

   SeeAlso: GetStatus, SetStatus, elementstates
 */

class Node
{
  friend class NodeManager;
  friend class ModelManager;
  friend class proxynode;
  friend class Model;
  friend class SimulationManager;

  Node& operator=( const Node& ); //!< not implemented

public:
  Node();
  Node( Node const& );
  virtual ~Node();

  /**
   * Virtual copy constructor.
   * This function should create a new object by
   * calling the derived class' copy constructor and
   * return its pointer.
   */
  virtual Node*
  clone() const
  {
    return nullptr;
  }

  /**
   * Returns true if the node has proxies on remote threads. This is
   * used to discriminate between different types of nodes, when adding
   * new nodes to the network.
   */
  virtual bool has_proxies() const;

  /**
   * Returns true if the node supports the Urbanczik-Senn plasticity rule
   */
  virtual bool supports_urbanczik_archiving() const;

  /**
   * Base Node class does not support postponed delivery. Postponed delivery has to be implemented in derived classes if
   * they need to support this feature.
   */
  inline virtual bool
  supports_postponed_delivery() const
  {
    return false;
  }

  /**
   * Returns true if the node only receives events from nodes/devices
   * on the same thread.
   */
  virtual bool local_receiver() const;

  /**
   * Returns true if the node exists only once per process, but does
   * not have proxies on remote threads. This is used to
   * discriminate between different types of nodes, when adding new
   * nodes to the network.
   *
   * TODO: Is this true for *any* model at all? Maybe MUSIC related?
   */
  virtual bool one_node_per_process() const;

  /**
   * Returns true if the node sends/receives off-grid events. This is
   * used to discriminate between different types of nodes when adding
   * new nodes to the network.
   */
  virtual bool is_off_grid() const;

  /**
   * Returns true if the node is a proxy node. This is implemented because
   * the use of RTTI is rather expensive.
   */
  virtual bool is_proxy() const;

  /**
   * Return class name.
   * Returns name of node model (e.g. "iaf_psc_alpha") as string.
   * This name is identical to the name that is used to identify
   * the model in the interpreter's model dictionary.
   */
  std::string get_name() const;

  /**
   * Return the element type of the node.
   * The returned Name is a free label describing the class of network
   * elements a node belongs to. Currently used values are "neuron",
   * "recorder", "stimulator", and "other", which are all defined as
   * static Name objects in the names namespace.
   * This function is overwritten with a corresponding value in the
   * derived classes
   */
  virtual Name get_element_type() const;

  /**
   * Return global Network ID.
   * Returns the global network ID of the Node.
   * Each node has a unique network ID which can be used to access
   * the Node comparable to a pointer.
   *
   * The smallest valid node ID is 1.
   */
  index get_node_id() const;

  /**
   * Return lockpointer to the NodeCollection that created this node.
   */
  NodeCollectionPTR get_nc() const;

  /**
   * Return model ID of the node.
   * Returns the model ID of the model for this node.
   * Model IDs start with 0.
   * @note The model ID is not stored in the model prototype instance.
   *       It is only set when actual nodes are created from a prototype.
   */
  int get_model_id() const;

  /**
   * Prints out one line of the tree view of the network.
   */
  virtual std::string
  print_network( int, int, std::string = "" )
  {
    return std::string();
  }

  /**
   * Returns true if node is frozen, i.e., shall not be updated.
   */
  bool is_frozen() const;

  /**
   * Returns true if the node uses the waveform relaxation method
   */
  bool node_uses_wfr() const;

  /**
   * Sets node_uses_wfr_ member variable
   * (to be able to set it to "true" for any class derived from Node)
   */
  void set_node_uses_wfr( const bool );

  /**
   * Initialize node prior to first simulation after node has been created.
   *
   * init() allows the node to configure internal data structures prior to
   * being simulated. The method has an effect only the first time it is
   * called on a given node, otherwise it returns immediately. init() calls
   * virtual functions init_state_() and init_buffers_().
   */
  virtual void init();

  /**
   * Re-calculate dependent parameters of the node.
   * This function is called each time a simulation is begun/resumed.
   * It must re-calculate all internal Variables of the node required
   * for spike handling or updating the node.
   *
   */
  virtual void pre_run_hook() = 0;

  /**
   * Re-calculate time-based properties of the node.
   * This function is called after a change in resolution.
   */
  virtual void
  calibrate_time( const TimeConverter& )
  {
  }

  /**
   * Cleanup node after Run. Override this function if a node needs to
   * "wrap up" things after a call to Run, i.e., before
   * SimulationManager::run() returns. Typical use-cases are devices
   * that need to flush buffers.
   */
  virtual void
  post_run_cleanup()
  {
  }

  /**
   * Finalize node.
   * Override this function if a node needs to "wrap up" things after a full simulation, i.e., a cycle of Prepare, Run,
   * Cleanup. Typical use-cases are devices that need to close files.
   */
  virtual void finalize();

  /**
   * Bring the node from state $t$ to $t+n*dt$.
   *
   * n->update(T, from, to) performs the update steps beginning at T+from .. T+to-1, ie, emitting events with time
   * stamps T+from+1 .. T+to.
   *
   * @param Time   network time at beginning of time slice.
   * @param long initial step inside time slice
   * @param long post-final step inside time slice
   *
   */
  virtual void update( const Time&, const long, const long ) = 0;

  /**
   * Prepare the node for the next update cycle.
   */
  virtual void
  prepare_update( const Time origin )
  {
  }

  /**
   * Bring the node from state $t$ to $t+n*dt$, sends SecondaryEvents
   * (e.g. GapJunctionEvent) and resets state variables to values at $t$.
   *
   * n->wfr_update(T, from, to) performs the update steps beginning
   * at T+from .. T+to-1.
   *
   * Does not emit spikes, does not log state variables.
   *
   * throws UnexpectedEvent if not reimplemented in derived class
   *
   * @param Time   network time at beginning of time slice.
   * @param long initial step inside time slice
   * @param long post-final step inside time slice
   *
   */
  virtual bool wfr_update( const Time&, const long, const long );

  /**
   * @defgroup status_interface Configuration interface.
   * Functions and infrastructure, responsible for the configuration
   * of Nodes from the SLI Interpreter level.
   *
   * Each node can be configured from the SLI level through a named
   * parameter interface. In order to change parameters, the user
   * can specify name value pairs for each parameter. These pairs
   * are stored in a data structure which is called Dictionary.
   * Likewise, the user can query the configuration of any node by
   * requesting a dictionary with name value pairs.
   *
   * The configuration interface consists of four functions which
   * implement storage and retrieval of named parameter sets.
   */

  /**
   * Change properties of the node according to the
   * entries in the dictionary.
   * @param d Dictionary with named parameter settings.
   * @ingroup status_interface
   */
  virtual void set_status( const DictionaryDatum& ) = 0;

  /**
   * Export properties of the node by setting
   * entries in the status dictionary.
   * @param d Dictionary.
   * @ingroup status_interface
   */
  virtual void get_status( DictionaryDatum& ) const = 0;

  DictionaryDatum get_connection_status( const synindex syn_id, const index lcid, DictionaryDatum& dict ) const;

  DictionaryDatum get_connection_status( const synindex syn_id,
    const index lcid,
    const size_t dendritic_delay_id,
    DictionaryDatum& dict ) const;

  void set_connection_status( const synindex syn_id,
    const index lcid,
    const DictionaryDatum& dict,
    const ConnectorModel& cm );

  void set_connection_status( const synindex syn_id,
    const index lcid,
    const size_t dendritic_delay_id,
    const DictionaryDatum& dict,
    ConnectorModel& cm );


public:
  /**
   * @defgroup event_interface Communication.
   * Functions and infrastructure, responsible for communication
   * between Nodes.
   *
   * Nodes communicate by sending an receiving events. The
   * communication interface consists of two parts:
   * -# Functions to handle incoming events.
   * -# Functions to check if a connection between nodes is possible.
   *
   * @see Event
   */

  /**
   * Send an event to the receiving_node passed as an argument.
   * This is required during the connection handshaking to test,
   * if the receiving_node can handle the event type and receptor_type sent
   * by the source node.
   *
   * // TODO JV (pt): Still needed?
   * If dummy_target is true, this indicates that receiving_node is derived from
   * ConnTestDummyNodeBase and used in the first call to send_test_event().
   * This can be ignored in most cases, but Nodes sending DS*Events to their
   * own event hooks and then *Events to their proper targets must send
   * DS*Events when called with the dummy target, and *Events when called with
   * the real target, see #478.
   */
  virtual port send_test_event( Node& receiving_node, const rport receptor_type, synindex syn_id );

  /**
   * Check if the node can handle a particular event and receptor type.
   * This function is called upon connection setup by send_test_event().
   *
   * handles_test_event() function is used to verify that the receiver
   * can handle the event. It can also be used by the receiver to
   * return information to the sender in form of the returned port.
   * The default implementation throws an IllegalConnection
   * exception.  Any node class should define handles_test_event()
   * functions for all those event types it can handle.
   *
   * See Kunkel et al, Front Neuroinform 8:78 (2014), Sec 3.
   *
   * @note The semantics of all other handles_test_event() functions is
   * identical.
   * @ingroup event_interface
   * @throws IllegalConnection
   */
  virtual port handles_test_event( SpikeEvent&, rport receptor_type );
  virtual port handles_test_event( WeightRecorderEvent&, rport receptor_type );
  virtual port handles_test_event( RateEvent&, rport receptor_type );
  virtual port handles_test_event( DataLoggingRequest&, rport receptor_type );
  virtual port handles_test_event( CurrentEvent&, rport receptor_type );
  virtual port handles_test_event( ConductanceEvent&, rport receptor_type );
  virtual port handles_test_event( DoubleDataEvent&, rport receptor_type );
  virtual port handles_test_event( GapJunctionEvent&, rport receptor_type );
  virtual port handles_test_event( InstantaneousRateConnectionEvent&, rport receptor_type );
  virtual port handles_test_event( DiffusionConnectionEvent&, rport receptor_type );
  virtual port handles_test_event( DelayedRateConnectionEvent&, rport receptor_type );

  /**
   * Required to check, if source neuron may send a SecondaryEvent.
   * This base class implementation throws IllegalConnection
   * and needs to be overwritten in the derived class.
   * @ingroup event_interface
   * @throws IllegalConnection
   */
  virtual void sends_secondary_event( GapJunctionEvent& ge );

  /**
   * Required to check, if source neuron may send a SecondaryEvent.
   * This base class implementation throws IllegalConnection
   * and needs to be overwritten in the derived class.
   * @ingroup event_interface
   * @throws IllegalConnection
   */
  virtual void sends_secondary_event( InstantaneousRateConnectionEvent& re );

  /**
   * Required to check, if source neuron may send a SecondaryEvent.
   * This base class implementation throws IllegalConnection
   * and needs to be overwritten in the derived class.
   * @ingroup event_interface
   * @throws IllegalConnection
   */
  virtual void sends_secondary_event( DiffusionConnectionEvent& de );

  /**
   * Required to check, if source neuron may send a SecondaryEvent.
   * This base class implementation throws IllegalConnection
   * and needs to be overwritten in the derived class.
   * @ingroup event_interface
   * @throws IllegalConnection
   */
  virtual void sends_secondary_event( DelayedRateConnectionEvent& re );

  /**
   * Register a STDP connection
   *
   * @throws IllegalConnection
   */
  virtual void register_stdp_connection( const delay, const delay, const synindex );

  /**
   * Change the number of different connection types to this node.
   */
  void
  resize_connections( const size_t size )
  {
    connections_.resize( size );
    connections_from_devices_.resize( size );
  }

  /**
   * Get the number of connections to this neuron of a specific connection type.
   */
  size_t
  get_num_conn_type_sources( const synindex syn_id ) const
  {
    if ( connections_[ syn_id ] )
    {
      return connections_[ syn_id ]->size();
    }
    return 0;
  }

  /**
   * Get the total number of connections to this node.
   */
  size_t
  get_num_connections() const
  {
    return std::accumulate( connections_.cbegin(),
      connections_.cend(),
      0,
      []( size_t sum, auto& sources_syn_id )
      {
        if ( sources_syn_id )
        {
          return sum + sources_syn_id->size();
        }
        else
        {
          return sum;
        }
      } );
  }

  /**
   * Get information about the source node of a specific connection.
   */
  index
  get_source( const synindex syn_id, const index local_connection_id )
  {
    assert( connections_[ syn_id ] );

    return connections_[ syn_id ]->get_source( local_connection_id );
  }

  std::vector< index >
  get_sources( const synindex syn_id )
  {
    assert( connections_[ syn_id ] );

    return connections_[ syn_id ]->get_sources();
  }

  index
  get_source_from_devices( const synindex syn_id, const index local_connection_id )
  {
    assert( connections_[ syn_id ] );

    return connections_[ syn_id ]->get_source( local_connection_id );
  }

  std::vector< index >
  get_sources_from_devices( const synindex syn_id )
  {
    assert( connections_[ syn_id ] );

    return connections_[ syn_id ]->get_sources();
  }

  void
  get_connections( std::deque< ConnectionID > conns,
    const synindex syn_id,
    const index source_node_id = DISABLED_NODE_ID,
    const long connection_label = UNLABELED_CONNECTION )
  {
    for ( const index lcid : get_connection_indices( syn_id, source_node_id, connection_label ) )
    {
      conns.push_back( ConnectionDatum(
        ConnectionID( connections_[ syn_id ]->get_source( lcid ), node_id_, thread_, syn_id, lcid ) ) );
    }
  }

  std::vector< index >
  get_connection_indices( const synindex syn_id,
    const index source_node_id = DISABLED_NODE_ID,
    const long connection_label = UNLABELED_CONNECTION ) const
  {
    if ( connections_[ syn_id ] )
    {
      if ( source_node_id != DISABLED_NODE_ID )
      {
        return connections_[ syn_id ]->get_connection_indices( source_node_id, connection_label );
      }
      else
      {
        return connections_[ syn_id ]->get_connection_indices( connection_label );
      }
    }
    return std::vector< index >();
  }

  void
  get_connections_from_devices( std::deque< ConnectionID > conns,
    const synindex syn_id,
    const index source_node_id = DISABLED_NODE_ID,
    const long connection_label = UNLABELED_CONNECTION )
  {
    for ( const index lcid : get_connection_indices_from_devices( syn_id, source_node_id, connection_label ) )
    {
      conns.push_back( ConnectionDatum(
        ConnectionID( connections_from_devices_[ syn_id ]->get_source( lcid ), node_id_, thread_, syn_id, lcid ) ) );
    }
  }

  std::vector< index >
  get_connection_indices_from_devices( const synindex syn_id,
    const index source_node_id = DISABLED_NODE_ID,
    const long connection_label = UNLABELED_CONNECTION ) const
  {
    if ( connections_from_devices_[ syn_id ] )
    {
      if ( source_node_id != DISABLED_NODE_ID )
      {
        return connections_from_devices_[ syn_id ]->get_connection_indices( source_node_id, connection_label );
      }
      else
      {
        return connections_from_devices_[ syn_id ]->get_connection_indices( connection_label );
      }
    }
    return std::vector< index >();
  }

  /**
   * Remove source information of all connections to this node.
   */
  void
  clear_sources()
  {
    for ( std::unique_ptr< ConnectorBase >& connections_per_syn_type : connections_ )
    {
      if ( connections_per_syn_type )
      {
        // connections_per_syn_type->clear_sources();  // TODO JV (pt)
      }
    }
    for ( std::unique_ptr< ConnectorBase >& connections_per_syn_type : connections_from_devices_ )
    {
      if ( connections_per_syn_type )
      {
        // connections_per_syn_type->clear_sources();  // TODO JV (pt)
      }
    }
  }

  /**
   * Disables a connection to this node, without removing it.
   */
  virtual void disable_connection( const synindex syn_id, const index local_connection_id );

  /**
   * Deletes all connections which are marked disabled.
   */
  virtual void remove_disabled_connections();

  /**
   * Completely removes a connection to this node.
   */
  void delete_connections();

  /**
   * This function calls check_connection() on the sender to check if the
   * receiver accepts the event type and receptor type requested by the sender.
   * \param source The source node
   * \param syn_id Connection type of the connection as index
   * \param receptor The ID of the requested receptor type
   */
  template < typename ConnectionT >
  void check_connection( Node& source, ConnectionT& connection, const synindex syn_id, const rport receptor_type,
    const delay total_delay, const typename ConnectionT::CommonPropertiesType& cp );

  /**
   * Adds a connection to this node of a specific connection type.
   */
  template < typename ConnectionT >
  const std::pair< index, size_t > add_connection( Node& source_node,
    const synindex syn_id,
    ConnectionT& connection,
    const rport receptor_type,
    const bool is_primary,
    const ConnectionType from_device,
    const delay axonal_delay,
    const delay dendritic_delay );

  /**
   * When receiving an event from a device, forward it to the corresponding connection and handle the event previously
   * updated by the connection.
   */
  template < typename EventT >
  void deliver_event_from_device( const thread tid,
    const synindex syn_id,
    const index local_target_connection_id,
    const delay d,
    const ConnectorModel* cm,
    EventT& e );

  /**
   * When receiving an incoming spike event, forward it to the corresponding connection and handle the event previously
   * updated by the connection.
   */
  virtual void deliver_event( const synindex syn_id,
    const index local_target_connection_id,
    const size_t dendritic_delay_id,
    const ConnectorModel* cm,
    const Time lag,
    const delay d,
    const double offset );

  /**
   * Handle incoming spike events.
   * @param e Event object.
   *
   * This handler has to be implemented if a Node should accept spike events.
   * @see class SpikeEvent
   * @ingroup event_interface
   */
  virtual void handle( SpikeEvent& e );

  /**
   * Handle incoming weight recording events.
   * @param e Event object.
   *
   * This handler has to be implemented if a Node should accept weight recording events.
   * @see class WeightRecordingEvent
   * @ingroup event_interface
   */
  virtual void handle( WeightRecorderEvent& e );

  /**
   * Handler for rate events.
   * @see handle(SpikeEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( RateEvent& e );

  /**
   * Handler for universal data logging request.
   * @see handle(SpikeEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( DataLoggingRequest& e );

  /**
   * Handler for universal data logging request.
   * @see handle(SpikeEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   * @note There is no connect_sender() for DataLoggingReply, since
   *       this event is only used as "back channel" for DataLoggingRequest.
   */
  virtual void handle( DataLoggingReply& e );

  /**
   * Handler for current events.
   * @see handle(thread, SpikeEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( CurrentEvent& e );

  /**
   * Handler for conductance events.
   * @see handle(thread, SpikeEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( ConductanceEvent& e );

  /**
   * Handler for DoubleData events.
   * @see handle(thread, SpikeEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( DoubleDataEvent& e );

  /**
   * Handler for gap junction events.
   * @see handle(thread, GapJunctionEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( GapJunctionEvent& e );

  /**
   * Handler for rate neuron events.
   * @see handle(thread, InstantaneousRateConnectionEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( InstantaneousRateConnectionEvent& e );

  /**
   * Handler for rate neuron events.
   * @see handle(thread, InstantaneousRateConnectionEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( DiffusionConnectionEvent& e );

  /**
   * Handler for delay rate neuron events.
   * @see handle(thread, DelayedRateConnectionEvent&)
   * @ingroup event_interface
   * @throws UnexpectedEvent
   */
  virtual void handle( DelayedRateConnectionEvent& e );

  // TODO JV (pt): Remove this again, is needed for devices to work, but should be cleaned up
  virtual void handle( SecondaryEvent& e );

  /**
   * @defgroup SP_functions Structural Plasticity in NEST.
   * Functions related to accessibility and setup of variables required for
   * the implementation of a model of Structural Plasticity in NEST.
   *
   */

  /**
   * Return the Ca_minus value at time Ca_t which corresponds to the time of
   * the last update in Calcium concentration which is performed each time
   * a Node spikes.
   * Return 0.0 if not overridden
   * @ingroup SP_functions
   */
  virtual double
  get_Ca_minus() const
  {
    return 0.0;
  }

  /**
   * Get the number of synaptic element for the current Node at Ca_t which
   * corresponds to the time of the last spike.
   * Return 0.0 if not overridden
   * @ingroup SP_functions
   */
  virtual double
  get_synaptic_elements( Name ) const
  {
    return 0.0;
  }

  /**
   * Get the number of vacant synaptic element for the current Node
   * Return 0 if not overridden
   * @ingroup SP_functions
   */
  virtual int
  get_synaptic_elements_vacant( Name ) const
  {
    return 0;
  }

  /**
   * Get the number of connected synaptic element for the current Node
   * Return 0 if not overridden
   * @ingroup SP_functions
   */
  virtual int
  get_synaptic_elements_connected( Name ) const
  {
    return 0;
  }

  /**
   * Get the number of all synaptic elements for the current Node at time t
   * Return an empty map if not overridden
   * @ingroup SP_functions
   */
  virtual std::map< Name, double >
  get_synaptic_elements() const
  {
    return std::map< Name, double >();
  }

  /**
   * Triggers the update of all SynapticElements
   * stored in the synaptic_element_map_. It also updates the calcium
   * concentration.
   * @param t double time when the update is being performed
   * @ingroup SP_functions
   */
  virtual void update_synaptic_elements( double ) {};

  /**
   * Is used to reduce the number of synaptic elements in the node through
   * time. This amount is defined by tau_vacant.
   * @ingroup SP_functions
   */
  virtual void decay_synaptic_elements_vacant() {};

  /**
   * Is used to update the number of connected
   * synaptic elements (SynapticElement::z_connected_) when a synapse
   * is formed or deleted.
   * @param type Name, name of the synaptic element to connect
   * @param n int number of new connections of the given type
   * @ingroup SP_functions
   */
  virtual void connect_synaptic_element( Name, int ) {};

  virtual double get_trace( const double pre_spike_time, const double dendritic_delay, const synindex syn_id );

  /**
   * return the Kminus value at t (in ms).
   * @throws UnexpectedEvent
   */
  // TODO JV (pt): Remove and update models
  virtual double get_K_value( const double dendritic_delay, const double t, const synindex syn_id );

  virtual double get_LTD_value( double t );

  /**
   * write the Kminus, nearest_neighbor_Kminus, and Kminus_triplet
   * values at t (in ms) to the provided locations.
   * @throws UnexpectedEvent
   */
  virtual void get_K_values( double t, double& Kminus, double& nearest_neighbor_Kminus, double& Kminus_triplet );

  /**
   * return the spike history for (t1,t2].
   * @throws UnexpectedEvent
   */
  // TODO JV (pt): Remove and update models
  virtual void get_history( double t1,
    double t2,
    std::deque< ArchivedSpikeTrace >::iterator* start,
    std::deque< ArchivedSpikeTrace >::iterator* finish );

  // for Clopath synapse
  virtual void get_LTP_history( double t1,
    double t2,
    std::deque< ArchivedSpikeGeneric >::iterator* start,
    std::deque< ArchivedSpikeGeneric >::iterator* finish );
  // for Urbanczik synapse
  virtual void get_urbanczik_history( double t1,
    double t2,
    std::deque< ArchivedSpikeGeneric >::iterator* start,
    std::deque< ArchivedSpikeGeneric >::iterator* finish,
    int );
  // make neuron parameters accessible in Urbanczik synapse
  virtual double get_C_m( int comp );
  virtual double get_g_L( int comp );
  virtual double get_tau_L( int comp );
  virtual double get_tau_s( int comp );
  virtual double get_tau_syn_ex( int comp );
  virtual double get_tau_syn_in( int comp );

  /**
   * Store the number of the thread to which the node is assigned.
   * The assignment is done after node creation by the Network class.
   * @see: NodeManager::add_node().
   */
  void set_thread( thread );

  /**
   * Retrieve the number of the thread to which the node is assigned.
   */
  thread get_thread() const;

  /**
   * Store the number of the virtual process to which the node is assigned.
   * This is assigned to the node in NodeManager::add_node().
   */
  void set_vp( thread );

  /**
   * Retrieve the number of the virtual process to which the node is assigned.
   */
  thread get_vp() const;

  /** Set the model id.
   * This method is called by NodeManager::add_node() when a node is created.
   * @see get_model_id()
   */
  void set_model_id( int );

  /** Execute post-initialization actions in node models.
   * This method is called by NodeManager::add_node() on a node once
   * is fully initialized, i.e. after node ID, nc, model_id, thread, vp is
   * set.
   */
  void set_initialized();

  /**
   * @returns type of signal this node produces
   * used in check_connection to only connect neurons which send / receive
   * compatible information
   */
  virtual SignalType
  sends_signal() const
  {
    return SPIKE;
  }

  /**
   * @returns type of signal this node consumes
   * used in check_connection to only connect neurons which send / receive
   * compatible information
   */
  virtual SignalType
  receives_signal() const
  {
    return SPIKE;
  }

  /**
   *  Return a dictionary with the node's properties.
   *
   *  get_status_base() first gets a dictionary with the basic
   *  information of an element, using get_status_dict_(). It then
   *  calls the custom function get_status(DictionaryDatum) with
   *  the created status dictionary as argument.
   */
  DictionaryDatum get_status_base();

  /**
   * Set status dictionary of a node.
   *
   * Forwards to set_status() of the derived class.
   * @internal
   */
  void set_status_base( const DictionaryDatum& );

  /**
   * Returns true if node is model prototype.
   */
  bool is_model_prototype() const;

  /**
   * set thread local index

   */
  void set_thread_lid( const index );

  /**
   * get thread local index
   */
  index get_thread_lid() const;

  /**
   * Sets the local device id.
   * Throws an error if used on a non-device node.
   * @see get_local_device_id
   */
  virtual void set_local_device_id( const index lsdid );

  /**
   * Gets the local device id.
   * Throws an error if used on a non-device node.
   * @see set_local_device_id
   */
  virtual index get_local_device_id() const;

  /**
   * Return if the node has any incoming stdp connections.
   */
  virtual bool has_stdp_connections() const;

  /**
   * Member of DeprecationWarning class to be used by models if parameters are
   * deprecated.
   */
  DeprecationWarning deprecation_warning;

private:
  void set_node_id_( index ); //!< Set global node id

  void set_nc_( NodeCollectionPTR );

  /** Return a new dictionary datum .
   *
   * This function is called by get_status_base() and returns a new
   * empty dictionary by default.  Some nodes may contain a
   * permanent status dictionary which is then returned by
   * get_status_dict_().
   */
  virtual DictionaryDatum get_status_dict_();

protected:
  /**
   * Configure state variables depending on runtime information.
   *
   * Overload this method if the node needs to adapt state variables prior to
   * first simulation to runtime information, e.g., the number of incoming
   * connections.
   */
  virtual void init_state_();

  /**
   * Configure persistent internal data structures.
   *
   * Let node configure persistent internal data structures, such as input buffers or ODE solvers, to runtime
   * information prior to first simulation.
   */
  virtual void init_buffers_();

  virtual void set_initialized_();

  Model& get_model_() const;

  //! Mark node as frozen.
  void
  set_frozen_( bool frozen )
  {
    frozen_ = frozen;
  }

  /**
   * Auxiliary function to downcast a Node to a concrete class derived from Node.
   * @note This function is used to convert generic Node references to specific
   *       ones when intializing parameters or state from a prototype.
   */
  template < typename ConcreteNode >
  const ConcreteNode& downcast( const Node& );

private:
  /**
   * Model ID.
   * It is only set for actual node instances, not for instances of class Node
   * representing model prototypes. Model prototypes always have model_id_==-1.
   * @see get_model_id(), set_model_id()
   */
  int model_id_;

  thread thread_;      //!< thread node is assigned to
  thread vp_;          //!< virtual process node is assigned to
  bool frozen_;        //!< node shall not be updated if true
  bool initialized_;   //!< state and buffers have been initialized
  bool node_uses_wfr_; //!< node uses waveform relaxation method

  NodeCollectionPTR nc_ptr_;

protected:
  /**
   * Global Element ID (node ID).
   *
   * The node ID is unique within the network. The smallest valid node ID is 1.
   */
  index node_id_;

  /**
   * Local id of this node in the thread-local vector of nodes.
   */
  index thread_lid_;

  /**
   * A structure to hold the Connector objects which in turn hold the connection information of all incoming connections
   * to this node. Corresponds to a two dimensional structure:
   * synapse types|connections
   */
  std::vector< std::unique_ptr< ConnectorBase > >
    connections_; // TODO JV (pt): Only store entries of syn_ids that are in use

  /**
   * A structure to hold the Connector objects which in turn hold the connection information of all incoming connections
   * from devices to this node. Corresponds to a two dimensional structure:
   * synapse types|connections
   */
  std::vector< std::unique_ptr< ConnectorBase > > connections_from_devices_;
};

inline bool
Node::is_frozen() const
{
  return frozen_;
}

inline bool
Node::node_uses_wfr() const
{
  return node_uses_wfr_;
}

inline bool
Node::supports_urbanczik_archiving() const
{
  return false;
}

inline void
Node::set_node_uses_wfr( const bool uwfr )
{
  node_uses_wfr_ = uwfr;
}

inline bool
Node::has_proxies() const
{
  return true;
}

inline bool
Node::local_receiver() const
{
  return false;
}

inline bool
Node::one_node_per_process() const
{
  return false;
}

inline bool
Node::is_off_grid() const
{
  return false;
}

inline bool
Node::is_proxy() const
{
  return false;
}

inline Name
Node::get_element_type() const
{
  return names::neuron;
}

inline index
Node::get_node_id() const
{
  return node_id_;
}

inline NodeCollectionPTR
Node::get_nc() const
{
  return nc_ptr_;
}

inline void
Node::set_node_id_( index i )
{
  node_id_ = i;
}


inline void
Node::set_nc_( NodeCollectionPTR nc_ptr )
{
  nc_ptr_ = nc_ptr;
}

inline int
Node::get_model_id() const
{
  return model_id_;
}

inline void
Node::set_model_id( int i )
{
  model_id_ = i;
}

inline bool
Node::is_model_prototype() const
{
  return vp_ == invalid_thread;
}

inline void
Node::set_thread( thread t )
{
  thread_ = t;
}

inline thread
Node::get_thread() const
{
  return thread_;
}

inline void
Node::set_vp( thread vp )
{
  vp_ = vp;
}

inline thread
Node::get_vp() const
{
  return vp_;
}

template < typename ConcreteNode >
const ConcreteNode&
Node::downcast( const Node& n )
{
  ConcreteNode const* tp = dynamic_cast< ConcreteNode const* >( &n );
  assert( tp != 0 );
  return *tp;
}

inline void
Node::set_thread_lid( const index tlid )
{
  thread_lid_ = tlid;
}

inline index
Node::get_thread_lid() const
{
  return thread_lid_;
}

inline void
Node::deliver_event( const synindex syn_id,
  const index local_target_connection_id,
  const size_t dendritic_delay_id,
  const ConnectorModel* cm,
  const Time lag,
  const delay total_delay,
  const double offset )
{
  SpikeEvent se;
  se.set_stamp( lag );
  se.set_offset( offset ); // TODO JV (help): Why can't offset be incorporated into lag?
  se.set_sender_node_id_info( thread_, syn_id, node_id_, local_target_connection_id );

  // Send the event to the connection over which this event is transmitted to the node. The connection modifies the
  // event by adding a weight and optionally updates its internal state as well.
  connections_[ syn_id ]->send( thread_, node_id_, local_target_connection_id, total_delay, cm, se );

  // TODO JV (pt): Optionally, the rport can be set here (somehow). For example by just handing it as a parameter to
  //  handle, or just handing the entire local connection id to the handle function (and storing an array of rports
  //  which can be indexed by the local connection id).

  handle( se );
}

} // namespace

#endif
