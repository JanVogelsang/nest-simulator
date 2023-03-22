/*
 *  static_synapse.h
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

#ifndef STATICSYNAPSE_H
#define STATICSYNAPSE_H

// Includes from nestkernel:
#include "connection.h"

namespace nest
{

/* BeginUserDocs: synapse, static

Short description
+++++++++++++++++

Synapse type for static connections

Description
+++++++++++

``static_synapse`` does not support any kind of plasticity. It simply stores
the parameters target, weight, delay and receiver port for each connection.

Transmits
+++++++++

SpikeEvent, RateEvent, CurrentEvent, ConductanceEvent,
DoubleDataEvent, DataLoggingRequest

See also
++++++++

tsodyks_synapse, stdp_synapse

EndUserDocs */

class static_synapse : public Connection
{
  double weight_;

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;

  typedef Connection ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  static_synapse()
    : ConnectionBase()
    , weight_( 1.0 )
  {
  }

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  static_synapse( const static_synapse& rhs ) = default;
  static_synapse& operator=( const static_synapse& rhs ) = default;

  using ConnectionBase::get_dendritic_delay_steps;


  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport ) override
    {
      return invalid_port;
    }
    port
    handles_test_event( RateEvent&, rport ) override
    {
      return invalid_port;
    }
    port
    handles_test_event( DataLoggingRequest&, rport ) override
    {
      return invalid_port;
    }
    port
    handles_test_event( CurrentEvent&, rport ) override
    {
      return invalid_port;
    }
    port
    handles_test_event( ConductanceEvent&, rport ) override
    {
      return invalid_port;
    }
    port
    handles_test_event( DoubleDataEvent&, rport ) override
    {
      return invalid_port;
    }
    port
    handles_test_event( DSSpikeEvent&, rport ) override
    {
      return invalid_port;
    }
    port
    handles_test_event( DSCurrentEvent&, rport ) override
    {
      return invalid_port;
    }
  };

  void
  check_connection( Node&, Node&, const rport, const synindex, const delay, const delay, const CommonPropertiesType& )
  {
  }

  void
  send( Event& e, const thread, const double axonal_delay, const CommonSynapseProperties&, Node* )
  {
    e.set_weight( weight_ );
    e.set_delay_steps( get_dendritic_delay_steps() + Time::delay_ms_to_steps( axonal_delay ) );
    e();
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void
  set_weight( double w )
  {
    weight_ = w;
  }
};

void
static_synapse::get_status( DictionaryDatum& d ) const
{

  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

void
static_synapse::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef STATICSYNAPSE_H */
