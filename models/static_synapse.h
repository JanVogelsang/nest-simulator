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

  void
  check_connection( Node&, Node&, rport, const synindex, const delay, const CommonPropertiesType& )
  {
  }

  void
  send( Event& e, const thread tid, const double, const CommonSynapseProperties& )
  {
#ifdef TIMER_DETAILED
    if ( tid == 0 )
    {
      kernel().event_delivery_manager.sw_deliver_conn_.stop();
      kernel().event_delivery_manager.sw_static_delivery_.start();
    }
#endif
    e.set_weight( weight_ );
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, const ConnectorModel& cm );

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
static_synapse::set_status( const DictionaryDatum& d, const ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef STATICSYNAPSE_H */
