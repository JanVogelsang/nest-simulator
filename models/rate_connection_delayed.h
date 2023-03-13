/*
 *  rate_connection_delayed.h
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


#ifndef RATE_CONNECTION_DELAYED_H
#define RATE_CONNECTION_DELAYED_H

#include "connection.h"

namespace nest
{

/* BeginUserDocs: synapse, connection with delay, rate

Short description
+++++++++++++++++

Synapse type for rate connections with delay

Description
+++++++++++

``rate_connection_delayed`` is a connector to create connections with delay
between rate model neurons.

To create instantaneous rate connections please use
the synapse type ``rate_connection_instantaneous``.

See also [1]_.

Transmits
+++++++++

DelayedRateConnectionEvent

References
++++++++++

.. [1] Hahne J, Dahmen D, Schuecker J, Frommer A, Bolten M, Helias M,
       Diesmann M (2017). Integration of continuous-time dynamics in a
       spiking neural network simulator. Frontiers in Neuroinformatics, 11:34.
       DOI: https://doi.org/10.3389/fninf.2017.00034

See also
++++++++

rate_connection_instantaneous, rate_neuron_ipn, rate_neuron_opn

EndUserDocs */

/**
 * Class representing a delayed rate connection. A rate_connection_delayed
 * has the properties weight, delay and receiver port.
 */
class RateConnectionDelayed : public Connection
{

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection ConnectionBase;
  typedef DelayedRateConnectionEvent EventType;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  RateConnectionDelayed()
    : ConnectionBase()
    , weight_( 1.0 )
  {
  }

  void
  check_connection( Node& s,
    Node& t,
    const rport receptor_type,
    const synindex syn_id,
    const delay dendritic_delay,
    const delay axonal_delay,
    const CommonPropertiesType& )
  {
    EventType ge;

    s.sends_secondary_event( ge );
    ge.set_sender( s );
  }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   */
  void
  send( Event& e, const thread, const delay axonal_delay, const delay dendritic_delay, const CommonSynapseProperties&, Node* )
  {
    e.set_weight( weight_ );
    e.set_delay_steps( dendritic_delay );
    e();
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  double weight_; //!< connection weight
};

void
RateConnectionDelayed::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

void
RateConnectionDelayed::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef RATE_CONNECTION_DELAYED_H */
