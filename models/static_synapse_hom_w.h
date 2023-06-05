/*
 *  static_synapse_hom_w.h
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

#ifndef STATICSYNAPSE_HOM_W_H
#define STATICSYNAPSE_HOM_W_H

// Includes from nestkernel:
#include "common_properties_hom_w.h"
#include "connection.h"

namespace nest
{

/* BeginUserDocs: synapse, static

Short description
+++++++++++++++++

Synapse type for static connections with homogeneous weight

Description
+++++++++++

``static_synapse_hom_w`` does not support any kind of plasticity. It simply
stores the parameters delay, target, and receiver port for each connection
and uses a common weight for all connections.

The common weight for all connections of this model must be set by
``SetDefaults`` on the model. If you create copies of this model using
``CopyModel``, each derived model can have a different weight.

Transmits
+++++++++

SpikeEvent, RateEvent, CurrentEvent, ConductanceEvent,
DataLoggingRequest, DoubleDataEvent

See also
++++++++

static_synapse

EndUserDocs */

class static_synapse_hom_w : public Connection
{

public:
  // this line determines which common properties to use
  typedef CommonPropertiesHomW CommonPropertiesType;
  typedef Connection ConnectionBase;

  using ConnectionBase::get_dendritic_delay_steps;

  void get_status( DictionaryDatum& d ) const;

  void
  check_connection( Node& s,
    Node& t,
    const rport receptor_type,
    const synindex syn_id,
    const delay dendritic_delay,
    const delay axonal_delay,
    const CommonPropertiesType& )
  {
  }

  /**
   * Checks to see if weight is given in syn_spec.
   */
  void
  check_synapse_params( const DictionaryDatum& syn_spec ) const
  {
    if ( syn_spec->known( names::weight ) )
    {
      throw BadProperty(
        "Weight cannot be specified since it needs to be equal "
        "for all connections when static_synapse_hom_w is used." );
    }
  }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param tid Thread ID of the target
   * \param cp Common properties-object of the synapse
   */
  void
  send( Event& e, const thread, const double axonal_delay, const CommonPropertiesHomW& cp, Node* )
  {
    e.set_weight( cp.get_weight() );
    e.set_delay_steps( get_dendritic_delay_steps() + Time::delay_ms_to_steps( axonal_delay ) );
  }

  void
  set_weight( double )
  {
    throw BadProperty(
      "Setting of individual weights is not possible! The common weights can "
      "be changed via "
      "CopyModel()." );
  }
};


void
static_synapse_hom_w::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< long >( d, names::size_of, sizeof( *this ) );
}

} // namespace

#endif /* #ifndef STATICSYNAPSE_HOM_W_H */
