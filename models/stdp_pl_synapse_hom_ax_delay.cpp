/*
 *  stdp_pl_synapse_hom_ax_delay.cpp
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

#include "stdp_pl_synapse_hom_ax_delay.h"

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"


namespace nest
{

//
// Implementation of class STDPPLHomAxDelayCommonProperties.
//

STDPPLHomAxDelayCommonProperties::STDPPLHomAxDelayCommonProperties()
  : CommonSynapseProperties()
  , tau_plus_( 20.0 )
  , tau_plus_inv_( 1. / tau_plus_ )
  , lambda_( 0.1 )
  , alpha_( 1.0 )
  , mu_( 0.4 )
{
}

void
STDPPLHomAxDelayCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );

  def< double >( d, names::tau_plus, tau_plus_ );
  def< double >( d, names::lambda, lambda_ );
  def< double >( d, names::alpha, alpha_ );
  def< double >( d, names::mu, mu_ );
}

void
STDPPLHomAxDelayCommonProperties::set_status( const DictionaryDatum& d, const ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );

  updateValue< double >( d, names::tau_plus, tau_plus_ );
  if ( tau_plus_ > 0. )
  {
    tau_plus_inv_ = 1. / tau_plus_;
  }
  else
  {
    throw BadProperty( "tau_plus > 0. required." );
  }
  updateValue< double >( d, names::lambda, lambda_ );
  updateValue< double >( d, names::alpha, alpha_ );
  updateValue< double >( d, names::mu, mu_ );
}


stdp_pl_synapse_hom_ax_delay::stdp_pl_synapse_hom_ax_delay()
  : ConnectionBase()
  , weight_( 1.0 )
  , Kplus_( 0.0 )
  , t_lastspike_( 0.0 )
{
}

void
stdp_pl_synapse_hom_ax_delay::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );

  // own properties, different for individual synapse
  def< double >( d, names::Kplus, Kplus_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

void
stdp_pl_synapse_hom_ax_delay::set_status( const DictionaryDatum& d, const ConnectorModel& cm )
{
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );

  updateValue< double >( d, names::Kplus, Kplus_ );
}

} // of namespace nest
