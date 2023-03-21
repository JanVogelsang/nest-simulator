/*
 *  axonal_delay_connection.h
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

#ifndef NEST_AXONAL_DELAY_CONNECTION_H
#define NEST_AXONAL_DELAY_CONNECTION_H


// Includes from nestkernel:
#include "connection.h"
#include "exceptions.h"
#include "nest_names.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{

/**
 * TODO JV
 */
class AxonalDelayConnection : public Connection
{
  double axonal_delay_; //!< Axonal delay in ms
public:
  using Connection::get_status;
  using Connection::set_status;

  AxonalDelayConnection()
    : axonal_delay_( 0 )
  {
  }

  void get_status( DictionaryDatum& d ) const;
  //  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Set the proportion of the transmission delay attributed to the axon.
   */
  void set_axonal_delay( const double axonal_delay );

  /**
   * Get the proportion of the transmission delay attributed to the axon.
   */
  double get_axonal_delay() const;

  bool supports_axonal_delay() const;
};

inline void
AxonalDelayConnection::set_axonal_delay( const double axonal_delay )
{
  if ( axonal_delay < 0.0 ) // consistency with overall delay is checked in check_connection()
  {
    throw BadProperty( "Axonal delay should not be negative." );
  }
  axonal_delay_ = axonal_delay;
}

inline double
AxonalDelayConnection::get_axonal_delay() const
{
  return axonal_delay_;
}

// TODO JV (pt): This does probably not belong here but in the ConnectorModel instead
inline bool
AxonalDelayConnection::supports_axonal_delay() const
{
  return true;
}

inline void
AxonalDelayConnection::get_status( DictionaryDatum& d ) const
{
  Connection::get_status( d );

  def< double >( d, names::axonal_delay, axonal_delay_ );
}

// inline void
// AxonalDelayConnection::set_status( const DictionaryDatum& d, ConnectorModel& cm )
//{
//   Connection::set_status( d, cm );
// }

}
#endif // NEST_AXONAL_DELAY_CONNECTION_H
