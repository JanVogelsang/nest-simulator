//
// Created by vogelsang on 25.02.23.
//

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
  using Connection::get_status;
  using Connection::set_status;

  double axonal_delay_; //!< Axonal delay in ms
public:
  AxonalDelayConnection()
    : axonal_delay_( 0 )
  {
  }

  void get_status( DictionaryDatum& d ) const;
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Get the proportion of the transmission delay attributed to the dendrite.
   */
  double
  get_dendritic_delay() const
  {
    return get_delay() - get_axonal_delay();
  }

  /**
   * Get the proportion of the transmission delay attributed to the axon.
   */
  double get_axonal_delay() const;

  bool supports_axonal_delay() const;
};

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

inline void
AxonalDelayConnection::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  Connection::set_status( d, cm );

  double axonal_delay;
  if ( updateValue< double >( d, names::axonal_delay, axonal_delay ) )
  {
    if ( axonal_delay < 0.0 ) // consistency with overall delay is checked in check_connection()
    {
      throw BadProperty( "Axonal delay should not be negative." );
    }
    axonal_delay_ = axonal_delay;
  }
}

}
#endif // NEST_AXONAL_DELAY_CONNECTION_H
