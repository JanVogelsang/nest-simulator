/*
 *  connector_model.h
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

#ifndef CONNECTOR_MODEL_H
#define CONNECTOR_MODEL_H

// C++ includes:
#include <cmath>
#include <string>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "connection_type_enum.h"
#include "nest_time.h"
#include "nest_types.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{
class CommonSynapseProperties;
class Event;
class Node;
class SecondaryEvent;
class TimeConverter;

class ConnectorModel
{

public:
  ConnectorModel( const std::string,
    const bool is_primary,
    const bool has_delay,
    const bool requires_symmetric,
    const bool supports_wfr,
    const bool requires_clopath_archiving,
    const bool requires_urbanczik_archiving );
  ConnectorModel( const ConnectorModel&, const std::string );
  virtual ~ConnectorModel()
  {
  }

  /**
   * Adds a connection.
   *
   * @param src Source node
   * @param tgt Target node
   * @param hetconn Connector vector
   * @param syn_id Synapse id
   * @param d Parameter dictionary to configure the synapse
   * @param delay Delay of the connection
   * @param weight Weight of the connection
   *
   * Delay and weight have the default value NAN, a special value, which
   * describes double values that are not a number. If delay or weight is
   * omitted, NAN indicates this and weight/delay are set only if they are
   * valid.
   */
  virtual const std::tuple< index, delay, delay > add_connection( Node& src,
    Node& tgt,
    const synindex syn_id,
    const DictionaryDatum& d,
    const double delay,
    const double axonal_delay,
    const double weight,
    const bool is_primary,
    const ConnectionType connection_type ) = 0;

  virtual ConnectorModel* clone( std::string, synindex syn_id ) const = 0;

  virtual void calibrate( const TimeConverter& tc ) = 0;

  virtual void get_status( DictionaryDatum& ) const = 0;
  virtual void set_status( const DictionaryDatum& ) = 0;

  virtual const CommonSynapseProperties& get_common_properties() const = 0;

  /**
   * Checks to see if illegal parameters are given in syn_spec.
   */
  virtual void check_synapse_params( const DictionaryDatum& ) const = 0;

  virtual SecondaryEvent* get_event() const = 0;

  virtual void set_syn_id( synindex syn_id ) = 0;

  virtual SecondaryEvent* create_event() const = 0;

  std::string
  get_name() const
  {
    return name_;
  }

  bool
  is_primary() const
  {
    return is_primary_;
  }

  bool
  has_delay() const
  {
    return has_delay_;
  }

  bool
  requires_symmetric() const
  {
    return requires_symmetric_;
  }

  bool
  requires_clopath_archiving() const
  {
    return requires_clopath_archiving_;
  }

  bool
  requires_urbanczik_archiving() const
  {
    return requires_urbanczik_archiving_;
  }

  bool
  supports_wfr() const
  {
    return supports_wfr_;
  }

protected:
  //! name of the ConnectorModel
  std::string name_;
  //! indicates whether the default delay must be checked
  bool default_delay_needs_check_;
  //! indicates whether this ConnectorModel belongs to a primary connection
  bool is_primary_;
  //! indicates whether ConnectorModel has a delay
  bool has_delay_;
  //! indicates that ConnectorModel requires symmetric connections
  bool requires_symmetric_;
  //! indicates whether connection can be used during wfr update
  bool supports_wfr_;
  //! indicates that ConnectorModel requires Clopath archiving
  bool requires_clopath_archiving_;
  //! indicates that ConnectorModel requires Urbanczik archiving
  bool requires_urbanczik_archiving_;

}; // ConnectorModel
template < typename ConnectionT >
class GenericConnectorModel : public ConnectorModel
{
private:
  typename ConnectionT::CommonPropertiesType cp_;
  //! used to create secondary events that belong to secondary connections
  typename ConnectionT::EventType* pev_;

  ConnectionT default_connection_;
  delay default_delay_;
  delay default_axonal_delay_;
  rport receptor_type_;

public:
  GenericConnectorModel( const std::string name,
    bool is_primary,
    bool has_delay,
    bool requires_symmetric,
    bool supports_wfr,
    bool requires_clopath_archiving,
    bool requires_urbanczik_archiving )
    : ConnectorModel( name,
      is_primary,
      has_delay,
      requires_symmetric,
      supports_wfr,
      requires_clopath_archiving,
      requires_urbanczik_archiving )
    , default_axonal_delay_( 0.0 )
    , receptor_type_( 0 )
  {
    default_delay_ = Time::delay_ms_to_steps( 1.0 );
  }

  GenericConnectorModel( const GenericConnectorModel& cm, const std::string name )
    : ConnectorModel( cm, name )
    , cp_( cm.cp_ )
    , pev_( cm.pev_ )
    , default_connection_( cm.default_connection_ )
    , default_delay_( cm.default_delay_ )
    , default_axonal_delay_( cm.default_axonal_delay_ )
    , receptor_type_( cm.receptor_type_ )
  {
  }

  const std::tuple< index, delay, delay > add_connection( Node& src,
    Node& tgt,
    const synindex syn_id,
    const DictionaryDatum& d,
    const double delay,
    const double axonal_delay,
    const double weight,
    const bool is_primary,
    const ConnectionType connection_type ) override;

  ConnectorModel* clone( std::string, synindex ) const override;

  void calibrate( const TimeConverter& tc ) override;

  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum& ) override;

  void
  check_synapse_params( const DictionaryDatum& syn_spec ) const override
  {
    default_connection_.check_synapse_params( syn_spec );
  }

  typename ConnectionT::CommonPropertiesType const&
  get_common_properties() const override
  {
    return cp_;
  }

  void set_syn_id( synindex syn_id ) override;

  typename ConnectionT::EventType*
  get_event() const override
  {
    assert( false );
    return 0;
  }

  ConnectionT const&
  get_default_connection() const
  {
    return default_connection_;
  }

  SecondaryEvent*
  create_event() const override
  {
    // Must not be called for a ConnectorModel belonging to a primary
    // connection. Only required for secondary connection types.
    assert( false );
    return nullptr; // make the compiler happy
  }

private:
  void used_default_delay();

}; // GenericConnectorModel

template < typename ConnectionT >
class GenericSecondaryConnectorModel : public GenericConnectorModel< ConnectionT >
{
private:
  //! used to create secondary events that belong to secondary connections
  typename ConnectionT::EventType* pev_;

public:
  GenericSecondaryConnectorModel( const std::string name,
    const bool has_delay,
    const bool requires_symmetric,
    const bool supports_wfr )
    : GenericConnectorModel< ConnectionT >( name,
      /*is _primary=*/false,
      has_delay,
      requires_symmetric,
      supports_wfr,
      /*requires_clopath_archiving=*/false,
      /*requires_urbanczik_archiving=*/false )
    , pev_( 0 )
  {
    pev_ = new typename ConnectionT::EventType();
  }

  GenericSecondaryConnectorModel( const GenericSecondaryConnectorModel& cm, const std::string name )
    : GenericConnectorModel< ConnectionT >( cm, name )
  {
    pev_ = new typename ConnectionT::EventType( *cm.pev_ );
  }
  ConnectorModel* clone( std::string name, synindex syn_id ) const;

  SecondaryEvent*
  create_event() const
  {
    return new typename ConnectionT::EventType();
  }
  ~GenericSecondaryConnectorModel()
  {
    if ( pev_ != 0 )
    {
      delete pev_;
    }
  }

  typename ConnectionT::EventType*
  get_event() const
  {
    return pev_;
  }
};

} // namespace nest

#endif /* #ifndef CONNECTOR_MODEL_H */
