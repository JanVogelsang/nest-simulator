/*
 *  precise_weighted_spike_generator.cpp
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

#include "precise_weighted_spike_generator.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "nest_impl.h"

// Includes from sli:
#include "arraydatum.h"
#include "dict.h"
#include "dictutils.h"


void
nest::register_precise_weighted_spike_generator( const std::string& name )
{
  register_node_model< precise_weighted_spike_generator >( name );
}

namespace nest
{

/* ----------------------------------------------------------------
 * Default constructor defining default parameters
 * ---------------------------------------------------------------- */

precise_weighted_spike_generator::Parameters_::Parameters_()
  : spike_stamps_()
  , spike_offsets_()
  , spike_weights_()
{
}


/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void
precise_weighted_spike_generator::Parameters_::get( DictionaryDatum& d ) const
{
  const size_t n_spikes = spike_stamps_.size();
  auto* times_ms = new std::vector< double >();
  times_ms->reserve( n_spikes );

  for ( size_t n = 0; n < n_spikes; ++n )
  {
    times_ms->push_back( spike_stamps_[ n ].get_ms() );
    ( *times_ms )[ n ] -= spike_offsets_[ n ];
  }

  ( *d )[ names::spike_times ] = DoubleVectorDatum( times_ms );
  ( *d )[ names::spike_weights ] = DoubleVectorDatum( new std::vector< double >( spike_weights_ ) );
}

void
precise_weighted_spike_generator::Parameters_::assert_valid_spike_time_and_insert_( double t, const Time&, const Time& )
{
  Time t_spike = Time::ms_stamp( t );

  spike_stamps_.push_back( t_spike );

  // t_spike is created with ms_stamp() that aligns the time to the next resolution step, so the offset has to be
  // greater or equal to t by construction. Since subtraction of close-by floating point values is not stable, we have
  // to compare with a delta.
  double offset = t_spike.get_ms() - t;

  // The second part of the test handles subnormal values of offset.
  if ( ( std::fabs( offset ) < std::numeric_limits< double >::epsilon() * std::fabs( t_spike.get_ms() + t ) * 2.0 )
    or ( std::fabs( offset ) < std::numeric_limits< double >::min() ) )
  {
    // if difference is smaller than scaled epsilon it is zero
    offset = 0.0;
  }
  assert( offset >= 0.0 );
  spike_offsets_.push_back( offset );
}

void
precise_weighted_spike_generator::Parameters_::set( const DictionaryDatum& d,
  State_& s,
  const Time& origin,
  const Time& now,
  Node* )
{
  const bool updated_spike_times = d->known( names::spike_times );

  if ( updated_spike_times )
  {
    const std::vector< double > d_times = getValue< std::vector< double > >( d->lookup( names::spike_times ) );
    const size_t n_spikes = d_times.size();
    spike_stamps_.clear();
    spike_stamps_.reserve( n_spikes );
    spike_offsets_.clear();
    spike_offsets_.reserve( n_spikes );

    // Check spike times for ordering and grid compatibility and insert them
    if ( not d_times.empty() )
    {
      // handle first spike time, no predecessor to compare with
      auto prev = d_times.begin();
      assert_valid_spike_time_and_insert_( *prev, origin, now );

      // handle all remaining spike times, compare to predecessor
      for ( auto next = prev + 1; next != d_times.end(); ++next, ++prev )
      {
        if ( *prev > *next )
        {
          throw BadProperty( "Spike times must be sorted in non-descending order." );
        }
        else
        {
          assert_valid_spike_time_and_insert_( *next, origin, now );
        }
      }
    }
  }

  // spike_weights can be the same size as spike_times, or can be of size 0 to
  // only use the spike_times array
  bool updated_spike_weights = d->known( names::spike_weights );
  if ( updated_spike_weights )
  {
    std::vector< double > spike_weights = getValue< std::vector< double > >( d->lookup( names::spike_weights ) );

    if ( spike_weights.empty() )
    {
      spike_weights_.clear();
    }
    else
    {
      if ( spike_weights.size() != spike_stamps_.size() )
      {
        throw BadProperty(
          "spike_weights must have the same number of elements as spike_times,"
          " or 0 elements to clear the property." );
      }

      spike_weights_.swap( spike_weights );
    }
  }

  // Set position to start if something changed
  if ( updated_spike_times or updated_spike_weights or d->known( names::origin ) )
  {
    s.position_ = 0;
  }
}


/* ----------------------------------------------------------------
 * Default constructor defining default state
 * ---------------------------------------------------------------- */

precise_weighted_spike_generator::State_::State_()
  : position_( 0 )
{
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

precise_weighted_spike_generator::precise_weighted_spike_generator()
  : Node()
  , StimulationDevice()
  , P_()
  , S_()
{
}

precise_weighted_spike_generator::precise_weighted_spike_generator( const precise_weighted_spike_generator& n )
  : Node( n )
  , StimulationDevice( n )
  , P_( n.P_ )
  , S_( n.S_ )
{
}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
precise_weighted_spike_generator::init_state_()
{
  StimulationDevice::init_state();
}

void
precise_weighted_spike_generator::init_buffers_()
{
  StimulationDevice::init_buffers();
}

void
precise_weighted_spike_generator::pre_run_hook()
{
  StimulationDevice::pre_run_hook();
}


/* ----------------------------------------------------------------
 * Other functions
 * ---------------------------------------------------------------- */
void
precise_weighted_spike_generator::update( Time const& sliceT0, const long from, const long to )
{
  if ( P_.spike_stamps_.empty() )
  {
    return;
  }

  assert( P_.spike_stamps_.size() == P_.spike_offsets_.size() );
  assert( P_.spike_weights_.empty() or P_.spike_stamps_.size() == P_.spike_weights_.size() );

  const Time tstart = sliceT0 + Time::step( from );
  const Time tstop = sliceT0 + Time::step( to );
  const Time& origin = StimulationDevice::get_origin();

  // We fire all spikes with time stamps up to including sliceT0 + to
  while ( S_.position_ < P_.spike_stamps_.size() )
  {
    const Time tnext_stamp = origin + P_.spike_stamps_[ S_.position_ ];

    // this might happen due to wrong usage of the generator
    if ( tnext_stamp <= tstart )
    {
      ++S_.position_;
      continue;
    }
    if ( tnext_stamp > tstop )
    {
      break;
    }

    if ( StimulationDevice::is_active( tnext_stamp ) )
    {
      SpikeEvent se;

      se.set_offset( P_.spike_offsets_[ S_.position_ ] );

      // we need to subtract one from stamp which is added again in send()
      long lag = Time( tnext_stamp - sliceT0 ).get_steps() - 1;

      // all spikes are sent locally, so offset information is always preserved
      kernel().event_delivery_manager.send( *this, se, lag );
    }

    ++S_.position_;
  }
}

void
precise_weighted_spike_generator::set_data_from_stimulation_backend( std::vector< double >& input_spikes )
{
  if ( not input_spikes.empty() )
  {
    Parameters_ ptmp = P_; // temporary copy in case of errors
    const Time& origin = StimulationDevice::get_origin();

    DictionaryDatum d = DictionaryDatum( new Dictionary );
    std::vector< double > times_ms;
    std::vector< double > weights;
    const size_t n_spikes = P_.spike_stamps_.size();
    times_ms.reserve( n_spikes + input_spikes.size() / 2 );
    weights.reserve( n_spikes + input_spikes.size() / 2 );
    for ( size_t n = 0; n < n_spikes; ++n )
    {
      times_ms.push_back( P_.spike_stamps_[ n ].get_ms() - P_.spike_offsets_[ n ] );
      weights.push_back( P_.spike_weights_[ n ] );
    }
    for ( size_t n = 0; n < input_spikes.size() / 2; ++n )
    {
      times_ms.push_back( input_spikes[ n * 2 ] );
      weights.push_back( input_spikes[ n * 2 + 1 ] );
    }

    ( *d )[ names::spike_times ] = DoubleVectorDatum( times_ms );
    ( *d )[ names::spike_weights ] = DoubleVectorDatum( weights );

    assert( times_ms.size() > 0 );
    ptmp.set( d, S_, origin, Time::step( times_ms[ times_ms.size() - 1 ] ), this );

    // if we get here, temporary contains consistent set of properties
    P_ = ptmp;
  }
}

}