/*
 *  device.h
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

#ifndef DEVICE_H
#define DEVICE_H


// Includes from nestkernel:
#include "nest_time.h"
#include "nest_types.h"
#include "node.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{

/**
 * @defgroup Devices
 * This group comprises stimulation and recording devices.
 */

/**
 * Class implementing common interface and properties common for all devices.
 *
 * This class provides a common interface for all derived device classes.
 * Each class derived from Node and implementing a device, should have a
 * member derived from class Device. This member will contribute the
 * implementation of device specific properties.
 *
 * This class manages the properties common to all devices, namely
 * origin, start and stop of the time window during which the device
 * is active and the optional device label. The precise semantics of
 * when the device is active depend on the type of device and are
 * defined in subclasses.
 *
 * @ingroup Devices
 *
 * @author HEP 2002-07-22, 2008-03-21, 2008-06-20
 */
class Device : public NodeBase
{
public:
  Device();
  Device( const Device& n );
  ~Device() override
  {
  }

  /** Reset dynamic state to that of model. */
  virtual void
  init_state()
  {
  }

  /** Reset buffers. */
  virtual void
  init_buffers()
  {
  }

  /** Set internal variables before calls to SimulationManager::run() */
  void pre_run_hook() override;

  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum& ) override;

  Name get_element_type() const override;

  bool has_proxies() const override;

  bool is_proxy() const override;

  /**
   *  Returns true if the device is active at the given time stamp.
   *  Semantics are implemented by subclasses.
   */
  virtual bool
  is_active( Time const& T ) const
  {
    return true;
  };

  Time const& get_origin() const;
  Time const& get_start() const;
  Time const& get_stop() const;

  /**
   * Modify Event object parameters during event delivery.
   *
   * Some Nodes want to perform a function on an event for each
   * of their targets. An example is the poisson_generator which
   * needs to draw a random number for each target. The DSSpikeEvent,
   * DirectSendingSpikeEvent, calls sender->event_hook(thread, *this)
   * in its operator() function instead of calling target->handle().
   * The default implementation of Node::event_hook() just calls
   * target->handle(DSSpikeEvent&). Any reimplementation must also
   * execute this call. Otherwise the event will not be delivered.
   * If needed, target->handle(DSSpikeEvent) may be called more than
   * once.
   */
  virtual void event_hook( DSSpikeEvent& );

  virtual void event_hook( DSCurrentEvent& );

  bool one_node_per_process() const override;

protected:
  /**
   * Return lower limit in steps.
   */
  long get_t_min_() const;

  /**
   * Return upper limit in steps.
   */
  long get_t_max_() const;

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    //! Origin of device time axis, relative to network time. Defaults to 0.
    Time origin_;

    //!< Start time, relative to origin. Defaults to 0.
    Time start_;

    //!< Stop time, relative to origin. Defaults to "infinity".
    Time stop_;

    Parameters_(); //!< Sets default parameter values

    //! Copy and recalibrate parameter set
    Parameters_( const Parameters_& );

    Parameters_& operator=( const Parameters_& );

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dictionary

  private:
    //! Update given Time parameter including error checking
    static void update_( const DictionaryDatum&, const Name&, Time& );
  };


  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {

    /**
     * Time step of device activation.
     *
     * t_min_ = origin_ + start_, in steps.
     * @note This is an auxiliary variable that is initialized to -1 in the
     * constructor and set to its proper value by calibrate. It should NOT
     * be returned by get_parameters().
     */
    long t_min_;

    /**
     * Time step of device deactivation.
     *
     * t_max_ = origin_ + stop_, in steps.
     * @note This is an auxiliary variable that is initialized to -1 in the
     * constructor and set to its proper value by calibrate. It should NOT
     * be returned by get_parameters().
     */
    long t_max_;
  };

  // ----------------------------------------------------------------

  Parameters_ P_;
  Variables_ V_;
};

inline void
Device::get_status( DictionaryDatum& d ) const
{
  NodeBase::get_status( d );

  P_.get( d );
}

inline void
Device::set_status( const DictionaryDatum& d )
{
  NodeBase::set_status( d );

  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
}

inline Time const&
Device::get_origin() const
{
  return P_.origin_;
}

inline Time const&
Device::get_start() const
{
  return P_.start_;
}

inline Time const&
Device::get_stop() const
{
  return P_.stop_;
}

inline long
Device::get_t_min_() const
{
  return V_.t_min_;
}

inline long
Device::get_t_max_() const
{
  return V_.t_max_;
}

inline Name
Device::get_element_type() const
{
  return names::device;
}

inline bool
Device::one_node_per_process() const
{
  return false;
}

inline bool
Device::has_proxies() const
{
  return false;
}

inline bool
Device::is_proxy() const
{
  return false;
}

}

#endif /* DEVICE_H */
