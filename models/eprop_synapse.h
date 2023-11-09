/*
 *  eprop_synapse.h
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

#ifndef EPROP_SYNAPSE_H
#define EPROP_SYNAPSE_H

// nestkernel
#include "connection.h"
#include "eprop_archiving_node.h"

namespace nest
{

/* BeginUserDocs: synapse, e-prop plasticity

Short description
+++++++++++++++++

Synapse type for e-prop plasticity

Description
+++++++++++

``eprop_synapse`` is a connector to create e-prop synapses between postsynaptic
neurons :math:`j` and presynaptic neurons and :math:`i` as defined in [1]_.  The
change of the synaptic weight :math:`\Delta W_{ji}` depends on the presynaptic
spikes :math:`z_i^{t-1}`, the pseudo-derivative of the postsynaptic membrane
voltage :math:`\psi_j^t` (which together form the eligibility trace
:math:`e_{ji}`), and the learning signal :math:`L_j^t` emitted by the readout
neurons. The eligibility trace and the presynaptic spike trains are
low-pass-filtered with kernels :math:`\mathcal{F}_\kappa` with
:math:`\kappa=\exp\left(\frac{-\delta t}{\tau_\text{out}}\right)` and
:math:`\mathcal{F}_\alpha` with
:math:`\kappa=\exp\left(\frac{-\delta t}{\tau_\text{m, out}}\right)`.
The overall weight update is scaled by the learning rate :math:`\eta`. The
general formula computing weight updates for eprop synapses projecting onto
recurrent neurons are thus given by:

.. math::
  \Delta W_{ji}^\text{rec} &= -\eta \sum_t L_j^t \bar{e}_{ji}^t \\
   &= -\eta \sum_t L_j^t \mathcal{F}_\kappa \left( \psi^t_j \bar{z}_i^{t-1}\right) \\
   &= -\eta \sum_t L_j^t \sum_{t'\leq t} \kappa^{t-t'} \psi^t_j \mathcal{F}_\alpha\left( z_i^{t-1}\right)\,.

If the postsynaptic neuron is an adaptive neuron, next to the membrane voltage, a
second hidden state variable, the threshold adaptation, is present, which
changes the eligibility trace:

.. math::
  e_{ji}^t &= \psi_j^t \left(\bar{z}^{t-1} - \beta \epsilon_{ji,a}^{t-1}\right)\,, \\
  \epsilon^{t-1}_{ji,\text{a}} &= \psi_j^{t-1}\bar{z}_i^{t-2} + \left( \rho - \psi_j^{t-1} \beta \right)
  \epsilon^{t-2}_{ji,a}\,, \\ \rho &= \exp\left(-\frac{\delta t}{\tau_\text{a}}\right)\,.

Furthermore, a firing rate regularization mechanism keeps the average firing
rate :math:`f^\text{av}_j` of the postsynaptic neuron close to a target firing rate
:math:`f^\text{target}`.

.. math::
  \Delta W_{ji}^\text{reg} = \eta c_\text{reg}
  \sum_t \frac{1}{Tn_\text{trial}} \left( f^\text{target}-f^\text{av}_j\right)e_{ji}^t,

whereby :math:`c_\text{reg}` scales the overall regularization and the average
is taken over the time that passed since the previous update, that is, the number of
trials :math:`n_\text{trials}` times the duration of an update interval :math:`T`.

The overall recurrent weight update is given by adding :math:`\Delta W_{ji}^\text{rec}`
and :math:`\Delta W_{ji}^\text{reg}`.

Since readout neurons :math:`k` are leaky integrators without a spiking
mechanism, the formula for computing weight updates for synapses projecting onto
readout neurons lacks the pseudo derivative and a firing regularization term.

.. math::
  \Delta W_{kj}^\text{out} = -\eta \sum_t L_j^t \mathcal{F}_\kappa \left(z_j^t\right).

The weights can also be optimized with the Adam scheme [2]_:

.. math::
  m_0 &= 0, v_0 = 0, t = 1 \\
  m_t &= \beta_1 \cdot m_{t-1} + \left(1-\beta_1\right) \cdot g_t \\
  v_t &= \beta_2 \cdot v_{t-1} + \left(1-\beta_2\right) \cdot g_t^2 \\
  \hat{m}_t &= \frac{m_t}{1-\beta_1^t} \\
  \hat{v}_t &= \frac{v_t}{1-\beta_2^t} \\
  \Delta W &= - \eta\frac{\hat{m_t}}{\sqrt{\hat{v}_t} + \epsilon}

E-prop synapses require archiving of continuous quantities. Therefore e-prop
synapses can only be connected to neuron models that are capable of doing this
archiving. So far, compatible models are ``eprop_iaf_psc_delta``,
``eprop_iaf_psc_delta_adapt``, and ``eprop_readout``.

For more information on e-prop plasticity, see the documentation on the other e-prop models:

    * :doc:`eprop_iaf_psc_delta<../models/eprop_iaf_psc_delta/>`
    * :doc:`eprop_readout<../models/eprop_readout/>`
    * :doc:`eprop_synapse<../models/eprop_synapse/>`
    * :doc:`eprop_learning_signal_connection<../models/eprop_learning_signal_connection/>`

Details on the event-based NEST implementation of e-prop can be found in [3]_.

.. warning::

   This synaptic plasticity rule does not take
   :ref:`precise spike timing <sim_precise_spike_times>` into
   account. When calculating the weight update, the precise spike time part
   of the timestamp is ignored.

Parameters
++++++++++

The following parameters can be set in the status dictionary.

===============  ========  ================  ================ ====================================================
**Common synapse properties**
------------------------------------------------------------------------------------------------------------------
Parameter        Unit      Math equivalent   Default          Description
===============  ========  ================  ================ ====================================================
adam_beta1                 :math:`\beta_1`   0.9              Beta1 parameter of Adam optimizer
adam_beta2                 :math:`\beta_2`   0.999            Beta2 parameter of Adam optimizer
adam_epsilon               :math:`\epsilon`  1e-8             Epsilon parameter of Adam optimizer
batch_size                                   1                Size of batch
optimizer                                    gradient_descent If adam, use Adam optimizer, if gd, gradient descent
===============  ========  ================  ================ ====================================================

=============  ====  =========================  =======  ===============================================================
**Individual synapse properties**
------------------------------------------------------------------------------------------------------------------------
Parameter      Unit  Math equivalent            Default  Description
=============  ====  =========================  =======  ===============================================================
adam_m               :math:`m`                      0.0  Initial value of first moment estimate m of Adam optimizer
adam_v               :math:`v`                      0.0  Initial value of second moment raw estimate v of Adam optimizer
delay          ms    :math:`d_{ji}`                 1.0  Dendritic delay
eta                  :math:`\eta`                  1e-4  Learning rate
tau_m_readout  ms    :math:`\tau_\text{m,out}`      0.0  Time constant for low-pass filtering of eligibility trace
weight         pA    :math:`W_{ji}`                 1.0  Synaptic weight
Wmax           pA    :math:`W_{ji}^\text{max}`    100.0  Maximal value for synaptic weight
Wmin           pA    :math:`W_{ji}^\text{min}`      0.0  Minimal value for synaptic weight
=============  ====  =========================  =======  ===============================================================

Recordables
+++++++++++

The following variables can be recorded.

  - synaptic weight ``weight``

Usage
+++++

This model can only be used in combination with the other e-prop models,
whereby the network architecture requires specific wiring, input, and output.
The usage is demonstrated in a
:doc:`supervised regression task <../auto_examples/eprop_plasticity/eprop_supervised_regression/>`
and a :doc:`supervised classification task <../auto_examples/eprop_plasticity/eprop_supervised_classification>`,
reproducing the original proof-of-concept tasks in [1]_.

Transmits
+++++++++

SpikeEvent, DSSpikeEvent

References
++++++++++

.. [1] Bellec G, Scherr F, Subramoney F, Hajek E, Salaj D, Legenstein R,
       Maass W (2020). A solution to the learning dilemma for recurrent
       networks of spiking neurons. Nature Communications, 11:3625.
       https://doi.org/10.1038/s41467-020-17236-y

.. [2] Kingma DP, Ba JL (2015). Adam: A method for stochastic optimization.
       Proceedings of International Conference on Learning Representations (ICLR).
       https://doi.org/10.48550/arXiv.1412.6980

.. [3] Korcsak-Gorzo A, Stapmanns J, Espinoza Valverde JA, Dahmen D,
       van Albada SJ, Bolten M, Diesmann M. Event-based implementation of
       eligibility propagation (in preparation)

See also
++++++++

EndUserDocs */

class EpropCommonProperties : public CommonSynapseProperties
{
public:
  EpropCommonProperties();

  void get_status( DictionaryDatum& d ) const;
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  double adam_beta1_;
  double adam_beta2_;
  double adam_epsilon_;
  long batch_size_;
  std::string optimizer_;
  bool average_gradient_;
};

void register_eprop_synapse( const std::string& name );

template < typename targetidentifierT >
class eprop_synapse : public Connection< targetidentifierT >
{

public:
  typedef EpropCommonProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  static constexpr ConnectionModelProperties properties = ConnectionModelProperties::HAS_DELAY
    | ConnectionModelProperties::IS_PRIMARY | ConnectionModelProperties::REQUIRES_EPROP_ARCHIVING
    | ConnectionModelProperties::SUPPORTS_HPC | ConnectionModelProperties::SUPPORTS_LBL;

  eprop_synapse();

  eprop_synapse( const eprop_synapse& ) = default;
  eprop_synapse& operator=( const eprop_synapse& ) = default;

  using ConnectionBase::get_delay;
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void send( Event& e, size_t thread, const EpropCommonProperties& cp );

  void optimize_via_gradient_descent( const long current_optimization_step, const EpropCommonProperties& cp );
  void optimize_via_adam( const long current_optimization_step, const EpropCommonProperties& cp );

  void ( eprop_synapse::*optimize )( const long current_optimization_step, const EpropCommonProperties& cp );

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    using ConnTestDummyNodeBase::handles_test_event;

    size_t
    handles_test_event( SpikeEvent&, size_t )
    {
      return invalid_port;
    }

    size_t
    handles_test_event( DSSpikeEvent&, size_t )
    {
      return invalid_port;
    }
  };

  void
  check_connection( Node& s, Node& t, size_t receptor_type, const CommonPropertiesType& cp )
  {
    // When we get here, delay has been set so we can check it.
    if ( get_delay_steps() != 1 )
    {
      throw KernelException( "eprop synapses currently require a delay of one simulation step" );
    }

    const bool is_source_recurrent_neuron =
      s.get_name() == "eprop_iaf_psc_delta" or s.get_name() == "eprop_iaf_psc_delta_adapt";
    const bool is_target_recurrent_neuron =
      t.get_name() == "eprop_iaf_psc_delta" or t.get_name() == "eprop_iaf_psc_delta_adapt";

    is_recurrent_to_recurrent_conn_ = is_source_recurrent_neuron and is_target_recurrent_neuron;

    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    EpropArchivingNode& t_arch = dynamic_cast< EpropArchivingNode& >( t );
    t_arch.init_update_history();

    const double dt = Time::get_resolution().get_ms();
    kappa_ = exp( -dt / tau_m_readout_ );

    const long update_interval = kernel().simulation_manager.get_eprop_update_interval().get_steps();
    t_next_update_ = update_interval;

    if ( cp.optimizer_ == "adam" )
    {
      optimize = &eprop_synapse::optimize_via_adam;
    }
    else if ( cp.optimizer_ == "gradient_descent" )
    {
      optimize = &eprop_synapse::optimize_via_gradient_descent;
    }
  }

  void
  set_weight( const double w )
  {
    weight_ = w;
  }

protected:
  double weight_;
  double eta_;
  double Wmin_;
  double Wmax_;
  long optimization_step_;
  long t_previous_spike_;
  long t_previous_update_;
  long t_next_update_;
  long t_previous_trigger_spike_;
  double tau_m_readout_; // time constant for low pass filtering of eligibility trace
  double kappa_;         // exp( -dt / tau_m_readout_ )
  double adam_m_;        // auxiliary variable for Adam optimizer
  double adam_v_;        // auxiliary variable for Adam optimizer
  double sum_grads_;     // sum of the gradients in one batch
  bool is_recurrent_to_recurrent_conn_;

  std::vector< long > presyn_isis_;
};

template < typename targetidentifierT >
constexpr ConnectionModelProperties eprop_synapse< targetidentifierT >::properties;

template < typename targetidentifierT >
inline void
eprop_synapse< targetidentifierT >::send( Event& e, size_t thread, const EpropCommonProperties& cp )
{
  EpropArchivingNode* target = dynamic_cast< EpropArchivingNode* >( get_target( thread ) );
  assert( target );

  const long t_spike = e.get_stamp().get_steps();
  const long update_interval = kernel().simulation_manager.get_eprop_update_interval().get_steps();
  const long shift = target->get_shift();

  const long interval_step = ( t_spike - shift ) % update_interval;

  if ( is_recurrent_to_recurrent_conn_ and interval_step == 0 )
  {
    return;
  }

  if ( t_previous_trigger_spike_ == 0 )
  {
    t_previous_trigger_spike_ = t_spike;
  }

  if ( t_previous_spike_ > 0 )
  {
    const long t = t_spike >= t_next_update_ + shift ? t_next_update_ + shift : t_spike;
    presyn_isis_.push_back( t - t_previous_spike_ );
  }

  if ( t_spike >= t_next_update_ + shift )
  {
    const long idx_current_update = ( t_spike - shift ) / update_interval;
    const long t_current_update = idx_current_update * update_interval;
    const long current_optimization_step = 1 + idx_current_update / cp.batch_size_;

    target->write_update_to_history( t_previous_update_, t_current_update );

    sum_grads_ += target->gradient_change(
      presyn_isis_, t_previous_update_, t_previous_trigger_spike_, kappa_, cp.average_gradient_ );

    if ( optimization_step_ < current_optimization_step )
    {
      sum_grads_ /= cp.batch_size_; // mean over batches

      ( this->*optimize )( current_optimization_step, cp );

      weight_ = std::max( Wmin_, std::min( weight_, Wmax_ ) );
    }

    t_previous_update_ = t_current_update;
    t_next_update_ = t_current_update + update_interval;

    t_previous_trigger_spike_ = t_spike;
  }

  t_previous_spike_ = t_spike;

  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
}

template < typename targetidentifierT >
inline void
eprop_synapse< targetidentifierT >::optimize_via_gradient_descent( const long current_optimization_step,
  const EpropCommonProperties& cp )
{
  weight_ -= eta_ * sum_grads_;
  optimization_step_ = current_optimization_step;
  sum_grads_ = 0.0;
}

template < typename targetidentifierT >
inline void
eprop_synapse< targetidentifierT >::optimize_via_adam( const long current_optimization_step,
  const EpropCommonProperties& cp )
{
  for ( ; optimization_step_ < current_optimization_step; ++optimization_step_ )
  {
    const double adam_beta1_factor = 1.0 - std::pow( cp.adam_beta1_, optimization_step_ );
    const double adam_beta2_factor = 1.0 - std::pow( cp.adam_beta2_, optimization_step_ );

    const double alpha_t = eta_ * std::sqrt( adam_beta2_factor ) / adam_beta1_factor;

    adam_m_ = cp.adam_beta1_ * adam_m_ + ( 1.0 - cp.adam_beta1_ ) * sum_grads_;
    adam_v_ = cp.adam_beta2_ * adam_v_ + ( 1.0 - cp.adam_beta2_ ) * sum_grads_ * sum_grads_;

    weight_ -= alpha_t * adam_m_ / ( std::sqrt( adam_v_ ) + cp.adam_epsilon_ );

    sum_grads_ = 0.0; // reset for following iterations
    // since more than 1 cycle through loop indicates past learning periods with vanishing gradients
  }
}

template < typename targetidentifierT >
eprop_synapse< targetidentifierT >::eprop_synapse()
  : ConnectionBase()
  , weight_( 1.0 )
  , eta_( 0.0001 )
  , Wmin_( 0.0 )
  , Wmax_( 100.0 )
  , optimization_step_( 1 )
  , t_previous_spike_( 0 )
  , t_previous_update_( 0 )
  , t_next_update_( 1000 )
  , t_previous_trigger_spike_( 0 )
  , tau_m_readout_( 10.0 )
  , kappa_( 0.0 )
  , adam_m_( 0.0 )
  , adam_v_( 0.0 )
  , sum_grads_( 0.0 )
{
}

template < typename targetidentifierT >
void
eprop_synapse< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::eta, eta_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< double >( d, names::tau_m_readout, tau_m_readout_ );
  def< long >( d, names::size_of, sizeof( *this ) );
  def< double >( d, names::adam_m, adam_m_ );
  def< double >( d, names::adam_v, adam_v_ );
}

template < typename targetidentifierT >
void
eprop_synapse< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::eta, eta_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::tau_m_readout, tau_m_readout_ );
  updateValue< double >( d, names::adam_m, adam_m_ );
  updateValue< double >( d, names::adam_v, adam_v_ );

  if ( weight_ < Wmin_ or weight_ > Wmax_ )
  {
    throw BadProperty( "Wmax >= weight >= Wmin must be satisfied." );
  }

  if ( tau_m_readout_ <= 0 )
  {
    throw BadProperty( "Membrane time constant of readout neuron constant must be > 0." );
  }
}

} // namespace nest

#endif // EPROP_SYNAPSE_H
