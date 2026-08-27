# -*- coding: utf-8 -*-
#
# test_axonal_delay_corrected_weights.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

"""
Check that retrospective spike correction leaves the *corrected* weight right.

``test_stdp_pl_synapse_hom`` compares the weights NEST hands to the weight recorder.
Those are the weights the synapse transmits, which are pre-correction by design, so that
test cannot see whether a correction landed on the right value -- it has to permit any
transmitted weight that a pending correction would change.

This test looks at the weight the synapse is left holding once every correction has been
applied, and compares it against ``send()`` evaluated with the whole postsynaptic history
available, which is what retrospective correction is meant to reproduce. Every spike time
is prescribed, so that reference is exact and the comparison is not statistical.

The parametrisation matters more than its size. A correction entry only has to be handed
on to the next spike on the same connection while two of them are in flight together, i.e.
while the presynaptic inter-spike interval is shorter than ``axonal - dendritic``; the
short-interval cases below are the ones that exercise that hand-off.
"""

from math import exp

import nest
import numpy as np
import pytest

RESOLUTION = 0.1
SIMTIME = 120.0
DRIVE_DELAY = 1.0
DRIVE_WEIGHT = 3000.0
N_PRE = 8
N_POST = 5
INITIAL_WEIGHT = 0.5
TAU_MINUS = 33.7
SYNAPSE_PARAMS = {"lambda": 0.1, "alpha": 0.1, "mu": 0.4, "tau_plus": 20.0}
EPS = 1e-6


def _facilitate(weight, kplus):
    return weight + SYNAPSE_PARAMS["lambda"] * pow(weight, SYNAPSE_PARAMS["mu"]) * kplus


def _depress(weight, kminus):
    new_weight = weight - SYNAPSE_PARAMS["alpha"] * SYNAPSE_PARAMS["lambda"] * weight * kminus
    return new_weight if new_weight > 0.0 else 0.0


def _causal_weight(pre_spikes, post_spikes, axonal_delay, dendritic_delay):
    """The weight ``send()`` arrives at when no postsynaptic spike is ever missed.

    Mirrors stdp_pl_synapse_hom_ax_delay::send(): facilitation walks the history over
    ``(t_last_pre + d_axon - d_dend, t_pre + d_axon - d_dend]``, closed at the upper end
    because ``get_history()`` is, and depression reads the postsynaptic trace strictly
    before ``t_pre + d_axon - d_dend`` because ``get_K_value()`` is strict. ``t_lastspike_``
    starts at ``-d_axon``, standing for a spike that reached the synapse at t=0.
    """
    weight = INITIAL_WEIGHT
    kplus = 0.0
    t_last_pre = -axonal_delay

    for t_pre in pre_spikes:
        t_arrival = t_pre + axonal_delay - dendritic_delay
        for t_post in post_spikes:
            if t_last_pre + axonal_delay - dendritic_delay + EPS < t_post <= t_arrival + EPS:
                minus_dt = t_last_pre + axonal_delay - (t_post + dendritic_delay)
                weight = _facilitate(weight, kplus * exp(minus_dt / SYNAPSE_PARAMS["tau_plus"]))
        kminus = sum(exp((t_post - t_arrival) / TAU_MINUS) for t_post in post_spikes if t_arrival - t_post > EPS)
        weight = _depress(weight, kminus)
        kplus = kplus * exp((t_last_pre - t_pre) / SYNAPSE_PARAMS["tau_plus"]) + 1.0
        t_last_pre = t_pre

    return weight


def _simulate(axonal_delay, dendritic_delay, pre_interval, num_threads):
    """Drive prescribed spike trains through an all-to-all ax-delay projection."""
    nest.ResetKernel()
    nest.set(resolution=RESOLUTION, local_num_threads=num_threads)

    stop = SIMTIME - 2.0 * (SIMTIME / 6.0)
    pre_drive = [np.round(np.arange(1.0 + 0.1 * i, stop, pre_interval + 0.1 * (i % 5)), 1) for i in range(N_PRE)]
    post_drive = [np.round(np.arange(2.0 + 0.1 * j, stop, 1.3 + 0.1 * (j % 7)), 1) for j in range(N_POST)]

    # parrot_neuron repeats its input exactly, so the presynaptic spike times are prescribed;
    # the postsynaptic times are whatever the neurons produce and are read back below.
    pre = nest.Create("parrot_neuron", N_PRE)
    post = nest.Create("iaf_psc_alpha", N_POST, params={"tau_minus": TAU_MINUS, "t_ref": 0.0})
    pre_gen = nest.Create("spike_generator", N_PRE, params=[{"spike_times": t} for t in pre_drive])
    post_gen = nest.Create("spike_generator", N_POST, params=[{"spike_times": t} for t in post_drive])

    nest.Connect(pre_gen, pre, "one_to_one", syn_spec={"synapse_model": "static_synapse", "delay": DRIVE_DELAY})
    nest.Connect(
        post_gen,
        post,
        "one_to_one",
        syn_spec={"synapse_model": "static_synapse", "weight": DRIVE_WEIGHT, "delay": DRIVE_DELAY},
    )

    recorder = nest.Create("spike_recorder")
    nest.Connect(pre + post, recorder, syn_spec={"synapse_model": "static_synapse"})

    nest.SetDefaults("stdp_pl_synapse_hom_ax_delay", SYNAPSE_PARAMS)
    nest.Connect(
        pre,
        post,
        "all_to_all",
        syn_spec={
            "synapse_model": "stdp_pl_synapse_hom_ax_delay",
            "weight": INITIAL_WEIGHT,
            "axonal_delay": axonal_delay,
            "dendritic_delay": dendritic_delay,
        },
    )

    nest.Simulate(SIMTIME)

    senders = recorder.events["senders"]
    times = recorder.events["times"]
    trains = {node: np.sort(times[senders == node]) for node in pre.tolist() + post.tolist()}

    return nest.GetConnections(synapse_model="stdp_pl_synapse_hom_ax_delay"), trains


@pytest.mark.parametrize(
    "axonal_delay, dendritic_delay, pre_interval",
    [
        (2.0, 0.5, 0.7),  # several spikes of one connection critical at once
        (2.0, 0.0, 0.7),  # purely axonal delay, same regime
        (3.0, 0.1, 0.7),  # axonal delay spanning many min-delay intervals
        (2.0, 0.5, 3.0),  # interval longer than axonal - dendritic: no hand-off needed
        (1.0, 0.5, 0.7),  # correction window shorter than one communication interval
        (1.5, 1.5, 0.7),  # axonal == dendritic: no correction can ever be required
    ],
)
# The 4-thread case needs a build with multithreading support; NEST raises on
# local_num_threads > 1 without it.
@pytest.mark.parametrize("num_threads", [1, pytest.param(4, marks=pytest.mark.skipif_missing_threads)])
def test_corrected_weights_match_causal_reference(axonal_delay, dendritic_delay, pre_interval, num_threads):
    """Every synapse must end on the weight an implementation that never misses a spike computes."""
    connections, trains = _simulate(axonal_delay, dendritic_delay, pre_interval, num_threads)

    sources = connections.source
    targets = connections.target
    weights = connections.weight
    if not isinstance(sources, (list, tuple)):
        sources, targets, weights = [sources], [targets], [weights]

    assert len(sources) > 0, "no local ax-delay connections to check"

    reference = {}
    for source, target, weight in zip(sources, targets, weights):
        if (source, target) not in reference:
            reference[(source, target)] = _causal_weight(trains[source], trains[target], axonal_delay, dendritic_delay)
        expected = reference[(source, target)]
        assert weight == pytest.approx(expected, rel=1e-9), (
            f"synapse {source}->{target} ended on {weight!r}, "
            f"but a causally correct implementation gives {expected!r}"
        )


def test_corrections_actually_fire():
    """Guard the guard: the parametrisation above has to enter the correction path at all.

    At ``axonal_delay <= dendritic_delay`` no correction is ever required, so a test that
    only ever ran such a case would pass without exercising anything.
    """
    _simulate(axonal_delay=2.0, dendritic_delay=0.5, pre_interval=0.7, num_threads=1)
    assert nest.kernel_status["num_corrections"] > 0

    _simulate(axonal_delay=1.5, dendritic_delay=1.5, pre_interval=0.7, num_threads=1)
    assert nest.kernel_status["num_corrections"] == 0
