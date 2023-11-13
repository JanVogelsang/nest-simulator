# -*- coding: utf-8 -*-
#
# test_metavision.py
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
This set of tests checks the integration of the Metavision SDK and the Metavision backend.
"""

import nest
import numpy.testing as nptest
import pytest


@pytest.fixture
def prepare_kernel():
    nest.ResetKernel()
    nest.local_num_threads = 1
    nest.resolution = 0.1


@pytest.mark.parametrize("h", [0.1, 0.2, 0.5, 1.0])
def test_spike_generator_precise_times_different_resolution(h, expected_spike_times):
    """
    Test the precise times of spikes for different resolutions.
    """
    nest.ResetKernel()
    nest.resolution = h

    sg = nest.Create("spike_generator", {})
    sr = nest.Create("spike_recorder")
    nest.Connect(sg, sr, syn_spec={"delay": 1.0, "weight": 1.0})

    nest.Simulate(7.0)

    actual_spike_times = sr.events["times"]
    nptest.assert_almost_equal(actual_spike_times, [])
