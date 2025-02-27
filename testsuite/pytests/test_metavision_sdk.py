# -*- coding: utf-8 -*-
#
# test_metavision_sdk.py
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
Test of Metavision SDK
"""

import unittest
from pathlib import Path

import nest


@nest.ll_api.check_stack
class MetavisionSDKTestCase(unittest.TestCase):
    """Tests integration of the Metavision SDK"""

    def test_EventsSpikes(self):
        """Spike Events"""
        nest.ResetKernel()
        nest.set_verbosity("M_WARNING")
        nest.resolution = 1.0
        nest.local_num_threads = 2

        nest.SetDefaults(
            "metavision", {"file_paths": [(Path(__file__).parent / "data/metavision_test.hdf5").resolve()]}
        )
        # nest.SetDefaults("metavision", {"serial_numbers": ["ABC123"]})

        cam = nest.Create("precise_weighted_spike_generator", 1280 * 720, params={"stimulus_source": "metavision"})
        sr = nest.Create("spike_recorder")

        nest.Connect(cam, sr, syn_spec={"delay": 1.0})

        nest.Simulate(0.095)

        d = nest.GetStatus(sr, "events")[0]

        self.assertGreater(len(d["times"]), 0)

        print(nest.kernel_status["time_simulate"])
        print(nest.kernel_status["local_spike_counter"])


def suite():
    suite = unittest.makeSuite(MetavisionSDKTestCase, "test")
    return suite


if __name__ == "__main__":
    # runner = unittest.TextTestRunner(verbosity=2)
    # runner.run(suite())
    MetavisionSDKTestCase().test_EventsSpikes()
