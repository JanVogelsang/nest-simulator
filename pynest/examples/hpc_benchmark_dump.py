import os
import sys
import time

import nest
import nest.raster_plot
import numpy as np
import scipy.special as sp

M_INFO = 10
M_ERROR = 30

###############################################################################
# Parameter section
# Define all relevant parameters: changes should be made here


params = {
    "num_threads": 8,  # total number of threads per process
    "scale": 0.01,  # scaling factor of the network size
    # total network size = scale*11250 neurons
    "simtime": 1000.0,  # total simulation time in ms
    "presimtime": 500.0,  # simulation time until reaching equilibrium
    "dt": 0.1,  # simulation step
    "record_spikes": False,  # switch to record spikes of excitatory
    # neurons to file
    "path_name": ".",  # path where all files will have to be written
    "log_file": "log",  # naming scheme for the log files
}


def convert_synapse_weight(tau_m, tau_syn, C_m):
    """
    Computes conversion factor for synapse weight from mV to pA.

    This function is specific to the leaky integrate-and-fire neuron model with alpha-shaped postsynaptic currents.
    """

    # compute time to maximum of V_m after spike input
    # to neuron at rest
    a = tau_m / tau_syn
    b = 1.0 / tau_syn - 1.0 / tau_m
    t_rise = 1.0 / b * (-lambertwm1(-np.exp(-1.0 / a) / a).real - 1.0 / a)

    v_max = (
        np.exp(1.0)
        / (tau_syn * C_m * b)
        * ((np.exp(-t_rise / tau_m) - np.exp(-t_rise / tau_syn)) / b - t_rise * np.exp(-t_rise / tau_syn))
    )
    return 1.0 / v_max


###############################################################################
# For compatibility with earlier benchmarks, we require a rise time of
# ``t_rise = 1.700759 ms`` and we choose ``tau_syn`` to achieve this for given
# ``tau_m``. This requires numerical inversion of the expression for ``t_rise``
# in ``convert_synapse_weight``. We computed this value once and hard-code
# it here.


tau_syn = 0.32582722403722841

brunel_params = {
    "NE": int(9000 * params["scale"]),  # number of excitatory neurons
    "NI": int(2250 * params["scale"]),  # number of inhibitory neurons
    "model_params": {  # Set variables for iaf_psc_alpha
        "E_L": 0.0,  # Resting membrane potential(mV)
        "C_m": 250.0,  # Capacity of the membrane(pF)
        "tau_m": 10.0,  # Membrane time constant(ms)
        "t_ref": 0.5,  # Duration of refractory period(ms)
        "V_th": 20.0,  # Threshold(mV)
        "V_reset": 0.0,  # Reset Potential(mV)
        # time const. postsynaptic excitatory currents(ms)
        "tau_syn_ex": tau_syn,
        # time const. postsynaptic inhibitory currents(ms)
        "tau_syn_in": tau_syn,
        "tau_minus": 30.0,  # time constant for STDP(depression)
        # V can be randomly initialized see below
        "V_m": 5.7,  # mean value of membrane potential
    },
    ####################################################################
    # Note that Kunkel et al. (2014) report different values. The values
    # in the paper were used for the benchmarks on K, the values given
    # here were used for the benchmark on JUQUEEN.
    "randomize_Vm": True,
    "mean_potential": 5.7,
    "sigma_potential": 7.2,
    "delay": 1.5,  # synaptic delay, all connections(ms)
    # synaptic weight
    "JE": 0.14,  # peak of EPSP
    "sigma_w": 3.47,  # standard dev. of E->E synapses(pA)
    "g": -5.0,
    "stdp_params": {
        "delay": 1.5,
        "alpha": 0.0513,
        "lambda": 0.1,  # STDP step size
        "mu": 0.4,  # STDP weight dependence exponent(potentiation)
        "tau_plus": 15.0,  # time constant for potentiation
    },
    "eta": 1.685,  # scaling of external stimulus
    "filestem": params["path_name"],
}

###############################################################################
# Function Section


def build_network(logger):
    """Builds the network including setting of simulation and neuron parameters, creation of neurons and connections.

    Requires an instance of Logger as argument.
    """
    tic = time.time()  # start timer on construction

    # unpack a few variables for convenience
    NE = brunel_params["NE"]
    NI = brunel_params["NI"]
    model_params = brunel_params["model_params"]
    stdp_params = brunel_params["stdp_params"]

    # set global kernel parameters
    nest.local_num_threads = params["num_threads"]
    nest.resolution = params["dt"]
    nest.overwrite_files = True
    nest.rng_seed = 55
    nest.set_verbosity(M_ERROR)

    nest.message(M_INFO, "build_network", "Creating excitatory population.")
    E_neurons = nest.Create("iaf_psc_alpha", NE, params=model_params)

    nest.message(M_INFO, "build_network", "Creating inhibitory population.")
    I_neurons = nest.Create("iaf_psc_alpha", NI, params=model_params)

    if brunel_params["randomize_Vm"]:
        nest.message(M_INFO, "build_network", "Randomizing membrane potentials.")

        random_vm = nest.random.normal(brunel_params["mean_potential"], brunel_params["sigma_potential"])
        nest.GetLocalNodeCollection(E_neurons).V_m = random_vm
        nest.GetLocalNodeCollection(I_neurons).V_m = random_vm

    # number of incoming excitatory connections
    CE = int(1.0 * NE / params["scale"])
    # number of incomining inhibitory connections
    CI = int(1.0 * NI / params["scale"])

    nest.message(M_INFO, "build_network", "Creating excitatory stimulus generator.")

    # Convert synapse weight from mV to pA
    conversion_factor = convert_synapse_weight(model_params["tau_m"], model_params["tau_syn_ex"], model_params["C_m"])
    JE_pA = conversion_factor * brunel_params["JE"]

    nu_thresh = model_params["V_th"] / (
        CE * model_params["tau_m"] / model_params["C_m"] * JE_pA * np.exp(1.0) * tau_syn
    )
    nu_ext = nu_thresh * brunel_params["eta"]

    E_stimulus = nest.Create("poisson_generator", 1, {"rate": nu_ext * CE * 1000.0})

    nest.message(M_INFO, "build_network", "Creating excitatory spike recorder.")

    if params["record_spikes"]:
        recorder_label = os.path.join(brunel_params["filestem"], "alpha_" + str(stdp_params["alpha"]) + "_spikes")
        E_recorder = nest.Create("spike_recorder", params={"record_to": "ascii", "label": recorder_label})

    BuildNodeTime = time.time() - tic

    logger.log(str(BuildNodeTime) + " # build_time_nodes")
    logger.log(str(memory_thisjob()) + " # virt_mem_after_nodes")

    tic = time.time()

    nest.SetDefaults("static_synapse_hpc", {"delay": brunel_params["delay"]})
    nest.CopyModel("static_synapse_hpc", "syn_ex", {"weight": JE_pA})
    nest.CopyModel("static_synapse_hpc", "syn_in", {"weight": brunel_params["g"] * JE_pA})

    stdp_params["weight"] = JE_pA
    nest.SetDefaults("stdp_pl_synapse_hom_hpc", stdp_params)

    nest.message(M_INFO, "build_network", "Connecting stimulus generators.")

    # Connect Poisson generator to neuron
    nest.Connect(E_stimulus, E_neurons, {"rule": "all_to_all"}, {"synapse_model": "syn_ex"})
    nest.Connect(E_stimulus, I_neurons, {"rule": "all_to_all"}, {"synapse_model": "syn_ex"})

    if not os.path.exists("connections_0.dat") and not os.path.exists("connections.sion"):
        nest.message(M_INFO, "build_network", "Connecting excitatory -> excitatory population.")

        nest.Connect(
            E_neurons,
            E_neurons,
            {"rule": "fixed_indegree", "indegree": CE, "allow_autapses": False, "allow_multapses": True},
            {"synapse_model": "stdp_pl_synapse_hom_hpc"},
        )

        nest.message(M_INFO, "build_network", "Connecting inhibitory -> excitatory population.")

        nest.Connect(
            I_neurons,
            E_neurons,
            {"rule": "fixed_indegree", "indegree": CI, "allow_autapses": False, "allow_multapses": True},
            {"synapse_model": "syn_in"},
        )

        nest.message(M_INFO, "build_network", "Connecting excitatory -> inhibitory population.")

        nest.Connect(
            E_neurons,
            I_neurons,
            {"rule": "fixed_indegree", "indegree": CE, "allow_autapses": False, "allow_multapses": True},
            {"synapse_model": "syn_ex"},
        )

        nest.message(M_INFO, "build_network", "Connecting inhibitory -> inhibitory population.")

        nest.Connect(
            I_neurons,
            I_neurons,
            {"rule": "fixed_indegree", "indegree": CI, "allow_autapses": False, "allow_multapses": True},
            {"synapse_model": "syn_in"},
        )
        nest.Prepare()
        exit(0)

    if params["record_spikes"]:
        if params["num_threads"] != 1:
            local_neurons_E = nest.GetLocalNodeCollection(E_neurons)
            local_neurons_I = nest.GetLocalNodeCollection(I_neurons)
            # GetLocalNodeCollection returns a stepped composite NodeCollection, which
            # cannot be sliced. In order to allow slicing it later on, we're creating a
            # new regular NodeCollection from the plain node IDs.
            local_neurons_E = nest.NodeCollection(local_neurons_E.tolist())
            local_neurons_I = nest.NodeCollection(local_neurons_I.tolist())
        else:
            local_neurons_E = E_neurons
            local_neurons_I = I_neurons

        nest.message(M_INFO, "build_network", "Connecting spike recorders.")
        nest.Connect(local_neurons_E, E_recorder, "all_to_all", "static_synapse_hpc")
        nest.Connect(local_neurons_I, E_recorder, "all_to_all", "static_synapse_hpc")

    # read out time used for building
    BuildEdgeTime = time.time() - tic

    logger.log(str(BuildEdgeTime) + " # build_edge_time")
    logger.log(str(memory_thisjob()) + " # virt_mem_after_edges")

    return None


def run_simulation():
    """Performs a simulation, including network construction."""

    # open log file
    with Logger(params["log_file"]) as logger:
        nest.ResetKernel()
        nest.set_verbosity(M_INFO)

        logger.log(str(memory_thisjob()) + " # virt_mem_0")

        sr = build_network(logger)

        tic = time.time()

        nest.Simulate(params["presimtime"])

        PreparationTime = time.time() - tic

        logger.log(str(memory_thisjob()) + " # virt_mem_after_presim")
        logger.log(str(PreparationTime) + " # presim_time")

        tic = time.time()

        nest.Simulate(params["simtime"])

        SimCPUTime = time.time() - tic

        logger.log(str(memory_thisjob()) + " # virt_mem_after_sim")
        logger.log(str(SimCPUTime) + " # sim_time")

        print(nest.kernel_status)


def compute_rate(sr):
    """Compute local approximation of average firing rate.

    This approximation is based on the number of local nodes, number of local spikes and total time. Since this also
    considers devices, the actual firing rate is usually underestimated.
    """

    n_local_spikes = nest.local_spike_counter  # sr.n_events
    n_local_neurons = brunel_params["Nrec"]
    simtime = params["simtime"]
    return 1.0 * n_local_spikes / (n_local_neurons * simtime) * 1e3


def memory_thisjob():
    """Wrapper to obtain current memory usage"""
    nest.ll_api.sr("memory_thisjob")
    return nest.ll_api.spp()


def lambertwm1(x):
    """Wrapper for LambertWm1 function"""
    # Using scipy to mimic the gsl_sf_lambert_Wm1 function.
    return sp.lambertw(x, k=-1 if x < 0 else 0).real


class Logger:
    """Logger context manager used to properly log memory and timing information from network simulations."""

    def __init__(self, file_name):
        # copy output to cout for ranks 0..max_rank_cout-1
        self.max_rank_cout = 5
        # write to log files for ranks 0..max_rank_log-1
        self.max_rank_log = 30
        self.line_counter = 0
        self.file_name = file_name

    def __enter__(self):
        if nest.Rank() < self.max_rank_log:
            # convert rank to string, prepend 0 if necessary to make
            # numbers equally wide for all ranks
            rank = "{:0" + str(len(str(self.max_rank_log))) + "}"
            fn = "{fn}_{rank}.dat".format(fn=self.file_name, rank=rank.format(nest.Rank()))

            self.f = open(fn, "w")

        return self

    def log(self, value):
        if nest.Rank() < self.max_rank_log:
            line = "{lc} {rank} {value} \n".format(lc=self.line_counter, rank=nest.Rank(), value=value)
            self.f.write(line)
            self.line_counter += 1

        if nest.Rank() < self.max_rank_cout:
            print(str(nest.Rank()) + " " + value + "\n", file=sys.stdout)
            print(str(nest.Rank()) + " " + value + "\n", file=sys.stderr)

    def __exit__(self, exc_type, exc_val, traceback):
        if nest.Rank() < self.max_rank_log:
            self.f.close()


if __name__ == "__main__":
    run_simulation()
