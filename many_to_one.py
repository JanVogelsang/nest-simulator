import sys
import nest


tau_syn = 0.32582722403722841
neuron_params = {
    "E_L": 0.0,
    "C_m": 250.0,
    "tau_m": 10.0,
    "t_ref": 0.5,
    "V_th": 20.0,
    "V_reset": 0.0,
    "tau_syn_ex": tau_syn,
    "tau_syn_in": tau_syn,
    "tau_minus": 30.0,
    "V_m": 5.7
}

stdp_params = {
    "delay": {"dend": 0.2, "mixed": 0.1, "ax": 0.}[sys.argv[1]],
    "axonal_delay": {"dend": 0., "mixed": 0.1, "ax": 0.2}[sys.argv[1]],
    # "delay": 0. if sys.argv[1] == "ax" else 2.5,
    # "axonal_delay": 5.0 if sys.argv[1] == "ax" else 2.5,
    "alpha": 0.0513,
    "mu": 0.4,
    "tau_plus": 15.0,
    "weight": 45.
}

T = 10.
NE = 100
eta = 0.


def run():
    stdp_params['lambda'] = 0.

    nest.ResetKernel()
    nest.local_num_threads = 1
    nest.rng_seed = 42
    nest.set_verbosity("M_ERROR")

    # E_ext = nest.Create("poisson_generator", 1, {"rate": eta * 1000.0})
    E_pg = nest.Create("poisson_generator", 1, params={"rate": 10.0})
    I_pg = nest.Create("poisson_generator", 1, params={"rate": 10.0 * NE / 5})
    E_neurons = nest.Create("parrot_neuron", NE)
    I_neurons = nest.Create("parrot_neuron", 1)
    post_neuron = nest.Create("iaf_psc_alpha_ax_delay", 1, params=neuron_params)
    sr = nest.Create("spike_recorder")

    nest.SetDefaults("stdp_pl_synapse_hom_ax_delay", stdp_params)

    nest.Connect(E_pg, E_neurons, syn_spec={"delay": 0.1})
    nest.Connect(I_pg, I_neurons, syn_spec={"delay": 0.1})
    nest.Connect(E_neurons, post_neuron, syn_spec={"synapse_model": "stdp_pl_synapse_hom_ax_delay"})
    nest.Connect(I_neurons, post_neuron, syn_spec={"weight": 45. * -5, "delay": 0.1})
    # nest.Connect(E_ext, post_neuron, syn_spec={"weight": 45.})
    nest.Connect(post_neuron, sr, syn_spec={"delay": 0.1})

    nest.Simulate(T)
    print(nest.local_spike_counter)

    return sr, post_neuron


run()
