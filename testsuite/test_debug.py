# import nest
#
# nest.ResetKernel()
#
# nest.print_time = True
# nest.local_num_threads = 1
# # nest.use_compressed_spikes = True
#
# n1, n2 = nest.Create("iaf_psc_alpha_ax_delay", 2)
#
# nest.Connect(n1, n2, syn_spec={"synapse_model": "stdp_pl_synapse_hom_ax_delay", "weight": 10., "delay": 1.})
#
# poisson_generator = nest.Create("poisson_generator", 1, params={"rate": 2000.})
#
# nest.Connect(poisson_generator, n1, syn_spec={"synapse_model": "static_synapse", "weight": 9999.})
#
# nest.Simulate(5)

import nest


neurons = nest.Create("iaf_psc_alpha_ax_delay", 4)
nest.Connect(neurons, neurons,
             {'rule': 'fixed_indegree', 'indegree': 3,
              'allow_autapses': False, 'allow_multapses': True},
             {'synapse_model': 'static_synapse'})
nest.Connect(neurons, neurons,
             {'rule': 'fixed_indegree', 'indegree': 3,
              'allow_autapses': False, 'allow_multapses': True},
             {'synapse_model': 'stdp_pl_synapse_hom_ax_delay'})

nest.Simulate(100)
nest.Simulate(100)
