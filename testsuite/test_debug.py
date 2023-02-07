import nest

nest.ResetKernel()

nest.local_num_threads = 3

n1, n2, _ = nest.Create("iaf_psc_alpha_ax_delay", 3)

nest.Connect(n1, n2, syn_spec={"synapse_model": "stdp_pl_synapse_hom_ax_delay", "weight": 10., "delay": 1.})

poisson_generator = nest.Create("poisson_generator", 1, params={"rate": 20.})

nest.Connect(poisson_generator, n1, syn_spec={"synapse_model": "static_synapse", "weight": 9999.})

nest.Simulate(100)
