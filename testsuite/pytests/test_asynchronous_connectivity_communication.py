import nest

neurons = nest.Create('iaf_psc_alpha', 100)
nest.Connect(neurons, neurons, conn_spec={'rule': 'fixed_indegree', 'indegree': 10}, syn_spec={"synapse_model": "static_synapse"})
nest.Prepare()
