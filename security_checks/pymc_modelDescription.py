import pymc3

model = pymc3.Model()

obs1, obs2, obs3 = 2, 4, 5

d = {'a1': {'name': 'a1', 'lower': 0, 'upper': 10},
     'a2': {'name': 'a2', 'lower': 0, 'upper': 10},
     'a3': {'name': 'a3', 'lower': 0, 'upper': 10}}

relations = {'b1': {'a1': 2, 'a2': 1}, 'b2': {'a1': 4}, 'b3': {'a3': 5}}

correspondances_dict = {'b1': {'random_var_name': 'm1', 'observation': obs1},
                        'b2': {'random_var_name': 'm2', 'observation': obs2},
                        'b3': {'random_var_name': 'm3', 'observation': obs3}}

with model:
    priors = {prior_name: pymc3.Uniform(prior_name,
                                        lower=d[prior_name]['lower'],
                                        upper=d[prior_name]['upper']) for prior_name in list(d.keys())}

    intermediate_vars = {intermediate_var: sum([relations[intermediate_var][prior_name] * priors[prior_name]
                                                for prior_name in list(relations[intermediate_var].keys())])
                         for intermediate_var in list(relations.keys())}

    observed_vars = {correspondances_dict[intermediate_var]['random_var_name']:
                         pymc3.Normal(correspondances_dict[intermediate_var]['random_var_name'],
                                      mu=intermediate_vars[intermediate_var],
                                      sd=0.1,
                                      observed=correspondances_dict[intermediate_var]['observation'])
                     for intermediate_var in list(intermediate_vars.keys())}

    trace = pymc3.sample(1000)
