import numpy as np
from libs_qrem import DeltaFilter, SLSQPFilter, LeastNormFilter, MooneyEtalFilter, NationEtalFilter
from dummy_data import n, cal_matrices, hist
from pprint import pprint

delta_mitigator = DeltaFilter(n, cal_matrices)
slsqp_mitigator = SLSQPFilter(n, cal_matrices)
lnp_mitigator = LeastNormFilter(n, cal_matrices)
mooney_mitigator = MooneyEtalFilter(n, cal_matrices)
nation_mitigator = NationEtalFilter(n, cal_matrices)

pprint("proposed (delta)")
delta_hist = delta_mitigator.apply(hist)
pprint("sums")
pprint(delta_mitigator.sum_of_x())
pprint(delta_mitigator.sum_of_x_hat())
pprint(delta_mitigator.sum_of_x_tilde())
pprint("matrices")
pprint(np.asarray(delta_mitigator.reduced_A())[:10, :10])
pprint(np.asarray(delta_mitigator.normalized_reduced_A())[:10, :10])
pprint(np.asarray(delta_mitigator.reduced_inv_A())[:10, :10])
pprint("norms")
pprint(delta_mitigator.exact_one_norm_of_inv_reduced_A())
pprint(delta_mitigator.exact_one_norm_of_reduced_inv_A())
pprint(delta_mitigator.iterative_one_norm_of_inv_reduced_A())
pprint("vectors")
pprint(delta_mitigator.x_s()[:10])
pprint(delta_mitigator.x_hat()[:10])
pprint(delta_mitigator.x_tilde()[:10])
pprint("others")
pprint(delta_mitigator.indices_to_keys_vector()[:10])
pprint(delta_mitigator.times())
pprint(delta_mitigator.expval())
pprint(delta_mitigator.mitigation_overhead())

"""
print("proposed (slsqp)")
slsqp_hist = slsqp_mitigator.apply(hist)
pprint(slsqp_mitigator.sum_of_x())
pprint(slsqp_mitigator.sum_of_x_hat())
pprint(slsqp_mitigator.sum_of_x_tilde())
pprint(slsqp_mitigator.reduced_A()[:10][:10])
pprint(slsqp_mitigator.normalized_reduced_A()[:10][:10])
pprint(slsqp_mitigator.reduced_inv_A()[:10][:10])
pprint(slsqp_mitigator.exact_one_norm_of_inv_reduced_A())
pprint(slsqp_mitigator.exact_one_norm_of_reduced_inv_A())
pprint(slsqp_mitigator.iterative_one_norm_of_inv_reduced_A())
pprint(slsqp_mitigator.x_s()[:10])
pprint(slsqp_mitigator.x_hat()[:10])
pprint(slsqp_mitigator.x_tilde()[:10])
pprint(slsqp_mitigator.indices_to_keys_vector()[:10])
pprint(slsqp_mitigator.times())
pprint(slsqp_mitigator.expval())
pprint(slsqp_mitigator.mitigation_overhead())


print("proposed (least norm)")
lnp_hist = lnp_mitigator.apply(hist)
pprint(lnp_mitigator.sum_of_x())
pprint(lnp_mitigator.sum_of_x_hat())
pprint(lnp_mitigator.sum_of_x_tilde())
pprint(lnp_mitigator.reduced_A()[:10][:10])
pprint(lnp_mitigator.normalized_reduced_A()[:10][:10])
pprint(lnp_mitigator.reduced_inv_A()[:10][:10])
pprint(lnp_mitigator.exact_one_norm_of_inv_reduced_A())
pprint(lnp_mitigator.exact_one_norm_of_reduced_inv_A())
pprint(lnp_mitigator.iterative_one_norm_of_inv_reduced_A())
pprint(lnp_mitigator.x_s()[:10])
pprint(lnp_mitigator.x_hat()[:10])
pprint(lnp_mitigator.x_tilde()[:10])
pprint(lnp_mitigator.indices_to_keys_vector()[:10])
pprint(lnp_mitigator.times())
pprint(lnp_mitigator.expval())
pprint(lnp_mitigator.mitigation_overhead())


print("Mooney et al.")
mooney_hist = mooney_mitigator.apply(hist)
pprint(mooney_mitigator.sum_of_x())
pprint(mooney_mitigator.sum_of_x_hat())
pprint(mooney_mitigator.sum_of_x_tilde())
pprint(mooney_mitigator.reduced_A()[:10][:10])
pprint(mooney_mitigator.normalized_reduced_A()[:10][:10])
pprint(mooney_mitigator.reduced_inv_A()[:10][:10])
pprint(mooney_mitigator.exact_one_norm_of_inv_reduced_A())
pprint(mooney_mitigator.exact_one_norm_of_reduced_inv_A())
pprint(mooney_mitigator.iterative_one_norm_of_inv_reduced_A())
pprint(mooney_mitigator.x_s()[:10])
pprint(mooney_mitigator.x_hat()[:10])
pprint(mooney_mitigator.x_tilde()[:10])
pprint(mooney_mitigator.indices_to_keys_vector()[:10])
pprint(mooney_mitigator.times())
pprint(mooney_mitigator.expval())
pprint(mooney_mitigator.mitigation_overhead())

print("Nation et al.")
nation_hist = nation_mitigator.apply(hist)
pprint(nation_mitigator.sum_of_x_hat())
pprint(nation_mitigator.sum_of_x_tilde())
pprint(nation_mitigator.reduced_A()[:10][:10])
pprint(nation_mitigator.normalized_reduced_A()[:10][:10])
pprint(nation_mitigator.reduced_inv_A()[:10][:10])
pprint(nation_mitigator.exact_one_norm_of_inv_reduced_A())
pprint(nation_mitigator.exact_one_norm_of_reduced_inv_A())
pprint(nation_mitigator.iterative_one_norm_of_inv_reduced_A())
pprint(nation_mitigator.mitigated_hist()[:10])
pprint(nation_mitigator.x_s()[:10])
pprint(nation_mitigator.x_hat()[:10])
pprint(nation_mitigator.x_tilde()[:10])
pprint(nation_mitigator.indices_to_keys_vector()[:10])
pprint(nation_mitigator.times())
pprint(nation_mitigator.expval())
pprint(nation_mitigator.mitigation_overhead())
"""