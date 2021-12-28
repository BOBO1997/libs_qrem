import numpy as np
from libs_qrem import DeltaFilter, LeastNormFilter, MooneyEtalFilter, NationEtalFilter
from dummy_data import n, cal_matrices, hist
from pprint import pprint

delta_mitigator = DeltaFilter(n, cal_matrices)
# slsqp_mitigator = SLSQPFilter(n, cal_matrices)
# lnp_mitigator = LeastNormFilter(n, cal_matrices)
# mooney_mitigator = MooneyEtalFilter(n, cal_matrices)
# nation_mitigator = NationEtalFilter(n, cal_matrices)


pprint("proposed (delta)")
delta_hist = delta_mitigator.apply(hist)
pprint("sums")
pprint(delta_mitigator.sum_of_x())
pprint(delta_mitigator.sum_of_x_hat())
pprint(delta_mitigator.sum_of_x_tilde())
pprint("matrices")
pprint(np.asarray(delta_mitigator.reduced_A())[:5, :5])
pprint(np.asarray(delta_mitigator.normalized_reduced_A())[:5, :5])
pprint(np.asarray(delta_mitigator.reduced_inv_A())[:5, :5])
pprint("norms")
pprint(delta_mitigator.exact_one_norm_of_inv_reduced_A())
pprint(delta_mitigator.exact_one_norm_of_reduced_inv_A())
pprint(delta_mitigator.iterative_one_norm_of_inv_reduced_A())
pprint("vectors")
pprint(delta_mitigator.indices_to_keys_vector()[:5])
pprint(delta_mitigator.x_s()[:5])
pprint(delta_mitigator.x_hat()[:5])
pprint(delta_mitigator.x_tilde()[:5])
pprint("others")
pprint(delta_mitigator.times())
pprint(delta_mitigator.expval())
pprint(delta_mitigator.mitigation_stddev())

dct = {# "exact_one_norm_of_reduced_inv_A": delta_mitigator.exact_one_norm_of_reduced_inv_A(),
       "mitigated_hist": list(delta_mitigator.mitigated_hist().keys())[:10],
       "x_s": delta_mitigator.x_s()[:10],
       "x_hat": delta_mitigator.x_hat()[:10],
       "x_tilde": delta_mitigator.x_tilde()[:10],
       "sum_of_x": delta_mitigator.sum_of_x(),
       "sum_of_x_hat": delta_mitigator.sum_of_x_hat(),
       "sum_of_x_tilde": delta_mitigator.sum_of_x_tilde(),
       "indices_to_keys_vector": delta_mitigator.indices_to_keys_vector(),
       "times": delta_mitigator.times(),
       "expval": delta_mitigator.expval(),
       "mitigation_stddev": delta_mitigator.mitigation_stddev(norm_type="exact"),
}

print("abort?????")

# pprint(dct)

print("segfo?????????")
# del delta_mitigator
print("segfo?????????")
