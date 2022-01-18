import numpy as np
from libs_qrem import DeltaFilter, LeastNormFilter, MooneyEtalFilter, NationEtalFilter
from dummy_data12 import n, cal_matrices, hist
from pprint import pprint

# mitigator = DeltaFilter(n, cal_matrices)
mitigator = MooneyEtalFilter(n, cal_matrices)


hist = mitigator.apply(hist, threshold = 0.1)
pprint("sums")
pprint(mitigator.sum_of_x())
pprint(mitigator.sum_of_x_hat())
pprint(mitigator.sum_of_x_tilde())
pprint("matrices")
pprint(np.asarray(mitigator.reduced_A())[:5, :5])
pprint(np.asarray(mitigator.normalized_reduced_A())[:5, :5])
pprint(np.asarray(mitigator.reduced_inv_A())[:5, :5])
pprint("norms")
pprint(mitigator.exact_one_norm_of_inv_reduced_A())
pprint(mitigator.exact_one_norm_of_reduced_inv_A())
pprint(mitigator.iterative_one_norm_of_inv_reduced_A())
pprint("vectors")
pprint(mitigator.indices_to_keys_vector()[:5])
pprint(mitigator.x_s()[:5])
pprint(mitigator.x_hat()[:5])
pprint(mitigator.x_tilde()[:5])
pprint("others")
pprint(mitigator.times())
pprint(mitigator.expval())
pprint(mitigator.mitigation_stddev())

dct = {# "exact_one_norm_of_reduced_inv_A": mitigator.exact_one_norm_of_reduced_inv_A(),
       "mitigated_hist": list(mitigator.mitigated_hist().keys())[:10],
       "x_s": mitigator.x_s()[:10],
       "x_hat": mitigator.x_hat()[:10],
       "x_tilde": mitigator.x_tilde()[:10],
       "sum_of_x": mitigator.sum_of_x(),
       "sum_of_x_hat": mitigator.sum_of_x_hat(),
       "sum_of_x_tilde": mitigator.sum_of_x_tilde(),
       "indices_to_keys_vector": mitigator.indices_to_keys_vector(),
       "times": mitigator.times(),
       "expval": mitigator.expval(),
       "mitigation_stddev": mitigator.mitigation_stddev(norm_type="exact"),
}

print("abort?????")

# pprint(dct)

print("segfo?????????")
# del mitigator
print("segfo?????????")
