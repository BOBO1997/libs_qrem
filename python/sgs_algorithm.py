from typing import Union
import numpy as np
from pprint import pprint

import priority_queue
from priority_queue import priority_queue

def sgs_algorithm(x: Union[dict, np.array], silent = False) -> dict:
    """
    The negative cancellation algorithm by Smolin, Gambetta, and Smith.
    O(NlogN) time, O(N) memory to the size of x: N

    Arguments
        x: sum 1 probability vecotor with negative values
    Returns
        x_tilde: physically correct probability vector
    """

    # compute the number and the sum of negative values
    pq = priority_queue(key_index=1)
    sum_of_x = 0
    for state_idx in x:  # O(N) time
        if x[state_idx] > 0:
            pq.push((state_idx, x[state_idx]))  # O(log(N)) time
            sum_of_x += x[state_idx]
    if not silent:
        print("number of positive values: ", pq.size())

    negative_accumulator = 1 - sum_of_x
    if negative_accumulator >= 0:
        print("accumulator is positive, we might even ignoring the necessal positive values.")

    while pq.size() > 0:  # O(N) time
        _, x_hat_i = pq.top()
        if x_hat_i + negative_accumulator / pq.size() < 0:
            negative_accumulator += x_hat_i
            _, _ = pq.pop()  # O(log(N)) time
            continue
        else:
            break

    x_tilde = {}
    denominator = pq.size()
    while pq.size() > 0:  # O(N) time
        state_idx, x_hat_i = pq.pop()  # O(log(N))
        x_tilde[state_idx] = x_hat_i + negative_accumulator / denominator

    return x_tilde