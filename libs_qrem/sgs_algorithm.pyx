import numpy as np
from pprint import pprint
from typing import Union, List
import copy
import heapq
from heapq import heappush, heappop


class priority_queue(object):
    """
    Priority queue wrapper which enables to compare the specific elements of container as keys.
    """

    def __init__(self, key_index=0):
        """
        Arguments
            key_index: the index of elements as keys
        """
        self.key = lambda item: item[key_index]
        self.index = 0
        self.data = []

    def size(self):
        """
        Return the size of heap
        """
        return len(self.data)

    def push(self, item):
        """
        Push a container to heap list
        
        Arguments
            item: container
        """
        heapq.heappush(self.data, (self.key(item), self.index, item))
        self.index += 1

    def pop(self):
        """
        Pop the smallest element of heap
        """
        if len(self.data) > 0:
            return heapq.heappop(self.data)[2]
        else:
            return None

    def top(self):
        """
        Refer the smallest element of heap
        """
        if self.size() > 0:
            return self.data[0][2]
        else:
            return None 

def sgs_algorithm_py(x, make_sum_to_one = False):
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
    negative_accumulator = 0

    if make_sum_to_one:
        for state_idx in x:  # O(N) time
            if x[state_idx] > 0:
                pq.push((state_idx, x[state_idx]))  # O(log(N)) time
                sum_of_x += x[state_idx]
        negative_accumulator = 1 - sum_of_x
    else:
        for state_idx in x:  # O(N) time
            pq.push((state_idx, x[state_idx]))  # O(log(N)) time
            sum_of_x += x[state_idx]

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