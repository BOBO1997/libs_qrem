from typing import Union, List
import copy
import numpy as np
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