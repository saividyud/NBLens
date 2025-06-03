import sys
import numpy as np

class Test:
    def __init__(self):
        self.arr = [1, 2, 3]
        self.arr2 = np.array([5, 6, 7])

hello = Test()

for var in vars(hello).values():
    print(sys.getsizeof(var))