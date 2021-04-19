import numpy as np

def compare(a, b):
    return a - b


def lin_search(val, arr):
    return bin_search(val, arr)
    i = 0
    while i < len(arr):
        if (compare(val, arr[i]) <= 0):
            return i
        i += 1
    return i

def bin_search(val, arr):
    l = 0
    r = len(arr) - 1
    while (l <= r):
        m = l + (r - l) // 2
        if (compare(arr[m], val) == 0):
            return m
        if (compare(arr[m], val) < 0):
            l = m + 1
        else:
            r = m - 1
    return l

r = np.random.rand(10000)
r.sort()
x_arr = np.random.rand(3000)
for x in x_arr:
	if lin_search(x, r) != bin_search(x, r):
		print("An error occured for {} and {} at value of {}".format(lin_search(x, r), bin_search(x, r), x))
