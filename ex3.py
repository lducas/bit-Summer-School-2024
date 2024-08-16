import matplotlib.pyplot as plt
import numpy as np
from numpy import zeros, array
from sol2 import LLL
from sys import argv
from math import sqrt
from numpy import zeros, array
from numpy.linalg import norm
from random import randint
from time import time

from sol1 import Gram_Schmidt_orth, nearest_plane
from sol2 import lagrange_reduce, size_reduce
from sol2 import LLL as slow_LLL
from ver2 import verify_ex3 as verify_ex1

#############
# Exercise 1:
# Make LLL faster. Note that we could do several non-overlapping Lagrange 
# steps between each call to Gram-Schmidt and Size-Reduction
#############

gamma_2 = sqrt(4/3)

def fast_LLL(B, epsilon=0.01, anim=True):
	pass

if argv[0].endswith("sol3.py") or argv[0].endswith("ex3.py"):
	verify_ex1(fast_LLL, Gram_Schmidt_orth)

print("Comparing speed ... ")

def gen_basis(n, q):
	h = n//2
	B = np.identity(n, dtype=int)
	for i in range(0, h):
		B[i, i] = q
		for j in range(0, h):
			B[h+j, i] = randint(0, q)
	return B

def timing(red, n, q):
	B = gen_basis(n, q)
	startT = time()
	L = red(B)
	if L is not None:
		L = list(L)
	T = time() - startT
	return T

ns = range(10, 31, 2)
q = 99

Tfast = [timing(fast_LLL, n, q) for n in ns]
Tslow = [timing(slow_LLL, n, q) for n in ns]

plt.plot(ns, Tfast, label="Fast LLL")
plt.plot(ns, Tslow, label="slow LLL")
plt.xlabel('dimension')
plt.ylabel('Time in s')
plt.legend()
plt.show()


#############
# Exercise 2 (probably to bring back home ...)
# The implementation so far remains quite naive. Serious implementation
# never really compute the Gram-Schmidt basis directly, but rather R
# part of the QR decomposition, and update it locally rather than 
# recomputing it. See for example the pseudo-code given here:
# https://www.ionica.nl/wp-content/uploads/2019/04/the-LLL-Algorithms.pdf#page=160
#############
