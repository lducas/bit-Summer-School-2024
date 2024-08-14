import numpy as np
from numpy import zeros, array
from random import randint
from math import sqrt

class ExerciseFailed(Exception):
    pass

def iarray(x):
	return array(x, dtype=int)

bases = [
array([[1],]),
array([[2],]),
array([[200],]),
iarray([[4,0],[3,1]]),
iarray([[50,0],[0,1]]),
iarray([[50,0],[40,1]]),
iarray([[50, -30, 14], [0, 200, -4], [350, 600, -12]]),
iarray([[50, 33, -30, 14], [200, 40, -45, -4], [1, 35, 0, -1200], [-15, 3, 8, -7]])
]


def verify_ex1(LR):
	for i in range(100):
		B = iarray([[999, 0], [randint(0, 999), 1]])
		B_ = np.copy(B)
		U = LR(B)
		if not isinstance(U[0,0], np.int64):
			raise ExerciseFailed("Exercise 1 Failed : U is not an integer matrix")
		if not round(abs(np.linalg.det(U))) == 1:
			raise ExerciseFailed("Exercise 1 Failed : U is not unimodular (det not equal to +/- 1)")
		if not np.all(B == U.dot(B_)):
			raise ExerciseFailed("Exercise 1 Failed : B_out not equal to U * B_in")
		if np.linalg.norm(B[1]) < np.linalg.norm(B[0]):
			raise ExerciseFailed("Exercise 1 Failed : Not Lagrange-reduced")
		if abs(B[0].dot(B[1]) / B[0].dot(B[0])) > .500001:
			raise ExerciseFailed("Exercise 1 Failed : Not Lagrange-reduced")

	print("Exercise 1 succeeded")


def in_lattice(B, t):
	x = np.linalg.solve(B.transpose(), t)
	xr = np.round(x)
	return np.allclose(x, xr)

def gen_same_lat(B, C):
	for b in B:
		if not in_lattice(C, b):
			return False
	for c in C:
		if not in_lattice(B, c):
			return False
	return True

def orth_proj(x, y):
	""" Returns the orthogonal projection of vector x orthognally to vector y"""
	return x - (x.dot(y)/(y.dot(y)))*y

def Gram_Schmidt_orth(B):
	""" Output the Gram-Schmidt Orthogonalization of a matrix B"""

	# Makes a copy of B, but change the type to float
	Bs = array(B, dtype=float)
	# Get the dimension of the square basis
	n,_ = B.shape
	for i in range(n):
		for j in range(i):
			Bs[i] = orth_proj(Bs[i], Bs[j])
	return Bs

def verify_ex2(SR,GSO):
	for B_ in bases:
		n,_ = B_.shape
		Bs = GSO(B_)
		B = np.copy(B_)
		SR(B, Bs)
		if not gen_same_lat(B, B_):
			raise ExerciseFailed("Exercise 2 Failed : Output basis does not generate the same lattice as input")
		for i in range(n):
			for j in range(i):
				if abs(B[i].dot(Bs[j]) / Bs[j].dot(Bs[j])) > .500001:
					raise ExerciseFailed("Exercise 2 Failed : Not Size-Reduced")
	print("Exercise 2 succeeded")

gamma_2 = sqrt(4/3)


def verify_ex3(LLL, GSO):
	epsilon=0.01
	for n in range(2, 20):
		q = 999999
		B_ = np.identity(n, dtype=int)
		B_[0, 0] = q
		for i in range(1, n):
			B_[i, 0] = randint(0, q)

		B = np.copy(B_)
		X = LLL(B,epsilon=epsilon, anim=False)
		X = list(X) if X is not None else None
		if not gen_same_lat(B, B_):
			raise ExerciseFailed("Exercise 3 Failed : Output basis does not generate the same lattice as input")
		Bs = GSO(B)
		for i in range(n-1):
			if np.linalg.norm(Bs[i]) > (gamma_2 + epsilon) * np.linalg.norm(Bs[i+1]):
				raise ExerciseFailed("Exercise 3 Failed : Not Weak-LLL reduced")
		for i in range(n):
			for j in range(i):
				if abs(B[i].dot(Bs[j]) / Bs[j].dot(Bs[j])) > .500001:
					raise ExerciseFailed("Exercise 3 Failed : Not Size-Reduced")

	print("Exercise 3 succeeded")	


