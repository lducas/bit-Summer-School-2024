import numpy as np
from numpy import zeros, array

class ExerciseFailed(Exception):
    pass

def iarray(x):
	return array(x, dtype=int)


bases = [
iarray([[1],]),
iarray([[2],]),
iarray([[200],]),
iarray([[1,0],[0,1]]),
iarray([[50,0],[0,1]]),
iarray([[50,0],[20,1]]),
iarray([[50, -30, 14], [0, 20, -4], [35, 0, -12]]),
iarray([[50, 33, -30, 14], [0, 20, 4, -4], [1, 35, 0, -12], [-15, 3, 8, -7]])
]

def verify_ex1(SR):
	for B in bases:
		n,_ = B.shape
		for it in range(20):
			x = 50 * (np.random.rand(n) - np.random.rand(n))
			xr = np.round(x)
			# Skip borderline cases to tolerate numerical stability errors
			if max(abs(x - xr)) > .499:
				pass
			t = x.dot(B)
			res = SR(B, t)
			if res is None:
				raise ExerciseFailed("Exercise 1 Failed : nothing returned")
			if not np.allclose(SR(B, t), xr.dot(B)):
				raise ExerciseFailed("Exercise 1 Failed")
	print("Exercise 1 succeeded")


def verify_ex2a(OP):
	for n in [2, 3, 5, 10, 50]:
		for it in range(20):
			x = 50 * (np.random.rand(n) - np.random.rand(n))
			y = 50 * (np.random.rand(n) - np.random.rand(n))
			z = OP(x, y)
			test1 = [z.dot(y), z.dot(x)]
			test2 = [0, z.dot(z)]
		if not np.allclose(test1, test2):
			raise ExerciseFailed("Exercise 2a Failed")
	print("Exercise 2a succeeded")

def verify_ex2b(GSO):
	for B in bases:
		Bs = GSO(B)
		n,_ = B.shape
		if not np.allclose(np.linalg.det(B), np.linalg.det(Bs)):
			raise ExerciseFailed("Exercise 2b Failed")

		D = Bs.dot(Bs.transpose())
		T = B.dot(Bs.transpose())

		if not np.allclose(np.diagonal(D), np.diagonal(T)):
			raise ExerciseFailed("Exercise 2b Failed")

		for i in range(n):
			for j in range(i):
				if not ((abs(D[j, i]) < 1e-8) and (abs(D[i, j]) < 1e-8) and (abs(T[j, i]) < 1e-8)):
					raise ExerciseFailed("Exercise 2b Failed")
	print("Exercise 2b succeeded")

def verify_ex3(IL, GSO, NP):
	for B in bases:
		n, _ = B.shape
		Bs = GSO(B)

		for it in range(20):
			x = 50 * (np.random.rand(n) - np.random.rand(n))
			t = x.dot(B)
			v = NP(B, Bs, np.copy(t))
			if not IL(B, v):
				raise ExerciseFailed("Exercise 3 Failed (not in the lattice)")
			for bs in Bs:
				if abs((t - v).dot(bs)) > .501 * bs.dot(bs):
					raise ExerciseFailed("Exercise 3 Failed (not the right lattice point)")
	print("Exercise 3 succeeded")

