import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import zeros, array
from numpy.linalg import norm
from sol1 import Gram_Schmidt_orth, nearest_plane
from ver2 import verify_ex1, verify_ex2, verify_ex3
from test_bases import B2, B4, B24
from math import sqrt, log
from random import randint
from sys import argv

############
# Exercise 1
# Implement Lagrange reduction Algorithm. The algorithm should output
# the 2x2 Unimodular transformation matrix U. To do so, start with U
# being the identity matrix, and apply the same transformation to U
# than to B all along the algorithm.
############

def lagrange_reduce(B):
	""" Given a basis with two rows as input, apply lagrange reduction to B (modified 
	in place). 
	Also output the transformation matrix U, sending the initial basis to the final basis
	"""
	U = np.identity(2, dtype=int)
	first = True
	while first or np.linalg.norm(B[1]) < np.linalg.norm(B[0]):
		B[[0, 1]] = B[[1, 0]]
		U[[0, 1]] = U[[1, 0]]
		first = False
		k = int(round(B[0].dot(B[1]) / B[0].dot(B[0])))
		B[1] -= k * B[0]
		U[1] -= k * U[0]

	return U

if argv[0].endswith("sol2.py") or argv[0].endswith("ex2.py"):
	verify_ex1(lagrange_reduce)


############
# Exercise 2
# Implement a the Size-Reduction Algorithm (in place) on a basis input B. 
# As a  by-product, provide as output the Gram-Schimdt orthogonalisation of B.
#
# Note: In the lecture notes, the NearestPlane is applied to a projection
# pi_C'(b_n). This projection is unecessary (but make the proof nicer).
# Ignore it in your implementation.
# The lecture note also give the algorithm in a reccursive manner, but for
# implementation in python, an iterative version will be much simpler
############

def size_reduce(B, Bs):
	""" Given a basis B and its Gram-Schmidt Bs, 
	apply size reduction to B (modified in place). 
	No return value.
	"""
	# Get the dimension number of vectors in the basis
	n,_ = B.shape
	for i in range(n):
		v = nearest_plane(B[:i], Bs[:i], B[i])
		B[i] -= v
	return

if argv[0].endswith("sol2.py") or argv[0].endswith("ex2.py"):
	verify_ex2(size_reduce, Gram_Schmidt_orth)

############
# Exercise 3
# Implement the LLL algorithm. Beware that The Gram-Schmidt basis needs to be 
# updated after every modification of B !
# 
# Note: at some point, you will need to apply lagrange to the basis 2-dimensional
# basis [pi_i(b_i), pi_i(b_{i+1})]. This can be re-constructed rather cheaply from
# B and B* by noting that:
# pi_i(b_i) = b*_i and 
# pi_i(b_{i+1}) = b*_{i+1} + (<b_i+1, b*_i> / ||b*_i||^2) * b*_i
############

gamma_2 = sqrt(4/3)

def LLL(B, epsilon=0.01, anim=True):
	n,_ = B.shape
	Bs = Gram_Schmidt_orth(B)
	size_reduce(B, Bs)

	while True:
		yield [log(norm(x)) for x in Bs]
		for i in range(n):
			if i==n-1:
				return
			if norm(Bs[i]) > (gamma_2+epsilon) * norm(Bs[i+1]):
				break
		# We have found an index i to be Lagrange-reduced
		x1 = np.copy(Bs[i])
		x2 = Bs[i+1] + (B[i+1].dot(Bs[i]) / (Bs[i].dot(Bs[i]))) * Bs[i]
		X = array([x1, x2])
		U = lagrange_reduce(X)
		B[i:i+2] = U.dot(B[i:i+2])
		Bs = Gram_Schmidt_orth(B)
		size_reduce(B, Bs)

if argv[0].endswith("sol2.py") or argv[0].endswith("ex2.py"):
	verify_ex3(LLL, Gram_Schmidt_orth)

#############
# Exercise 4
# Add the following line to the beginning of your LLL loop above, and enjoy the show
#
# yield [log(norm(x)) for x in Bs]
#
#############

def anim_LLL(n, q):
	B = np.identity(n, dtype=int)
	B[0, 0] = q
	for i in range(1, n):
		B[i, 0] = randint(0, q)

	try:
		data = list(LLL(B, anim=True))
	except:
		print("Exercise 4 Failed")
		exit()

	for i in range(40):
		data.append(data[-1])

	fig, ax = plt.subplots()
	line2 = ax.plot(range(n), data[0], label="basis profile")[0]
	ax.set(xlim=[0, n], ylim=[0, log(q)/4])
	ax.legend()

	def update_anim(frame):
		line2.set_ydata(data[frame])
		return line2

	ani = animation.FuncAnimation(fig=fig,func=update_anim, frames=len(data), interval=50)
	plt.show()
	#ani.save("lll.gif", writer="pillow")


if argv[0].endswith("sol2.py") or argv[0].endswith("ex2.py"):
	anim_LLL(20, 999999)