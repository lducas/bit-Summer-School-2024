import numpy as np
import matplotlib.pyplot as plt
from numpy import zeros, array
from numpy.linalg import norm
from sys import argv
from ver1 import verify_ex1, verify_ex2a, verify_ex2b, verify_ex3
from test_bases import B2, B4, B24

# You never need to unroll vector of matrix operation:
# - x+y is the simple way to add vectors
# - A.dot(B) to multiply matrices or vectors, including inner products
# - Most linear algebraic operations are provided in numpy, or numpy.linalg
#   https://numpy.org/doc/stable/reference/routines.linalg.html
# - coefficient-wise operation are also provided, such as np.round(x) to round each coordinate
# - Doing so will be simpler, more readable, and much faster than writing a for loop

#   ^
#  /|\    When implementing, it is more common to use *row* notations (and faster...)
# / | \   Do so. And beware on which side you are multiplying vectors and matrices !
#/__*__\

# The exercises comprises of function to be implemented (except Exercise 0):
# Replace the keyword "pass" with your implementation of the desired function


############ 
# Exercise 0
# Read the following utility function and test to get use to python an numpy
# syntax and objects
############

# By default, arrays are float in numpy. The following is a shortener to
# construct integer arrays.
def iarray(x):
	return array(x, dtype=int)

def in_lattice(B, t):
	""" Check wether the vector v belongs to the lattice generated by B,
	up to some very small tolerance to accomodate for numerical errors"""

	# Write v in base B as x: v = x * B
	x = np.linalg.solve(B.transpose(), t)
	# Round to make it an integer vector. If v this indeed a lattice point,
	# this step should merely fix floating-point numerical errors.
	xr = np.round(x)
	# Check equality up to some small tolerance
	# https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
	return np.allclose(x, xr)

############
# Exercise 1
# Implement the Simple Rounding Algorithm 
############

def simple_rounding(B, t):
	""" Returns a lattice vector close to v given a basis B"""
	x = np.linalg.solve(B.transpose(), t)
	xr = np.round(x)
	return xr.dot(B)

if argv[0].endswith("sol1.py") or argv[0].endswith("ex1.py"):
	verify_ex1(simple_rounding)

############
# Exercise 2
# Implement Gram-Schmidt Orthogonalization Algorihm
# Use a projection Subroutine
############

def orth_proj(x, y):
	""" Returns the orthogonal projection of vector x orthognally to vector y"""
	return x - (x.dot(y)/(y.dot(y)))*y

if argv[0].endswith("sol1.py") or argv[0].endswith("ex1.py"):
	verify_ex2a(orth_proj)

def Gram_Schmidt_orth(B):
	""" Output the Gram-Schmidt Orthogonalization of a matrix B"""

	# Get the dimension of the square basis
	n,_ = B.shape
	# Makes a copy of B, but change the type to float
	Bs = array(B, dtype=float)

	for i in range(n):
		for j in range(i):
			Bs[i] = orth_proj(Bs[i], Bs[j])
	return Bs

if argv[0].endswith("sol1.py") or argv[0].endswith("ex1.py"):
	verify_ex2b(Gram_Schmidt_orth)

############
# Exercise 3
# Implement The Nearest-Plane Rounding Algorithm (Babai)
############

def nearest_plane(B, Bs, t):
	""" Given a basis B, its Gram-Schmidt Orthogonalization Bs, 
	and a target t, 	returns a lattice point close to t """

	n,d = B.shape
	# copy t, to avoid modying the value of the function caller
	tt = np.copy(t)

	v = zeros(d, dtype=int)
	for i in reversed(range(n)):
		x = tt.dot(Bs[i]) / Bs[i].dot(Bs[i])
		k = int(round(x))
		tt -= k * B[i]
		v += k * B[i]
	return v
if argv[0].endswith("sol1.py") or argv[0].endswith("ex1.py"):
	verify_ex3(in_lattice, Gram_Schmidt_orth, nearest_plane)

############
# Exercise 4
# Compare the distribution of norm of points in the fundamental
# domains P(B) and P(Bs) given a basis B, by plotting their histograms.
# The plotting function is provided.
#
# There is no need to actually run simple_rounding nor nearest_plane.
############

def plot_two_hist(data_SR, data_NP, n, save=False):
	"""Take is input two lists and plot two histograms"""
	_, bins, _ = plt.hist(data_SR, bins=100, density=True, label="Simple Rounding")
	_ = plt.hist(data_NP, bins=bins, alpha=0.5, density=True, label="Nearest Plane")
	plt.title("Length of random points in Fundamental Parallelepiped \n Basis dimension: %d"%n)
	plt.legend()
	if save:
		plt.savefig("ParallelepipedDistDim%d.png"%n)
	else:
		plt.show()
	plt.clf()
	plt.close()


def compare_norm_distrib(B, samples):
	n, _ = B.shape
	Bs = Gram_Schmidt_orth(B)
	#data_SR = [norm((np.random.rand(n) - array(n*[.5])).dot(B)) for x in range(samples)]
	data_SR = [norm(np.random.uniform(-.5,.5, n).dot(B)) for x in range(samples)]
	data_NP = [norm((np.random.rand(n) - .5).dot(Bs)) for x in range(samples)]
	plot_two_hist(data_SR, data_NP, n)


if argv[0].endswith("sol1.py") or argv[0].endswith("ex1.py"):
	compare_norm_distrib(B2, 50000)
	compare_norm_distrib(B4, 50000)
	compare_norm_distrib(B24, 5000)

