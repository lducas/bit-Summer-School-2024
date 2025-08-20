import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from sol1 import Gram_Schmidt_orth, nearest_plane
from math import sqrt, log
from sys import exit

############
# Exercise 1
# Implement Lagrange reduction Algorithm. The algorithm should output
# the 2x2 Unimodular transformation matrix U. To do so, start with U
# being the identity matrix, and apply the same transformation to U
# than to B all along the algorithm.
############

def lagrange_reduce(B):
    """
    Apply Lagrange reduction to a 2D lattice basis `B` in place.

    :param B: A NumPy array of shape (2, n) representing the lattice basis.
              Each row is a basis vector.
    :type B: numpy.ndarray

    :returns: A 2x2 unimodular integer matrix `U` such that the reduced basis satisfies:
              B_reduced = U @ B_original.
              The input basis `B` is modified in place to its reduced form.
    :rtype: numpy.ndarray

    :raises ValueError: If the input array does not have exactly two row vectors.

    :notes: Make use of in-place swapping of rows, and int() conversion.
    """
    if B.shape[0] != 2:
        raise ValueError("Input basis B must have exactly two vectors (2 rows).")

    pass


# def lagrange_reduce(B):
# 	""" Given a basis with two rows as input, apply lagrange reduction to B (modified 
# 	in place). 
# 	Also output the transformation matrix U, sending the initial basis to the final basis
# 	"""
# 	U = np.identity(2, dtype="int64")
# 	first = True
# 	while first or np.linalg.norm(B[1]) < np.linalg.norm(B[0]):
# 		B[[0, 1]] = B[[1, 0]]
# 		U[[0, 1]] = U[[1, 0]]
# 		first = False
# 		k = int(round(B[0].dot(B[1]) / B[0].dot(B[0])))
# 		B[1] -= k * B[0]
# 		U[1] -= k * U[0]

# 	return U


############
# Exercise 2
# Implement a the Size-Reduction Algorithm (in place) on a basis input B. 
# As a  by-product, provide as output the Gram-Schimdt orthogonalisation of B.
############

def size_reduce(B, Bs):
	"""
    Apply size reduction to a lattice basis `B` using its Gram-Schmidt orthogonalization `Bs`.

    :param B: A NumPy array of shape (n, m), representing the lattice basis.
              Each row is a basis vector. This array is modified in place.
    :type B: numpy.ndarray

    :param Bs: A NumPy array of shape (n, m), representing the Gram-Schmidt orthogonalization
               of `B`. Each row corresponds to the orthogonalized vector of the respective
               row in `B`. This array is not modified.
    :type Bs: numpy.ndarray

    :return: This function modifies `B` in place and does not return a value.
    :rtype: None

    :notes: In the lecture notes, the NearestPlane is applied to a projection
	pi_C'(b_n). This projection is unecessary (but make the proof nicer).
	Ignore it in your implementation.
	The lecture note also give the algorithm in a reccursive manner, but for
	implementation in python, an iterative version will be much simpler
    """

	pass


############
# Exercise 3
# Implement the LLL algorithm. Beware that The Gram-Schmidt basis needs to be 
# updated after every modification of B !
############

def LLL(B, epsilon=0.01, gamma_2=sqrt(4/3), max_iter=1000, animate=True):
	"""
	Perform LLL (Lenstra–Lenstra–Lovász) lattice basis reduction.

	:param B: A NumPy array of shape (n, m), where each row is a basis vector.
	          This matrix is modified in place.
	:type B: numpy.ndarray

	:param epsilon: A small positive float to relax the Lovász condition.
	                Must satisfy 0 < epsilon < 1.
	:type epsilon: float

	:param gamma_2: The delta value in the Lovász condition (typically 0.75).
	:type gamma_2: float

	:param max_iter: Maximum number of LLL iterations before stopping.
	:type max_iter: int

	:param animate: If True, yields the logarithm of the norms of the orthogonal
	                vectors at each iteration.
	:type animate: bool

	:return: None
	:rtype: None

	:notes: Note: at some point, you will need to apply lagrange to the basis 2-dimensional
	basis [pi_i(b_i), pi_i(b_{i+1})]. This can be re-constructed rather cheaply from
	B and B* by noting that:
	pi_i(b_i) = b*_i and 
	pi_i(b_{i+1}) = b*_{i+1} + (<b_i+1, b*_i> / ||b*_i||^2) * b*_i
	"""
	
	pass



# def LLL(B, epsilon=0.01, anim=True):
# 	n,_ = B.shape
# 	Bs = Gram_Schmidt_orth(B)
# 	size_reduce(B, Bs)

# 	while True:
# 		yield [log(np.linalg.norm(x)) for x in Bs]
# 		for i in range(n):
# 			if i==n-1:
# 				return
# 			if np.linalg.norm(Bs[i]) > (gamma_2+epsilon) * np.linalg.norm(Bs[i+1]):
# 				break
# 		# We have found an index i to be Lagrange-reduced
# 		x1 = np.copy(Bs[i])
# 		x2 = Bs[i+1] + (B[i+1].dot(Bs[i]) / (Bs[i].dot(Bs[i]))) * Bs[i]
# 		X = np.array([x1, x2])
# 		U = lagrange_reduce(X)
# 		B[i:i+2] = U.dot(B[i:i+2])
# 		Bs = Gram_Schmidt_orth(B)
# 		size_reduce(B, Bs)

############
# Helper functions
############

def anim_LLL(n, q):
	B = np.identity(n, dtype=int)
	B[0, 0] = q
	for i in range(1, n):
		B[i, 0] = np.random.randint(0, q)

	try:
		data = list(LLL(B, animate=True))
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

	return animation.FuncAnimation(fig=fig,func=update_anim, frames=len(data), interval=50)

############
# Main runner
############
if __name__ == "__main__":
    ani = anim_LLL(n=20, q=999999)
    plt.show()