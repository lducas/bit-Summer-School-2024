import numpy as np
import matplotlib.pyplot as plt
from sol1 import simple_rounding
from time import sleep
# from itertools import product

# The exercises comprises of function to be implemented (except the first function,
# namely, generate_lattice_points that is already implemented). For the rest of the
# functions, replace the keyword "pass" with your implementation.

def generate_lattice_points(B, xlim, ylim):
    """
    Generate all lattice points lying inside a given axis-aligned rectangular
    region in Euclidean space, for a 2D lattice defined by a basis matrix.

    Given a 2×2 lattice basis `B`, each lattice point can be expressed as
    B @ [i, j] where i and j are integers (lattice coordinates). The function
    returns all such points whose Euclidean coordinates (x, y) lie within
    the rectangle defined by `xlim` and `ylim`.

    :param B: A 2×2 NumPy array whose columns are the basis vectors of the lattice.
    :type B: numpy.ndarray
    :param xlim: Tuple (xmin, xmax) specifying the horizontal extent of the box
                in Euclidean coordinates.
    :type xlim: tuple of float
    :param ylim: Tuple (ymin, ymax) specifying the vertical extent of the box
                in Euclidean coordinates.
    :type ylim: tuple of float

    :return: A NumPy array of shape (m, 2), where each row is the (x, y) coordinate
            of a lattice point inside the specified box.
    :rtype: numpy.ndarray

    :notes: This function works by mapping the bounding box corners into lattice
            coordinate space using the inverse of `B`, determining the integer
            index range that covers the box, and then filtering only those lattice
            points that fall inside the Euclidean bounds.
    """
    lattice_points = []
    
    # The corners of the bounding box in Euclidean space
    corners = [
        [xlim[0], ylim[0]],
        [xlim[0], ylim[1]],
        [xlim[1], ylim[0]],
        [xlim[1], ylim[1]],
    ]

    # Map the corners to lattice coordinates
    coords_corners = [np.linalg.solve(B.transpose(), c) for c in corners]

    # Determine min/max in lattice coordinates and expand to integer grid
    min_i = int(np.floor(min(c[0] for c in coords_corners)))
    max_i = int(np.ceil(max(c[0] for c in coords_corners)))
    min_j = int(np.floor(min(c[1] for c in coords_corners)))
    max_j = int(np.ceil(max(c[1] for c in coords_corners)))

    # Generate points and filter to actual Euclidean box
    for i in range(min_i, max_i + 1):
        for j in range(min_j, max_j + 1):
            p = np.array([i, j]) @ B
            if xlim[0] <= p[0] <= xlim[1] and ylim[0] <= p[1] <= ylim[1]:
                lattice_points.append(p)

    return np.array(lattice_points)



############
# Exercise 1
# Enumerate lattice vectors within a ball of radius r
############

def enumerate(x, r):
    """
    Return all integer vectors in Z^n whose coordinates differ from `x`
    by at most `r` (per coordinate).

    :param t: A NumPy array representing the vector `x`.
    :type t: numpy.ndarray
    :param l: An integer or float representing the distance `r`.
    :type l: int or float

    :return: A list of lattice vectors in Z^n whose coordinates are
    within distance `r` (per coordinate) from `x`.
    :rtype: list of numpy.ndarray

    :notes: Make use of numpy concatenate() function and the built-in append function.
    """

    # Base case: no dimensions left
    if len(x) == 0:
        return [np.array([], dtype=int)]
    
    results = []

    # Loop over all integer coordinates within ±r of t[0]
    lower = int(np.ceil(x[0] - r))
    upper = int(np.ceil(x[0] + r))

    for coord in range(lower, upper):
        # Recursively enumerate for the rest of the coordinates
        for sub_result in enumerate(x[1:], r):
            # Create a new vector by prepending the current coordinate
            vec = np.concatenate(([coord], sub_result))
            results.append(vec)

    return results

# def enumerate_Zn_ball(x, r):
#     """
#     Return all lattice vectors in Z^n that whose coordinates at the distance at most r
#     from the corresponding cooridnates of the target vector `t`.

#     :param t: A NumPy array representing the target vector.
#     :type t: numpy.ndarray
#     :param l: An integer or float representing the distance.
#     :type l: int or float

#     :return: A list of lattice vectors in Z^n whose coordinates are within distance 'r' from `t`.
#     :rtype: list of numpy.ndarray

#     """

#     # Search around the nearest integer vector
#     results = []
#     for offsets in product(range(-int(np.ceil(r)), int(np.ceil(r))), repeat=len(x)):
#         # Create a candidate vector by adding offsets to the target vector
#         y = x + np.array(offsets)

#         # Check if the candidate vector is in integer lattice
#         if np.allclose(y, np.round(y), atol=1e-9):
#             results.append(y)

#     return results

############
# Exercise 2
# Implement the Simple Enumeration Algorithm 
############

def simple_enumeration(B, t, l):
    """
	Return a lattice vector close to the target vector `t`, using the Simple Enumeration algorithm.

	:param B: A square (n x n) NumPy array representing a lattice basis.
	:type B: numpy.ndarray
	:param t: A NumPy array representing the target vector to approximate with a lattice point.
	:type t: numpy.ndarray
	:param l: An integer value representing the fundamental domain scaling.
	:type l: int

	:return: A lattice vector in the lattice generated by `B` that is close to `t`.
	:rtype: numpy.ndarray

	:notes: Make use of numpy.linalg function solve and numpy function round.
	"""

    # Initialize the best lattice vector found so far
    c = np.ndarray(t.shape, dtype=int)
    c.fill(np.iinfo(int).max)
	
    # Enumerate integer coordinate vectors v inside an l/2 ball around x
    # Note: make us of enumeration_Zn_ball() that was pereviously implemented
    x = np.linalg.solve(B.transpose(), t)
    enum = enumerate(x, l / 2)

    assert enum, "Enumeration returned None, expected a list of vectors."
	
    # For each v in the ball, project it into the lattice and check if the
    # corresponding lattice vetor is closer to t than the current best candidate
    for v in enum:
        y = v @ B
        if np.linalg.norm(y - t) < np.linalg.norm(c - t):
            c = y
	
	# Step 5: Return the closest lattice vector found
    return c


if __name__ == "__main__":

    # Slightly skewed (non-orthogonal) lattice basis
    B = np.array([[1.5, 0.9], [0.7, 1.5]], dtype=float)

    # Example target vector
    t = np.array([2, 1], dtype=float)

    # Example parameter l
    l = 4

    # Enumerating all integer vectors around target coordinates
    print("Target vector:", t)
    print("\nEnumerating integer vectors around target coordinates...")
    all_vecs = np.array(enumerate(t, l / 2))


    # Running simple rounding algorithm
    print("\nRunning simple rounding...")
    rounding_vec = simple_rounding(B, t)
    print("Closest lattice vector found using simple rounding:", rounding_vec)

    # Running simple enumeration algorithm
    print("\nRunning simple enumeration...")
    simple_enum_vec = simple_enumeration(B, t, l)
    print("Closest lattice vector found using simple enumeration:", simple_enum_vec)

    # Generate full lattice for a visible region
    lattice_vecs = generate_lattice_points(B, (-5, 5), (-5, 5))

    # Plot 1 — enumeration region within lattice
    plt.scatter(lattice_vecs[:, 0], lattice_vecs[:, 1], c="gray", alpha=0.4, label="Lattice Points")
    plt.scatter(all_vecs[:, 0], all_vecs[:, 1], c="lightblue", marker="s", label="Enumerated Points")
    plt.scatter(t[0], t[1], c="red", marker="x", s=50, label="Target")
    plt.legend()
    plt.title("Enumerated Region in Skewed Lattice")
    plt.show()

    sleep(2)  # Pause before next plot

    # Plot 2 — with results highlighted
    plt.scatter(lattice_vecs[:, 0], lattice_vecs[:, 1], c="gray", alpha=0.4, label="Lattice Points")
    plt.scatter(all_vecs[:, 0], all_vecs[:, 1], c="lightblue", marker="s", label="Enumerated Points")
    plt.scatter(t[0], t[1], c="red", marker="x", s=50, label="Target")
    plt.scatter(rounding_vec[0], rounding_vec[1], c="orange", label="Rounding Result")
    plt.scatter(simple_enum_vec[0], simple_enum_vec[1], c="green", label="Enumeration Result")
    plt.legend()
    plt.title("Comparison: Rounding vs Enumeration in Skewed Lattice")
    plt.show()

