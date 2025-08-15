import numpy as np

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

