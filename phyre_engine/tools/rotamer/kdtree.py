import numpy as np
from scipy.spatial import cKDTree
"""
Module containing a periodic wrapper of :py:class:`scipy.spatial.cKDTree` for
quick lookup of residues adopting similar backbone angles.
"""

class PeriodicKDTree:
    """
    Periodic wrapper around :py:class:`scipy.spatial.cKDTree`.

    The purpose of this class is to facilitate quick lookup of residues adopting
    similar backbone angles, but it can be used for any periodic k-nearest
    neighbours lookup.

    Scipy provides an implementation of a KDTree for looking up nearest
    neighbours, but it cannot account for periodic images. When we are working
    with residue dihedral angles, we want to be able to find the residues with
    the closest angles, but need to account for the periodic nature of the
    torsion angles.

    This class is not particularly intelligent: when looking up a point in the
    quadtree, the point is evaluated over nine positions. For an angle in
    degrees, we find the nearest neighbours of point :math:`x` at :math:`x -
    360, x, x + 360`. In two dimensions, this requires nine evaluations.
    """

    def __init__(self, coord_list, width, height):
        """
        Initialise a periodic quadtree containing each coordinate in
        ``coord_list``.

        :param coord_list: Numpy array of 2D points.
        :param width: Width of the system about which to loop.
        :param height: Height of the system.
        """
        self.kdtree = cKDTree(coord_list)
        self.width = width
        self.height = height

    def query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf):
        """
        Find the ``k`` points closest to ``x``.

        Parameters and return values are defined identically to
        :py:meth:`scipy.spatial.cKDTree.query`.

        The closest points are found by finding the closest points in each
        periodic image, then picking the closest ``k`` unique points.
        """

        distances = []
        indices = []
        for x_img in self.periodic_images(x):
            img_distances, img_indices = self.kdtree.query(
                x_img, k, eps, p, distance_upper_bound)
            distances.append(img_distances)
            indices.append(img_indices)

        distances = np.concatenate(tuple(distances))
        indices = np.concatenate(tuple(indices))

        order = distances.argsort()
        distances = distances[order]
        indices = indices[order]

        unique_indices = np.sort(np.unique(indices, return_index=True)[1])
        indices = indices[unique_indices]
        distances = distances[unique_indices]

        num_to_return = min(k, len(unique_indices))
        return distances[0:num_to_return], indices[0:num_to_return]

    def periodic_images(self, point):
        """Yield a numpy array containing periodic images of ``point``."""
        for i in range(-1, 2):
            for j in range(-1, 2):
                x = point[0] + (i * self.width)
                y = point[1] + (j * self.height)
                yield np.array([x, y])
