"""
Polygon objects - used for passing model layers and attributes to the Bott and Talwani alogrithms

Taken from: Fatiando a Terra: https://github.com/fatiando/fatiando/blob/master/fatiando/mesher/geometry.py

Recreated here because fatiando is a python 2.7 package only and hence incompatible with the gmg python3 project.

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of
the 12th Python in Science Conference, pp. 91 - 98.

"""
from __future__ import division, absolute_import
from future.builtins import object, super
import copy as cp
import numpy as np

# THE SMALLEST SEPARATION, IN METRES, ALLOWED BETWEEN AN OBSERVATION POINT AND A POLYGON VERTEX. EVERY FORWARD
# KERNEL TAKES THE LOGARITHM OF THAT SEPARATION, WHICH IS SINGULAR AT ZERO, SO A COINCIDENT VERTEX IS NUDGED AWAY
SMALL = 1.0e-6


def vertex_pair(x, z, v, nverts, xp, zp):
    """
    Return the v'th polygon side, as vertex coordinates relative to the observation points.

    The last vertex is paired with the first one so that the polygon is closed. A vertex that lands on an
    observation point is nudged clear of it: the forward kernels all evaluate the logarithm of the distance
    between the two, which is singular there and would otherwise put inf or nan into the calculated anomaly. That
    happens in practice whenever a node sits at the calculation elevation and its x coincides with one of the
    observation points, which is easily done since models are routinely drawn from the surface.

    :param x: 1D array (float) : X coordinates of the polygon vertices
    :param z: 1D array (float) : Z coordinates of the polygon vertices
    :param v: int : Index of the first vertex of the side
    :param nverts: int : Number of vertices in the polygon
    :param xp: 1D array (float) : X coordinates of the observation points
    :param zp: float or 1D array (float) : Z coordinates of the observation points

    :returns: tuple of 1D arrays (float) : x1, z1, x2, z2
    """
    x1 = x[v] - xp
    z1 = z[v] - zp

    if v == nverts - 1:
        # PAIR THE LAST VERTEX WITH THE FIRST ONE
        x2 = x[0] - xp
        z2 = z[0] - zp
    else:
        x2 = x[v + 1] - xp
        z2 = z[v + 1] - zp

    on_point = (np.abs(x1) < SMALL) & (np.abs(z1) < SMALL)
    x1 = np.where(on_point, SMALL, x1)
    z1 = np.where(on_point, SMALL, z1)

    on_point = (np.abs(x2) < SMALL) & (np.abs(z2) < SMALL)
    x2 = np.where(on_point, SMALL, x2)
    z2 = np.where(on_point, SMALL, z2)

    return x1, z1, x2, z2


class GeometricElement(object):
    """
    Base class for all geometric elements.
    """

    def __init__(self, props):
        self.props = {}
        if props is not None:
            for p in props:
                self.props[p] = props[p]

    def addprop(self, prop, value):
        """
        Add a physical property to this geometric element.

        If it already has the property, the given value will overwrite the
        existing one.

        Parameters:

        * prop : str
            Name of the physical property.
        * value : float
            The value of this physical property.

        """
        self.props[prop] = value

    def copy(self):
        """ Return a deep copy of the current instance."""
        return cp.deepcopy(self)


class Polygon(GeometricElement):
    """
    A polygon object (2D).

    .. note:: Most applications require the vertices to be **clockwise**!

    Parameters:

    * vertices : list of lists
        List of [x, y] pairs with the coordinates of the vertices.
    * props : dict
        Physical properties assigned to the polygon.
        Ex: ``props={'density':10, 'susceptibility':10000}``

    Examples::

        >>> poly = Polygon([[0, 0], [1, 4], [2, 5]], {'density': 500})
        >>> poly.props
        {'density': 500}
        >>> poly.nverts
        3
        >>> poly.vertices
        array([[0, 0],
               [1, 4],
               [2, 5]])
        >>> poly.x
        array([0, 1, 2])
        >>> poly.y
        array([0, 4, 5])

    """

    def __init__(self, vertices, props=None):
        super().__init__(props)
        self._vertices = np.asarray(vertices)

    @property
    def vertices(self):
        return self._vertices

    @property
    def nverts(self):
        return len(self.vertices)

    @property
    def x(self):
        return self.vertices[:, 0]

    @property
    def y(self):
        return self.vertices[:, 1]