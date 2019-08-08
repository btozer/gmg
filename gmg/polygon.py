"""
Polygon objects - used for passing model layers and attributes to the Bott and Talwani alogrithms

Taken from: Fatiando a Terra: https://github.com/fatiando/fatiando/blob/master/fatiando/mesher/geometry.py

Recreated here becuase fatiando is a python 2.7 package only and hence incompatiable with the gmg python3 project.

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of
the 12th Python in Science Conference, pp. 91 - 98.

"""
from __future__ import division, absolute_import
from future.builtins import object, super
import copy as cp
import numpy as np


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