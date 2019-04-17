"""
Calculate the gravitational attraction of a 2D body defined as a polygon of
infinte extent orthogonal to the cross section using the method of Bott (1965).

Use the :func:`~fatiando.mesher.Polygon` object to create polygons.

.. warning:: the vertices must be given clockwise! If not, the result will have
    an inverted sign.

**Components**

* :func:`~gmg.gravmag.bott.gz`

**References**

Bott, M. H. P. (1969). GRAVN. Durham geophysical computer specification No. 1.

Uieda, L, Oliveira Jr, V C, Ferreira, A, Santos, H B; Caparica Jr, J F (2014), Fatiando a Terra: a
Python package for modeling and inversion in geophysics. figshare. doi:10.6084/m9.figshare.1115194
"""

import numpy
from numpy import arctan2, sin, cos, log
from fatiando.constants import G, SI2MGAL


def Gz(xp, zp, polygons, dens=None):
    """
    Calculates the :math:`g_z` gravity acceleration component.

    .. note:: The coordinate system of the input parameters is z -> **DOWN**.

    .. note:: All input values in **SI** units(!) and output in **mGal**!

    Parameters:

    * xp, zp : arrays
        The x and z coordinates of the computation points.
    * polygons : list of :func:`~fatiando.mesher.Polygon`
        The density model used.
        Polygons must have the property ``'density'``. Polygons that don't have
        this property will be ignored in the computations. Elements of
        *polygons* that are None will also be ignored.
    * dens : float or None

        .. note:: The y coordinate of the polygons is used as z!
        .. note:: Data are numpy arrays
        .. note:: Uses numpy arrays to calculate gravity effect for slope of node pairs
                  at all Xp observation points simultaneously
    Returns:

    * g_z : array
        The :math:`g_z` component calculated on the computation points

    """

    g_z = numpy.zeros_like(xp)
    for polygon in polygons:
        if polygon is None or ('density' not in polygon.props and dens is None):
            continue
        if dens is None:
            density = polygon.props['density']
        else:
            density = dens
        x = polygon.x
        z = polygon.y
        nverts = polygon.nverts
        for v in range(nverts):
            # CHANGE THE COORDINATES OF THIS VERTICE
            xv = x[v] - xp
            zv = z[v] - zp

            # THE LAST VERTICE PAIRS WITH THE FIRST ONE
            if v == nverts - 1:
                xvp1 = x[0] - xp
                zvp1 = z[0] - zp
            else:
                xvp1 = x[v + 1] - xp
                zvp1 = z[v + 1] - zp

            # RUN BOTT ALGORITHM
            theta = -1 * arctan2((zvp1 - zv), (xvp1 - xv))
            phi_1 = arctan2(zv, xv)
            phi_2 = arctan2(zvp1, xvp1)
            r1 = numpy.sqrt(zv ** 2 + xv ** 2)
            r2 = numpy.sqrt(zvp1 ** 2 + xvp1 ** 2)

            tmp = (((xv * (sin(theta))) + (zv * (cos(theta)))) * ((sin(theta)) * (log(r2 / r1)) +
                                                                  (cos(theta)) * (phi_2 - phi_1)) +
                                                                    ((zvp1 * phi_2) - (zv * phi_1)))

            g_z = g_z + tmp * density
    g_z = g_z * SI2MGAL * 2.0 * G
    return g_z
