"""
Calculate the gravitational attraction of a 2D body with polygonal vertical
cross-section using the formula of Bott et al. (1965)

Use the :func:`~fatiando.mesher.Polygon` object to create polygons.

.. warning:: the vertices must be given clockwise! If not, the result will have
    an inverted sign.

**Components**

* :func:`~fatiando.gravmag.bott.gz`

**References**


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
        If not None, will use this value instead of the ``'density'`` property
        of the polygons. Use this, e.g., for sensitivity matrix building.

        .. note:: The y coordinate of the polygons is used as z!
        .. note:: Data are numpy arrays
        .. note:: Uses numpy arrays to calculate gravity effect for slope of node pairs
                  at all Xp observation points simultaneously
    Returns:

    * gz : array
        The :math:`g_z` component calculated on the computation points

    """

    anomaly = numpy.zeros_like(xp)
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
        for v in xrange(nverts):
            # %Change the coordinates of this vertice
            xv = x[v] - xp
            zv = z[v] - zp
            # %The last vertice pairs with the first one
            if v == nverts - 1:
                xvp1 = x[0] - xp
                zvp1 = z[0] - zp
            else:
                xvp1 = x[v + 1] - xp
                zvp1 = z[v + 1] - zp

            # %Bott ALGORITHM
            theta = -1*arctan2((zvp1 - zv), (xvp1 - xv))
            phi_1 = arctan2(zv, xv)
            phi_2 = arctan2(zvp1, xvp1)
            r1 = numpy.sqrt(zv**2 + xv**2)
            r2 = numpy.sqrt(zvp1**2 + xvp1**2)

            tmp = -(((xv * (sin(theta))) + (zv * (cos(theta)))) * ((sin(theta)) * (log(r2/r1)) +
                    (cos(theta)) * (phi_2-phi_1)) + ((zvp1*phi_2) - (zv*phi_1)))
            anomaly = anomaly + tmp * density
    anomaly = anomaly * SI2MGAL * 2.0 * G
    return anomaly
