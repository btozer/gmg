"""
Calculate the VGG of a 2D body defined as a polygon of infinite extent orthogonal to the cross section using the
method of Kim and Wessel (2016).

Use the :func:`~fatiando.mesher.Polygon` object to create polygons.

.. warning:: the vertices must be given clockwise! If not, the result will have an inverted sign.

**References**

Kim, S. S., & Wessel, P. (2016). New analytic solutions for modeling vertical gravity gradient anomalies.
Geochemistry, Geophysics, Geosystems, 17(5), 1915–1924. https://doi.org/10.1002/2016GC006263

Uieda, L, Oliveira Jr, V C, Ferreira, A, Santos, H B; Caparica Jr, J F (2014), Fatiando a Terra: a
Python package for modeling and inversion in geophysics. figshare. doi:10.6084/m9.figshare.1115194
"""

import numpy as np
from numpy import arctan2, sin, cos, log

# CONSTANTS
SI_TO_EOTVOS = 1.0e9  # Conversion factor from SI units to EOTVOS :math: `(m/s^2)/m = 1e9`
G = 6.673e-11  # The gravitational constant in :math:`m^3 kg^{-1} s^{-1}`


def gz(xp, zp, polygons):
    """
    Calculates the :math:`g_z` gravity acceleration component.

    .. note:: The coordinate system of the input parameters is z -> **DOWN**.

    .. note:: All input values in **SI** units(!) and output in **mGal**!

    Parameters:

    * xp : array
        The x coordinates of the computation points.

    * zp : array
        The z coordinates of the computation points. (Equals gmg.gravity_observation_elv)

    * polygons : list of :func:`~fatiando.mesher.Polygon`
        Polygons must have the property ``'density'``. Polygons that don't have this property will be ignored in the
        computations. Elements of*polygons* that are None will also be ignored.

        .. note:: The y coordinate of the polygons is used as z!

        .. note:: Data are np arrays

        .. note:: Uses np arrays to calculate gravity effect for slope of node pairs at all Xp observation points
                  simultaneously
    Returns:

    * g_z : array
        The :math:`g_z` component calculated on the computation points
    """

    # INITIALIZE OUTPUT ARRAY
    g_z = np.zeros_like(xp)

    # LOOP THROUGH THE MODEL POLYGONS
    for polygon in polygons:

        # CHECK IF THE CURRENT LAYER HAS A DENSITY CONTRAST SET. SKIP LAYER IF FALSE
        if polygon is None or 'density' not in polygon.props:
            continue
        else:
            density = polygon.props['density']

        # SET X AND Y NODES AND THE NUMBER OF VERTICES IN THE CURRENT LAYER
        x = polygon.x
        z = polygon.y
        nverts = polygon.nverts

        # LOOP THROUGH THE VERTEX PAIRS FOR THE CURRENT POLYGON
        for v in range(nverts):
            # SET THE CURRENT INPUT VERTEX
            xv = x[v] - xp
            zv = z[v] - zp

            # SET THE SECOND VERTEX
            if v == nverts - 1:
                # PAIR THE LAST VERTEX WITH THE FIRST ONE
                xvp1 = x[0] - xp
                zvp1 = z[0] - zp
            else:
                # SET THE SECOND NODE
                xvp1 = x[v + 1] - xp
                zvp1 = z[v + 1] - zp

            del_x = xvp1 - xv
            del_z = zvp1 - zv

            theata_1 = 2 * arctan2(zv, xv)
            theata_2 = 2 * arctan2(zvp1, xvp1)
            del_theata = theata_2 - theata_1

            r1 = xv**2 + zv**2
            r2 = xvp1**2 + zvp1**2
            del_r = r2 - r1

            # RUN ALGORITHM
            current_layer_anomaly = (del_z * (del_x * log(r1/r2) - del_z * del_theata)) / (del_x**2 + del_z**2) \
                                    + sin(theata_2) * log(zvp1) - sin(theata_1) * log(zv)

            # ADD CURRENT LAYER ANOMALY TO TOTAL ANOMALY (VIA SUPER POSITION)
            g_z = g_z + current_layer_anomaly * density * -G * SI_TO_EOTVOS

    # CONVERT TO EOTVOS AND RETURN ARRAY
    return g_z
