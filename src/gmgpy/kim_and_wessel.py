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
from numpy import arctan2, sin, log

# gmg.py PUTS ITS OWN DIRECTORY ON sys.path AND IMPORTS ITS MODULES FLAT, WHILE THE PACKAGED TESTS IMPORT gmgpy AS
# A PACKAGE, SO THE SHARED HELPER HAS TO BE REACHABLE BOTH WAYS
try:
    from .polygon import vertex_pair
except ImportError:
    from polygon import vertex_pair

# CONSTANTS
SI_TO_EOTVOS = 1.0e9  # Conversion factor from SI units to EOTVOS :math: `(m/s^2)/m = 1e9`
G = 6.673e-11  # The gravitational constant in :math:`m^3 kg^{-1} s^{-1}`


def _sin_log_z(theata, zv):
    """
    Evaluate the :math:`\\sin(2\\theta) \\log z` term of Kim and Wessel (2016) for any sign of z, including zero.

    Two things have to be handled that the bare expression does not.

    The **absolute value** is exact rather than a guard against a nan. The term enters as
    :math:`\\sin\\theta_{2} \\log z_{2} - \\sin\\theta_{1} \\log z_{1}` summed around a *closed* polygon; the
    analytic continuation to :math:`z < 0` is :math:`\\log|z| + i\\pi`, and the :math:`i\\pi` part telescopes to
    zero around the loop. By the same argument the arbitrary scale of the logarithm cancels, which is what makes
    taking the logarithm of a bare coordinate admissible here at all. Without it, a node at or above the
    calculation elevation returns nan and, because the polygons are summed, blanks the whole profile.

    :math:`z = 0` is a **removable** singularity: :math:`\\sin 2\\theta = 2xz/r^{2}` vanishes linearly in z while
    :math:`\\log|z|` diverges only logarithmically, so the product tends to zero. A node drawn exactly at the
    calculation elevation hits this at every observation point at once, not merely at the one it touches.

    :param theata: 1D array (float) : Twice the angle to the vertex, in radians
    :param zv: 1D array (float) : Depth of the vertex relative to the observation points

    :returns: 1D array (float) : The value of the term
    """
    # SUBSTITUTING 1 FOR A ZERO DEPTH BOTH KEEPS THE LOGARITHM AWAY FROM ZERO, SO NO DIVIDE BY ZERO IS EVER
    # EVALUATED, AND MAKES THE TERM log(1) = 0, WHICH IS EXACTLY THE LIMIT DESCRIBED ABOVE
    return sin(theata) * log(np.abs(np.where(zv == 0.0, 1.0, zv)))


def gz(xp, zp, polygons):
    """
    Calculates the :math:`g_z` gravity acceleration component.

    .. note:: The coordinate system of the input parameters is z -> **DOWN**.

    .. note:: All input values in **SI** units(!) and output in **mGal**!

    Parameters:

    * xp : array
        The x coordinates of the computation points.

    * zp : float
       The z elevation of the computation points (Equals gmg.vgg_observation_elv)
    
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
            # SET THE VERTEX PAIR, RELATIVE TO THE OBSERVATION POINTS AND CLEAR OF ANY OF THEM
            xv, zv, xvp1, zvp1 = vertex_pair(x, z, v, nverts, xp, zp)

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
                                    + _sin_log_z(theata_2, zvp1) - _sin_log_z(theata_1, zv)

            # ADD CURRENT LAYER ANOMALY TO TOTAL ANOMALY (VIA SUPER POSITION)
            g_z = g_z + current_layer_anomaly * density * -G * SI_TO_EOTVOS

    # RETURN ANOMALY IN EOTVOS
    return g_z
