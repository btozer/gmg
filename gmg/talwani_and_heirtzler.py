"""
Calculate the two-dimensional magnetic field strength of a n-sided polygon using the formula from Talwani, M., &
Heirtzler, J. R. (1964). Polygon is infinte orthogonal to the model cross section.

Use the :func:`~fatiando.mesher.Polygon` object to create polygons.

.. warning:: the vertices must be given clockwise! If not, the result will have
    an inverted sign.

**Components**

* :func:`~gmg.gravmag.talwani_and_heirtzler.ntz`

**References**

Talwani, M., & Heirtzler, J. R. (1964). Computation of magnetic anomalies caused by two dimensional structures of
arbitrary shape. In: Computers in the mineral industries, Part 1 (ed.) Parks G A, Stanford Univ. Publ., Geol. Sci.
9 464-480.

Uieda, L, Oliveira Jr, V C, Ferreira, A, Santos, H B; Caparica Jr, J F (2014), Fatiando a Terra: a Python package 
for modeling and inversion in geophysics. figshare. doi:10.6084/m9.figshare.1115194
----

"""

import numpy as np
import math as m

def ntT(xp, zp, polygons, profile_azimuth, declination, inclination, F, remanence, angle_a, angle_b):

    """
    Calculates the :math:`nt_T` magnetic field strength.

    .. note:: The coordinate system of the input parameters is z -> **DOWN**.

    .. note:: All input values in **SI** units(!) and output in **nT**!

    Parameters:

    * xp, zp : arrays
        The x and z coordinates of the computation points.
    * polygons : list of :func:`~fatiando.mesher.Polygon`
        The susceptibility model used.
        Polygons must have the property ``'susceptibility'``. Polygons that don't have
        this property will be ignored in the computations. Elements of
        *polygons* that are None will also be ignored.
    * susceptibility : float or None
    * angle_a : float or None
    * angle_b : float or None

        .. note:: The y coordinate of the polygons is used as z!
        .. note:: Data are numpy arrays
        .. note:: Uses numpy arrays to calculate magnetic effect for slope of node pairs
                  at all Xp observation points simultaneously
    Returns:

    * nt_T : array
        The :math:`nt_T` component calculated on the computation points

    """

    # %INITIALISE ARRAYS
    VASUM = np.zeros_like(xp)
    HASUM = np.zeros_like(xp)
    TASUM = np.zeros_like(xp)

    # %SET POLYGON NUMBER
    i = 0
    for polygon in polygons:

        if polygon is None or 'susceptibility' not in polygon.props:
            continue
        else:
            C = (polygon.props['susceptibility'])  # %CONVERT FROM SI INPUT TO e.m.u (used in original code)
        if C == 0.0:  # %IF C = 0 THEN SKIP THIS POLYGON; ELSE RUN THE CODE
            continue
        else:
            # %DETERMINE ANGLES IN RADIANS
            print remanence[i]
            if remanence[i] == 1:
                inclination = angle_a[i]
                declination = angle_b[i]

            CDIP = m.cos(m.radians(inclination))
            CDIPD = m.cos(m.radians(inclination))
            SDIP = m.sin(m.radians(inclination))
            SDIPD = m.sin(m.radians(inclination))

            SD = m.cos(m.radians(profile_azimuth-declination))
            SDD = m.cos(m.radians(profile_azimuth-declination))
            x = polygon.x
            z = polygon.y
            nverts = polygon.nverts

            PSUM  = np.zeros_like(xp)
            QSUM  = np.zeros_like(xp)

            # %SET VERTICE
            for v in range(0, nverts):
                xv = x[v] - xp
                zv = z[v] - zp
                # %SET NEXT VERT: THE LAST VERTICE PAIRS WITH THE FIRST ONE
                if v == nverts - 1:
                    xv2 = x[0] - xp
                    zv2 = z[0] - zp
                else:
                    xv2 = x[v + 1] - xp
                    zv2 = z[v + 1] - zp

                # %RUN ALGORITHM
                R1s = (xv**2)+(zv**2)
                theta = np.arctan2(zv, xv)
                R2s = (xv2**2)+(zv2**2)
                thetab = np.arctan2(zv2, xv2)
                thetad = theta-thetab
                X1_2 = xv-xv2
                Z2_1 = zv2-zv
                XSQ  = X1_2**2
                ZSQ  = Z2_1**2
                XZ   = Z2_1 * X1_2
                GL   =0.5 * np.log(R2s/R1s)

                P    = ((ZSQ/(XSQ+ZSQ))*thetad)+((XZ/(XSQ+ZSQ))*GL)
                Q    = ((XZ/(XSQ+ZSQ))*thetad)-((ZSQ/(XSQ+ZSQ))*GL)
                P[zv == zv2] = 0
                Q[zv == zv2] = 0

                PSUM = P + PSUM
                QSUM = Q + QSUM

            VASUM = VASUM + 2.*C*F*((CDIP*SD*QSUM)-(SDIP*PSUM))
            HASUM = HASUM + 2.*C*F*((CDIP*SD*PSUM)+(SDIP*QSUM))
        TASUM = (HASUM*CDIPD*SDD)+(VASUM*SDIPD)
        i += 1
    return TASUM
