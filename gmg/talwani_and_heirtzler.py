"""
Calculate the two-dimensional magnetic field strength of a n-sided polygon using the formula from Talwani, M., &
Heirtzler, J. R. (1964). Polygons are of infinite extent orthogonal to the model cross section.

Use the :func:`~fatiando.mesher.Polygon` object to create polygons.

.. warning:: The vertices must be given clockwise! If not, the result will have an inverted sign.

**References**

Talwani, M., & Heirtzler, J. R. (1964). Computation of magnetic anomalies caused by two dimensional structures of
arbitrary shape. In: Computers in the mineral industries, Part 1 (ed.) Parks G A, Stanford Univ. Publ., Geol. Sci.
9 464-480.

Uieda, L, Oliveira Jr, V C, Ferreira, A, Santos, H B; Caparica Jr, J F (2014), Fatiando a Terra: a Python package
for modeling and inversion in geophysics. figshare. doi:10.6084/m9.figshare.1115194
"""

import numpy as np
import math as m

def nt(xp, zp, polygons):
    """
    Calculates the :math:`ntz` magnetic field strength in nT.

    .. note:: The coordinate system of the input parameters is z -> **DOWN**.

    .. note:: All input values in **SI** units and output in **nT**

    Parameters:

    * xp : array
        The x coordinates of the computation points.

    * zp : array
        The z coordinates of the computation points.

    * polygons : list of :func:`~fatiando.mesher.Polygon`
        Polygons must have the property ``'susceptibility'``. Polygons that don't have
        this property will be ignored in the computations. Elements of *polygons* that are None will also be ignored.

    * f: float
        The total magnetic regional field intensity, in nT.

    * angle_a : float
        The vertical magnetisation vector angle, measured in the vertical plane from zero at the horizontal and positive
        downwards, in degrees. Will equal I (the inclination of earth's field) if magnetisation is induced only.

    * angle_b : float
        The angle between the horizontal projection of magnetisation vector and geographic north, measured in the
        horizontal plane, in a positive clockwise direction from geographic north, in degrees. Will equal D (the
        declination of earth's field) if magnetisation is induced only.

    * angle_c : float
        The angle between the positive x axis and geographic north, measured clockwise from geographic north, in
        degrees.

    * observation_elv : float
        The elevation at which to calculate the predeicted magnetic anomaly. Zero is default (For data reduced to MSL).
        Elevations may be positive for areomagnetic surveys.

    .. note:: The y coordinate of the polygons is used as z!

    .. note:: Data are numpy arrays

    .. note:: The algorithm uses numpy arrays to calculate magnetic anomly for each node pairing at all xp observation
              points simultaneously

    Returns:

        * nt : array
            The :math:`n_t` component calculated on the computation points (xp, zp) (TASUM in T&H 1964)
    """

    # INITIALISE TOTAL ANOMALY OUTPUT ARRAYS
    VASUM = np.zeros_like(xp)
    HASUM = np.zeros_like(xp)
    n_t = np.zeros_like(xp)

    # LOOP THROUGH THE MODEL POLYGONS
    for polygon in polygons:

        # CHECK IF THE CURRENT LAYER HAS A SUSCEPTIBILITY CONTRAST SET. SKIP LAYER IF FALSE
        if polygon is None or polygon.props['susceptibility'] == 0.0:
            continue
        else:
            f = polygon.props['f']
            k = polygon.props['susceptibility']
            if f != 1.0:
                # IF MAGNETIZATION IS INDUCED ONLY (f > 1.0) THEN CONVERT TO e.m.u UNITS; ELSE INPUT IS IN (A m^-1)
                k = k / (4 * m.pi)  # CONVERT FROM SI UNITS TO e.m.u (USED IN ORIGINAL CODE)

            # SET EARTH FIELD
            f = polygon.props['f']

            # SET X AND Y NODES AND THE NUMBER OF VERTICES IN THE CURRENT LAYER
            x = polygon.x
            z = polygon.y
            nverts = polygon.nverts

            # DETERMINE ANGLES IN RADIANS
            inclination = polygon.props['angle_a']
            CDIP = m.cos(m.radians(inclination))
            CDIPD = m.cos(m.radians(inclination))
            SDIP = m.sin(m.radians(inclination))
            SDIPD = m.sin(m.radians(inclination))

            declination = polygon.props['angle_b']
            model_azimuth = polygon.props['angle_c']
            SD = m.cos(m.radians(model_azimuth - declination))
            SDD = m.cos(m.radians(model_azimuth - declination))

            # INITIALIZE CURRENT LAYER OUTPUT ARRAYS
            PSUM = np.zeros_like(xp)
            QSUM = np.zeros_like(xp)

            # LOOP THROUGH THE VERTEX PAIRS FOR THE CURRENT POLYGON
            for v in range(0, nverts):
                # SET THE CURRENT INPUT VERTEX
                xv = x[v] - xp
                zv = z[v] - zp

                if v == nverts - 1:
                    # SET THE SECOND VERTEX
                    xv2 = x[0] - xp
                    zv2 = z[0] - zp
                else:
                    # PAIR THE LAST VERTEX WITH THE FIRST ONE
                    xv2 = x[v + 1] - xp
                    zv2 = z[v + 1] - zp

                # RUN T&H ALGORITHM
                R1s = (xv ** 2) + (zv ** 2)
                theta = np.arctan2(zv, xv)
                R2s = (xv2 ** 2) + (zv2 ** 2)
                thetab = np.arctan2(zv2, xv2)
                thetad = theta - thetab
                X1_2 = xv - xv2
                Z2_1 = zv2 - zv
                XSQ = X1_2 ** 2
                ZSQ = Z2_1 ** 2
                XZ = Z2_1 * X1_2
                GL = 0.5 * np.log(R2s / R1s)

                # CALCULATE P AND Q FOR EACH NODE PAIR
                P = ((ZSQ / (XSQ + ZSQ)) * thetad) + ((XZ / (XSQ + ZSQ)) * GL)
                Q = ((XZ / (XSQ + ZSQ)) * thetad) - ((ZSQ / (XSQ + ZSQ)) * GL)
                P[zv == zv2] = 0  # IF zv == zv2 THEN P = 0
                Q[zv == zv2] = 0  # IF zv == zv2 THEN Q = 0

                # ADD P AND Q ANOMALIES TO THE RUNNING TOTAL FOR THE CURRENT POLYGON
                PSUM = P + PSUM
                QSUM = Q + QSUM

            # ADD CURRENT POLYGON ANOMALY TO TOTAL ANOMALY (VIA SUPER POSITION)
            VASUM = VASUM + 2. * k * f * ((CDIP * SD * QSUM) - (SDIP * PSUM))
            HASUM = HASUM + 2. * k * f * ((CDIP * SD * PSUM) + (SDIP * QSUM))

        # CALCULATE TOTAL FIELD ANOMALY
        n_t = (HASUM * CDIPD * SDD) + (VASUM * SDIPD)
    return n_t
