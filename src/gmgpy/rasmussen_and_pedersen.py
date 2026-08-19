"""
Calculate the gravity and magnetic anomalies of a 2.75D body, that is, a body defined as a polygon in the model cross
section which has a *finite* extent in the strike direction (orthogonal to the section). The extent is given by the two
signed distances y1 and y2, measured from the plane of the section. 2.5D is the symmetric special case y1 = -y2.

The formulation is that of Rasmussen and Pedersen (1979), who showed that the well known two dimensional formulae of
Talwani et al. (1959) and Talwani and Heirtzler (1964) require only a modest change to represent a body whose tails in
the strike direction have been cut off. A single fundamental integral serves both the gravity and the magnetic case.

As y1 -> -infinity and y2 -> +infinity the expressions reduce to the two dimensional result. Rather than rely on that
limit numerically, a polygon whose ``'y1'`` or ``'y2'`` property is ``None`` is taken to be of infinite extent and is
passed straight to the existing two dimensional kernels (:mod:`~gmgpy.bott` and
:mod:`~gmgpy.talwani_and_heirtzler`), whose results are then summed with the 2.75D ones by superposition. A model may
therefore mix bodies of finite and infinite extent, and a model in which every body is of infinite extent reproduces
the two dimensional result exactly.

.. warning:: The vertices must be given in the same order as for :mod:`~gmgpy.bott` and
             :mod:`~gmgpy.talwani_and_heirtzler`. If not, the result will have an inverted sign.

.. note:: There is no 2.75D vertical gravity gradient here. Rasmussen and Pedersen give the gravity and the magnetic
          field but not the vertical gradient, so gmg calculates the VGG in two dimensions only
          (:mod:`~gmgpy.kim_and_wessel`). In a model containing bodies of finite extent the gravity and magnetic
          anomalies will honour that extent but the VGG will not.

**References**

Rasmussen, R., & Pedersen, L. B. (1979). End corrections in potential field modeling. Geophysical Prospecting, 27(4),
749-760. doi:10.1111/j.1365-2478.1979.tb00994.x

Zunino, A., Ghirotto, A., Armadillo, E., & Fichtner, A. (2022). Hamiltonian Monte Carlo probabilistic joint inversion
of 2D (2.75D) gravity and magnetic data. Geophysical Research Letters, 49, e2022GL099789. doi:10.1029/2022GL099789

Talwani, M., & Heirtzler, J. R. (1964). Computation of magnetic anomalies caused by two dimensional structures of
arbitrary shape. In: Computers in the mineral industries, Part 1 (ed.) Parks G A, Stanford Univ. Publ., Geol. Sci.
9 464-480.
"""

import numpy as np

# gmg.py PUTS ITS OWN DIRECTORY ON sys.path AND IMPORTS ITS MODULES FLAT, WHILE THE PACKAGED TESTS IMPORT gmgpy AS A
# PACKAGE. THIS IS THE FIRST MODULE HERE THAT NEEDS A SIBLING, SO IT HAS TO WORK BOTH WAYS.
try:
    from . import bott
    from . import talwani_and_heirtzler
except ImportError:
    import bott
    import talwani_and_heirtzler

G = bott.G
SI2MGAL = bott.SI2MGAL

# CONSTANTS
MU0_NT = 400.0 * np.pi  # PERMEABILITY OF FREE SPACE IN nT m A^-1 (I.E. mu_0 SCALED FROM T TO nT)

# THE SMALLEST SEPARATION, IN METRES, ALLOWED BETWEEN AN OBSERVATION POINT AND A POLYGON VERTEX. A VERTEX CLOSER THAN
# THIS IS NUDGED AWAY, SINCE THE LOGARITHMS BELOW ARE SINGULAR THERE. THE TWO DIMENSIONAL KERNELS HAVE NO SUCH GUARD
# AND RETURN inf/nan INSTEAD.
SMALL = 1.0e-6


def gz(xp, zp, polygons):
    """
    Calculates the :math:`g_z` gravity acceleration component of a set of 2.75D bodies.

    .. note:: The coordinate system of the input parameters is z -> **DOWN**.

    .. note:: All input values in **SI** units(!) and output in **mGal**!

    Parameters:

    * xp : array
        The x coordinates of the computation points.

    * zp : float or array
        The z elevation of the computation points (Equals gmg.gravity_observation_elv).

    * polygons : list of :func:`~gmgpy.polygon.Polygon`
        Polygons must have the property ``'density'``. Polygons that don't have this property will be ignored in the
        computations. Elements of *polygons* that are None will also be ignored.

        Polygons may also carry the properties ``'y1'`` and ``'y2'``, the signed extents of the body in the strike
        direction, in metres. If either is absent or None the body is treated as being of infinite extent and is
        passed to :func:`~gmgpy.bott.gz`.

        .. note:: The y coordinate of the polygons is used as z!

    Returns:

    * g_z : array
        The :math:`g_z` component calculated on the computation points
    """
    infinite, finite = _split_by_extent(polygons)

    # BODIES OF INFINITE EXTENT ARE THE EXISTING TWO DIMENSIONAL PROBLEM
    g_z = bott.gz(xp, zp, infinite)

    # LOOP THROUGH THE BODIES OF FINITE EXTENT
    for polygon in finite:
        if 'density' not in polygon.props:
            continue
        density = polygon.props['density']
        y1, y2 = _extents(polygon)

        # SET X AND Y NODES AND THE NUMBER OF VERTICES IN THE CURRENT LAYER
        x = polygon.x
        z = polygon.y
        nverts = polygon.nverts

        current_layer_anomaly = np.zeros_like(xp)

        # LOOP THROUGH THE VERTEX PAIRS FOR THE CURRENT POLYGON
        for v in range(nverts):
            x1, z1, x2, z2 = _vertex_pair(x, z, v, nverts, xp, zp)

            cos_phi, sin_phi, s, u1, u2, w = _rotated_frame(x1, z1, x2, z2)
            r1_y1, r1_y2, r2_y1, r2_y2, lor, difatan = _log_and_angle_terms(u1, u2, w, y1, y2)

            # log((u2 + R2(y)) / (u1 + R1(y))) EVALUATED AT EACH END OF THE BODY
            lor2_y1 = np.log((u2 + r2_y1) / (u1 + r1_y1))
            lor2_y2 = np.log((u2 + r2_y2) / (u1 + r1_y2))

            current_layer_anomaly = current_layer_anomaly + \
                (cos_phi * ((y2 * lor2_y2) - (y1 * lor2_y1)) +
                 (((x1 * z2) - (z1 * x2)) / s) * ((cos_phi * difatan) - (sin_phi * lor)))

        # ADD CURRENT LAYER ANOMALY TO TOTAL ANOMALY (VIA SUPER POSITION), CONVERTING TO MGAL
        # NB: THE SIGN HERE IS SET TO MATCH THE VERTEX ORDERING CONVENTION OF bott.gz
        g_z = g_z + (current_layer_anomaly * density * G * SI2MGAL)

    return g_z


def nt(xp, zp, polygons):
    """
    Calculates the :math:`n_t` total field magnetic anomaly of a set of 2.75D bodies, in nT.

    .. note:: The coordinate system of the input parameters is z -> **DOWN**.

    .. note:: All input values in **SI** units and output in **nT**

    Parameters:

    * xp : array
        The x coordinates of the computation points.

    * zp : float or array
        The z elevation of the computation points (Equals gmg.mag_observation_elv)

    * polygons : list of :func:`~gmgpy.polygon.Polygon`
        Polygons must have the properties ``'susceptibility'``, ``'f'``, ``'angle_a'``, ``'angle_b'`` and
        ``'angle_c'``, as documented in :func:`~gmgpy.talwani_and_heirtzler.nt`. Polygons that don't have these
        properties will be ignored in the computations, as will elements of *polygons* that are None.

        Polygons may also carry the properties ``'y1'`` and ``'y2'``, the signed extents of the body in the strike
        direction, in metres. If either is absent or None the body is treated as being of infinite extent and is
        passed to :func:`~gmgpy.talwani_and_heirtzler.nt`.

    .. note:: Magnetisation is taken to be induced only, matching :mod:`~gmgpy.talwani_and_heirtzler`: the angles A
              and B describe the magnetisation vector *and* set the direction the anomaly is projected onto. The
              along strike component of magnetisation, :math:`J_{y}`, contributes nothing in two dimensions but does
              once the body has ends, so it is carried here.

    Returns:

    * n_t : array
        The total field anomaly calculated on the computation points (xp, zp)
    """
    infinite, finite = _split_by_extent(polygons)

    # BODIES OF INFINITE EXTENT ARE THE EXISTING TWO DIMENSIONAL PROBLEM
    n_t = talwani_and_heirtzler.nt(xp, zp, infinite)

    # LOOP THROUGH THE BODIES OF FINITE EXTENT
    for polygon in finite:
        if polygon.props.get('susceptibility', 0.0) == 0.0:
            continue  # SKIP THIS POLYGON
        y1, y2 = _extents(polygon)

        # INDUCED MAGNETISATION CONTRAST, IN A/m
        k = polygon.props['susceptibility']
        f = polygon.props['f']  # EARTHS FIELD (nT)
        j = k * f / MU0_NT

        # DETERMINE ANGLES IN RADIANS
        angle_a = np.radians(polygon.props['angle_a'])  # MAGNETISATION VECTOR DIP
        angle_b = np.radians(polygon.props['angle_b'])  # MAGNETISATION VECTOR AZIMUTH
        angle_c = np.radians(polygon.props['angle_c'])  # MODEL PROFILE AZIMUTH

        # RESOLVE THE MAGNETISATION VECTOR ONTO THE MODEL FRAME. NB: J_y CONTRIBUTES NOTHING IN TWO DIMENSIONS, WHICH
        # IS WHY talwani_and_heirtzler CARRIES ONLY J_x AND J_z, BUT IT DOES CONTRIBUTE ONCE THE BODY HAS ENDS.
        # NB: THE MAGNETIC EXPRESSIONS BELOW MEASURE y IN THE OPPOSITE SENSE TO y1/y2 AND TO THE GRAVITY KERNEL, SO
        # THE ALONG STRIKE COMPONENT IS NEGATED HERE ON THE WAY IN AND AGAIN ON THE WAY OUT WHEN b_y IS PROJECTED.
        # THIS ONLY SHOWS UP FOR AN ASYMMETRIC BODY (y1 != -y2); SEE test_04_rasmussen_and_pedersen.py.
        j_x = j * np.cos(angle_a) * np.cos(angle_c - angle_b)
        j_y = -j * np.cos(angle_a) * np.sin(angle_c - angle_b)
        j_z = j * np.sin(angle_a)

        # SET X AND Y NODES AND THE NUMBER OF VERTICES IN THE CURRENT LAYER
        x = polygon.x
        z = polygon.y
        nverts = polygon.nverts

        # INITIALIZE CURRENT LAYER FIELD COMPONENT ARRAYS
        b_x = np.zeros_like(xp)
        b_y = np.zeros_like(xp)
        b_z = np.zeros_like(xp)

        # LOOP THROUGH THE VERTEX PAIRS FOR THE CURRENT POLYGON
        for v in range(nverts):
            x1, z1, x2, z2 = _vertex_pair(x, z, v, nverts, xp, zp)

            cos_phi, sin_phi, s, u1, u2, w = _rotated_frame(x1, z1, x2, z2)
            r1_y1, r1_y2, r2_y1, r2_y2, lor, difatan = _log_and_angle_terms(u1, u2, w, y1, y2)

            # log(((R2(y) - u2)(R1(y) + u1)) / ((R2(y) + u2)(R1(y) - u1))) EVALUATED AT EACH END OF THE BODY
            lor2_y1 = np.log(((r2_y1 - u2) * (r1_y1 + u1)) / ((r2_y1 + u2) * (r1_y1 - u1)))
            lor2_y2 = np.log(((r2_y2 - u2) * (r1_y2 + u1)) / ((r2_y2 + u2) * (r1_y2 - u1)))

            term_1 = (cos_phi * j_x) + (sin_phi * j_z)
            term_2 = (cos_phi * j_z) - (sin_phi * j_x)
            term_3 = 0.5 * (lor2_y1 - lor2_y2)
            term_4 = (term_1 * lor) - (term_2 * difatan) - (j_y * term_3)

            # ADD THIS SIDE'S CONTRIBUTION TO THE RUNNING TOTAL FOR THE CURRENT POLYGON
            b_x = b_x - (sin_phi * term_4)
            b_z = b_z + (cos_phi * term_4)
            b_y = b_y + (((sin_phi * j_x) - (cos_phi * j_z)) * term_3) + (j_y * difatan)

        # PROJECT ONTO THE MAGNETISATION DIRECTION TO GIVE THE TOTAL FIELD ANOMALY, AND CONVERT FROM A/m TO nT
        # NB: THE LEADING MINUS MATCHES THE VERTEX ORDERING CONVENTION OF talwani_and_heirtzler.nt
        total_field = ((b_x * np.cos(angle_a) * np.cos(angle_c - angle_b)) +
                       (b_z * np.sin(angle_a)) +
                       (-b_y * np.cos(angle_a) * np.sin(angle_c - angle_b)))
        n_t = n_t - (total_field * MU0_NT / (4.0 * np.pi))

    return n_t


def _split_by_extent(polygons):
    """
    Split *polygons* into those of infinite extent in the strike direction and those of finite extent.

    A polygon is of infinite extent when it lacks a ``'y1'``/``'y2'`` property or either of them is None. Elements
    that are None are dropped, as the two dimensional kernels also drop them.

    :param polygons: list of :func:`~gmgpy.polygon.Polygon`

    :returns: tuple (list, list) : the infinite extent polygons and the finite extent polygons
    """
    infinite = []
    finite = []

    for polygon in polygons:
        if polygon is None:
            continue
        elif polygon.props.get('y1') is None or polygon.props.get('y2') is None:
            infinite.append(polygon)
        else:
            finite.append(polygon)

    return infinite, finite


def _extents(polygon):
    """
    Return the strike direction extents of *polygon* as an ordered pair, in metres.

    :param polygon: :func:`~gmgpy.polygon.Polygon` carrying finite ``'y1'`` and ``'y2'`` properties

    :returns: tuple (float, float) : y1, y2 with y1 < y2
    """
    y1 = float(polygon.props['y1'])
    y2 = float(polygon.props['y2'])

    if y1 > y2:
        y1, y2 = y2, y1
    if y2 - y1 < SMALL:
        # A BODY OF ZERO THICKNESS HAS NO ANOMALY, BUT THE LOGARITHMS BELOW WOULD BE SINGULAR
        y2 = y1 + SMALL

    return y1, y2


def _vertex_pair(x, z, v, nverts, xp, zp):
    """
    Return the v'th polygon side, as vertex coordinates relative to the observation points.

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

    # NUDGE A VERTEX THAT SITS ON AN OBSERVATION POINT, WHERE THE LOGARITHMS BELOW ARE SINGULAR
    on_point = (np.abs(x1) < SMALL) & (np.abs(z1) < SMALL)
    x1 = np.where(on_point, SMALL, x1)
    z1 = np.where(on_point, SMALL, z1)

    on_point = (np.abs(x2) < SMALL) & (np.abs(z2) < SMALL)
    x2 = np.where(on_point, SMALL, x2)
    z2 = np.where(on_point, SMALL, z2)

    return x1, z1, x2, z2


def _rotated_frame(x1, z1, x2, z2):
    """
    Rotate a polygon side onto its own axis, following Rasmussen and Pedersen (1979).

    u is measured along the side and w perpendicular to it, so that the integral along the side has one variable of
    integration. NB: the side vector, and hence phi and s, do not depend on the observation point, so they are scalars
    while u1, u2 and w are arrays.

    :param x1: 1D array (float) : X coordinate of the first vertex, relative to the observation points
    :param z1: 1D array (float) : Z coordinate of the first vertex, relative to the observation points
    :param x2: 1D array (float) : X coordinate of the second vertex, relative to the observation points
    :param z2: 1D array (float) : Z coordinate of the second vertex, relative to the observation points

    :returns: tuple : cos_phi, sin_phi, s, u1, u2, w
    """
    x21 = x2 - x1
    z21 = z2 - z1

    # NUDGE A HORIZONTAL SIDE SO THAT THE SIGN OF phi IS WELL DEFINED
    z21 = np.where(z21 == 0.0, SMALL, z21)

    s = np.sqrt((x21 ** 2) + (z21 ** 2))
    phi = np.arctan2(z21, x21)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    u1 = (cos_phi * x1) + (sin_phi * z1)
    u2 = (cos_phi * x2) + (sin_phi * z2)
    w = (-sin_phi * x1) + (cos_phi * z1)

    return cos_phi, sin_phi, s, u1, u2, w


def _log_and_angle_terms(u1, u2, w, y1, y2):
    """
    The logarithmic and angular terms of the Rasmussen and Pedersen (1979) integral that the gravity and the magnetic
    kernels have in common.

    :param u1: 1D array (float) : Along side coordinate of the first vertex
    :param u2: 1D array (float) : Along side coordinate of the second vertex
    :param w: 1D array (float) : Perpendicular distance from the observation points to the side
    :param y1: float : Extent of the body in the strike direction, negative side, in metres
    :param y2: float : Extent of the body in the strike direction, positive side, in metres

    :returns: tuple : r1_y1, r1_y2, r2_y1, r2_y2, lor, difatan
    """
    d1 = (u1 ** 2) + (w ** 2)
    d2 = (u2 ** 2) + (w ** 2)

    r1 = np.sqrt(d1)
    r2 = np.sqrt(d2)

    # THE DISTANCES FROM THE OBSERVATION POINTS TO EACH END OF EACH VERTICAL EDGE OF THE BODY
    r1_y1 = np.sqrt(d1 + (y1 ** 2))
    r1_y2 = np.sqrt(d1 + (y2 ** 2))
    r2_y1 = np.sqrt(d2 + (y1 ** 2))
    r2_y2 = np.sqrt(d2 + (y2 ** 2))

    lor = np.log((r1 * (r2_y2 + y2)) / (r2 * (r1_y2 + y2))) + \
        np.log((r1 * (r2_y1 - y1)) / (r2 * (r1_y1 - y1)))

    difatan = _wrap_to_pi(np.arctan2(u2 * y2, w * r2_y2) - np.arctan2(u1 * y2, w * r1_y2)) - \
        _wrap_to_pi(np.arctan2(u2 * y1, w * r2_y1) - np.arctan2(u1 * y1, w * r1_y1))

    return r1_y1, r1_y2, r2_y1, r2_y2, lor, difatan


def _wrap_to_pi(angle):
    """
    Wrap a difference of two arctangents back into -pi -> pi.

    A polygon side that crosses the x axis of the observation point puts the difference outside that range, which
    would appear in the calculated anomaly as a jump of 2 pi.

    :param angle: 1D array (float) : Angle differences, in radians

    :returns: 1D array (float) : The same differences, in the range -pi -> pi
    """
    angle = np.where(angle < -np.pi, angle + (2.0 * np.pi), angle)
    angle = np.where(angle > np.pi, angle - (2.0 * np.pi), angle)

    return angle
