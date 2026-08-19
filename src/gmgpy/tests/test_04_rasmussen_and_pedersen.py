"""
Benchmark the Rasmussen & Pedersen (1979) 2.75D gravity and magnetic algorithms.

Three independent things are checked here:

1. That a body of infinite extent still gives the two dimensional answer, both by delegation (``y1``/``y2`` unset)
   and in the limit of a very large extent. This also pins the vertex ordering convention, since reversing it
   inverts the sign.
2. That a body of *finite* extent agrees with an exact three dimensional solution for the same body: Nagy (1966)
   for gravity, and a direct summation of point dipoles for magnetics.
3. Properties that must hold whatever the formulation: splitting a body along strike must sum to the whole, and
   shrinking the extent must reduce the anomaly.

The strike direction extents used are deliberately asymmetric about the section (``y1 != -y2``). A symmetric body
cannot detect an error in the sign of the along strike direction, because reflecting it leaves it unchanged.

**references**

Nagy, D. (1966). The gravitational attraction of a right rectangular prism. Geophysics, 31(2), 362-371.
doi:10.1190/1.1439779
"""
import numpy as np
import pytest

from gmgpy import bott
from gmgpy import rasmussen_and_pedersen
from gmgpy import talwani_and_heirtzler
from gmgpy.polygon import Polygon

DENSITY_CONTRAST = 500.0  # kg m^-3
SUSCEPTIBILITY = 0.01  # SI
FIELD_INTENSITY = 50000.0  # nT

# A BODY THIS LONG IS INDISTINGUISHABLE FROM ONE OF INFINITE EXTENT
EFFECTIVELY_INFINITE = 1.0e9

# MAGNETISATION GEOMETRIES: ANGLE A (DIP), ANGLE B (AZIMUTH), ANGLE C (PROFILE AZIMUTH), ALL IN DEGREES.
# THE CASES WITH C - B != 0 ARE THE ONES WITH A NON ZERO ALONG STRIKE MAGNETISATION COMPONENT.
MAGNETISATIONS = [(0.0, 0.0, 0.0),
                  (60.0, 0.0, 0.0),
                  (30.0, 10.0, 45.0),
                  (-45.0, 200.0, 300.0)]

# SEE test_01_bott_gravity.py FOR WHY BOTH A RELATIVE AND AN ABSOLUTE TERM ARE USED
RTOL = 1.0e-3

# THE DIPOLE SUM BELOW IS A MIDPOINT RULE, SO IT CARRIES ITS OWN O(1/N^2) DISCRETISATION ERROR AND NEEDS A LOOSER
# TOLERANCE THAN THE CLOSED FORM BENCHMARKS. AT N_DIPOLE_CELLS THE ERROR IS ABOUT 7e-5.
N_DIPOLE_CELLS = 48
RTOL_DIPOLE = 3.0e-4


def magnetic_props(angle_a, angle_b, angle_c, y1=None, y2=None):
    """
    Build the polygon properties for a magnetic calculation.

    :param angle_a: float : Magnetisation vector dip, in degrees
    :param angle_b: float : Magnetisation vector azimuth, in degrees
    :param angle_c: float : Model profile azimuth, in degrees
    :param y1: float : Extent in the strike direction, negative side, in metres. None for infinite extent
    :param y2: float : Extent in the strike direction, positive side, in metres. None for infinite extent

    :returns: dict : The polygon properties
    """
    return {'susceptibility': SUSCEPTIBILITY,
            'angle_a': angle_a,
            'angle_b': angle_b,
            'angle_c': angle_c,
            'f': FIELD_INTENSITY,
            'y1': y1,
            'y2': y2}


def analytical_prism_gz(xp, zp, prism, density_contrast):
    """
    Gravity anomaly of a right rectangular prism, after Nagy (1966).

    :param xp: 1D array (float) : X coordinates of the observation points
    :param zp: 1D array (float) : Z coordinates of the observation points
    :param prism: the prism geometry (x1, x2, z1, z2, y1, y2)
    :param density_contrast: float : Density contrast of the prism (kg m^-3)

    :returns: 1D array (float) : The calculated anomaly (mGal)
    """
    anomaly = np.zeros_like(xp)

    for n in range(len(xp)):
        # PRISM CORNERS RELATIVE TO THE CURRENT OBSERVATION POINT
        x = [prism.x1 - xp[n], prism.x2 - xp[n]]
        y = [prism.y1, prism.y2]
        z = [prism.z1 - zp[n], prism.z2 - zp[n]]

        # SUM THE EIGHT CORNERS WITH ALTERNATING SIGN
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    r = np.sqrt((x[i] ** 2) + (y[j] ** 2) + (z[k] ** 2))
                    anomaly[n] += ((-1.0) ** (i + j + k)) * \
                        ((x[i] * np.log(y[j] + r)) +
                         (y[j] * np.log(x[i] + r)) -
                         (z[k] * np.arctan2(x[i] * y[j], z[k] * r)))

    return anomaly * density_contrast * bott.G * bott.SI2MGAL


def dipole_sum_nt(xp, zp, prism, angle_a, angle_b, angle_c, n_cells):
    """
    Total field magnetic anomaly of a right rectangular prism, by direct summation of point dipoles.

    The prism is divided into n_cells**3 equal cells, each treated as a point dipole at its centre. This is a
    midpoint rule and converges on the exact answer as O(1/n_cells**2). It is used here in preference to a published
    closed form because superposition of point dipoles follows straight from the definition of magnetisation, so it
    is an independent check rather than a second transcription.

    :param xp: 1D array (float) : X coordinates of the observation points
    :param zp: 1D array (float) : Z coordinates of the observation points
    :param prism: the prism geometry (x1, x2, z1, z2, y1, y2)
    :param angle_a: float : Magnetisation vector dip, in degrees
    :param angle_b: float : Magnetisation vector azimuth, in degrees
    :param angle_c: float : Model profile azimuth, in degrees
    :param n_cells: int : Number of cells along each axis

    :returns: 1D array (float) : The calculated anomaly (nT)
    """
    a, b, c = np.radians([angle_a, angle_b, angle_c])

    # UNIT MAGNETISATION VECTOR IN THE MODEL FRAME, AND THE INDUCED MAGNETISATION MAGNITUDE IN A/m
    direction = np.array([np.cos(a) * np.cos(c - b), np.cos(a) * np.sin(c - b), np.sin(a)])
    magnetisation = SUSCEPTIBILITY * FIELD_INTENSITY / rasmussen_and_pedersen.MU0_NT

    # CELL CENTRES AND THE MOMENT OF ONE CELL
    edges = [np.linspace(lo, hi, n_cells, endpoint=False) + (0.5 * (hi - lo) / n_cells)
             for lo, hi in ((prism.x1, prism.x2), (prism.y1, prism.y2), (prism.z1, prism.z2))]
    cell_volume = ((prism.x2 - prism.x1) * (prism.y2 - prism.y1) * (prism.z2 - prism.z1)) / (n_cells ** 3)
    cx, cy, cz = [a.ravel() for a in np.meshgrid(*edges, indexing='ij')]
    moment = magnetisation * cell_volume * direction

    anomaly = np.zeros_like(xp)

    for n in range(len(xp)):
        # VECTOR FROM EACH DIPOLE TO THE OBSERVATION POINT, WHICH LIES IN THE PLANE OF THE SECTION (y = 0)
        rx = xp[n] - cx
        ry = -cy
        rz = zp[n] - cz
        r_squared = (rx ** 2) + (ry ** 2) + (rz ** 2)
        r = np.sqrt(r_squared)

        # B = (mu_0 / 4 pi) (3 (m.rhat) rhat - m) / r^3, WITH THE LEADING FACTOR IN nT PER A m^2 BEING 100
        m_dot_r = (moment[0] * rx) + (moment[1] * ry) + (moment[2] * rz)
        b_x = 100.0 * (((3.0 * m_dot_r * rx) / r_squared) - moment[0]) / (r * r_squared)
        b_y = 100.0 * (((3.0 * m_dot_r * ry) / r_squared) - moment[1]) / (r * r_squared)
        b_z = 100.0 * (((3.0 * m_dot_r * rz) / r_squared) - moment[2]) / (r * r_squared)

        # PROJECT ONTO THE MAGNETISATION DIRECTION, MATCHING THE INDUCED ONLY CONVENTION OF gmg
        anomaly[n] = np.sum((b_x * direction[0]) + (b_y * direction[1]) + (b_z * direction[2]))

    return anomaly


def assert_close(calculated, expected, rtol=RTOL):
    """Compare two anomalies with a combined relative and peak scaled absolute tolerance."""
    np.testing.assert_allclose(calculated, expected, rtol=rtol,
                               atol=rtol * np.max(np.abs(expected)))


def test_gravity_without_extents_delegates(cylinder, profile):
    """
    Test a polygon with no strike extents set is passed through to the 2D gravity algorithm unchanged
    """
    xp, zp = profile
    props = {'density': DENSITY_CONTRAST}

    calculated = rasmussen_and_pedersen.gz(xp, zp, [Polygon(cylinder.nodes, dict(props, y1=None, y2=None))])
    expected = bott.gz(xp, zp, [Polygon(cylinder.nodes, props)])

    # DELEGATION MUST BE EXACT, NOT MERELY CLOSE, SO THAT EXISTING MODELS ARE UNAFFECTED
    np.testing.assert_array_equal(calculated, expected)


def test_magnetic_without_extents_delegates(cylinder, profile):
    """
    Test a polygon with no strike extents set is passed through to the 2D magnetic algorithm unchanged
    """
    xp, zp = profile
    angle_a, angle_b, angle_c = MAGNETISATIONS[2]
    props = magnetic_props(angle_a, angle_b, angle_c)

    calculated = rasmussen_and_pedersen.nt(xp, zp, [Polygon(cylinder.nodes, dict(props))])
    del props['y1'], props['y2']
    expected = talwani_and_heirtzler.nt(xp, zp, [Polygon(cylinder.nodes, props)])

    # DELEGATION MUST BE EXACT, NOT MERELY CLOSE, SO THAT EXISTING MODELS ARE UNAFFECTED
    np.testing.assert_array_equal(calculated, expected)


def test_gravity_reduces_to_2d(cylinder, profile):
    """
    Test the 2.75D gravity algorithm reproduces the 2D algorithm for a body of very large strike extent
    """
    xp, zp = profile
    props = {'density': DENSITY_CONTRAST}

    calculated = rasmussen_and_pedersen.gz(
        xp, zp, [Polygon(cylinder.nodes, dict(props, y1=-EFFECTIVELY_INFINITE, y2=EFFECTIVELY_INFINITE))])
    expected = bott.gz(xp, zp, [Polygon(cylinder.nodes, props)])

    assert_close(calculated, expected)


@pytest.mark.parametrize('angle_a, angle_b, angle_c', MAGNETISATIONS)
def test_magnetic_reduces_to_2d(cylinder, profile, angle_a, angle_b, angle_c):
    """
    Test the 2.75D magnetic algorithm reproduces the 2D algorithm for a body of very large strike extent

    The along strike component of magnetisation contributes nothing in two dimensions, so this must hold for every
    magnetisation direction, not just those lying in the plane of the section.
    """
    xp, zp = profile
    props = magnetic_props(angle_a, angle_b, angle_c,
                           y1=-EFFECTIVELY_INFINITE, y2=EFFECTIVELY_INFINITE)

    calculated = rasmussen_and_pedersen.nt(xp, zp, [Polygon(cylinder.nodes, props)])
    expected = talwani_and_heirtzler.nt(
        xp, zp, [Polygon(cylinder.nodes, magnetic_props(angle_a, angle_b, angle_c))])

    assert_close(calculated, expected)


def test_gravity_finite_prism(prism, profile):
    """
    Test the 2.75D gravity algorithm reproduces the analytical solution for a finite length prism
    """
    xp, zp = profile
    polygons = [Polygon(prism.nodes, {'density': DENSITY_CONTRAST, 'y1': prism.y1, 'y2': prism.y2})]

    calculated = rasmussen_and_pedersen.gz(xp, zp, polygons)
    expected = analytical_prism_gz(xp, zp, prism, DENSITY_CONTRAST)

    assert_close(calculated, expected)


@pytest.mark.parametrize('angle_a, angle_b, angle_c', MAGNETISATIONS)
def test_magnetic_finite_prism(prism, coarse_profile, angle_a, angle_b, angle_c):
    """
    Test the 2.75D magnetic algorithm reproduces a direct dipole summation for a finite length prism
    """
    xp, zp = coarse_profile
    props = magnetic_props(angle_a, angle_b, angle_c, y1=prism.y1, y2=prism.y2)

    calculated = rasmussen_and_pedersen.nt(xp, zp, [Polygon(prism.nodes, props)])
    expected = dipole_sum_nt(xp, zp, prism, angle_a, angle_b, angle_c, N_DIPOLE_CELLS)

    assert_close(calculated, expected, rtol=RTOL_DIPOLE)


def test_gravity_extents_are_additive(cylinder, profile):
    """
    Test splitting a body along strike gives the same gravity anomaly as the whole body
    """
    xp, zp = profile
    props = {'density': DENSITY_CONTRAST}

    whole = rasmussen_and_pedersen.gz(xp, zp, [Polygon(cylinder.nodes, dict(props, y1=-3000.0, y2=5000.0))])
    parts = rasmussen_and_pedersen.gz(xp, zp, [Polygon(cylinder.nodes, dict(props, y1=-3000.0, y2=1000.0)),
                                               Polygon(cylinder.nodes, dict(props, y1=1000.0, y2=5000.0))])

    assert_close(parts, whole)


@pytest.mark.parametrize('angle_a, angle_b, angle_c', MAGNETISATIONS)
def test_magnetic_extents_are_additive(cylinder, profile, angle_a, angle_b, angle_c):
    """
    Test splitting a body along strike gives the same magnetic anomaly as the whole body
    """
    xp, zp = profile

    whole = rasmussen_and_pedersen.nt(
        xp, zp, [Polygon(cylinder.nodes, magnetic_props(angle_a, angle_b, angle_c, -3000.0, 5000.0))])
    parts = rasmussen_and_pedersen.nt(
        xp, zp, [Polygon(cylinder.nodes, magnetic_props(angle_a, angle_b, angle_c, -3000.0, 1000.0)),
                 Polygon(cylinder.nodes, magnetic_props(angle_a, angle_b, angle_c, 1000.0, 5000.0))])

    assert_close(parts, whole)


def test_gravity_converges_on_2d_from_below(cylinder, profile):
    """
    Test a shortening body gives a monotonically smaller gravity anomaly than an infinite one

    This is the whole point of 2.75D modelling: an infinite body over predicts the anomaly, which biases a fitted
    density contrast low.
    """
    xp, zp = profile
    props = {'density': DENSITY_CONTRAST}

    infinite = np.max(np.abs(bott.gz(xp, zp, [Polygon(cylinder.nodes, props)])))
    peaks = [np.max(np.abs(rasmussen_and_pedersen.gz(
        xp, zp, [Polygon(cylinder.nodes, dict(props, y1=-extent, y2=extent))])))
        for extent in (1.0e3, 5.0e3, 2.0e4, 1.0e5, 1.0e6)]

    assert np.all(np.diff(peaks) > 0.0), 'anomaly must grow as the body lengthens'
    assert np.all(np.array(peaks) < infinite), 'a finite body must not out predict an infinite one'
    assert peaks[-1] == pytest.approx(infinite, rel=RTOL)


def test_mixed_extents_superpose(cylinder, prism, profile):
    """
    Test a model holding both an infinite and a finite body is the sum of the two calculated separately
    """
    xp, zp = profile
    infinite_body = Polygon(cylinder.nodes, {'density': DENSITY_CONTRAST})
    finite_body = Polygon(prism.nodes, {'density': -DENSITY_CONTRAST, 'y1': prism.y1, 'y2': prism.y2})

    together = rasmussen_and_pedersen.gz(xp, zp, [infinite_body, finite_body])
    separately = (rasmussen_and_pedersen.gz(xp, zp, [infinite_body]) +
                  rasmussen_and_pedersen.gz(xp, zp, [finite_body]))

    assert_close(together, separately)
