"""
Check every forward kernel against bodies placed relative to the calculation elevation.

The other test modules put the whole body safely below the observation points. Real models do not: a layer
drawn from the surface sits *on* the calculation elevation, topography puts one *above* it, and an outcropping
body *crosses* it. Each of those used to produce nan.

Two things are asserted here:

1. That no kernel returns a non-finite value, or warns, for any of those geometries. One nan is enough to blank
   an entire anomaly curve, because the polygons are summed by superposition.
2. That the values are right, by comparing against references built independently of gmg:

   * gravity, against a sum of infinite horizontal **line masses**, ``g = 2*G*lambda*dz/(dx^2+dz^2)``, which is
     exact in the strike direction and valid for either sign of dz;
   * the VGG, against a central difference of :func:`~gmgpy.bott.gz`, which the gravity check above validates.

.. note:: Observation points lying *inside* a polygon are excluded from the magnetic comparisons. Both magnetic
          kernels are known to be wrong there and gmg marks those samples rather than trusting them; see
          ``docs/_sources/manual_magnetic_algorithm.txt``.
"""
import warnings

import numpy as np
import pytest

from gmgpy import bott
from gmgpy import kim_and_wessel
from gmgpy import rasmussen_and_pedersen
from gmgpy import talwani_and_heirtzler
from gmgpy.polygon import Polygon
from gmgpy.tests.conftest import get_rectangle_points

DENSITY_CONTRAST = 500.0  # kg m^-3
SUSCEPTIBILITY = 0.01  # SI
FIELD_INTENSITY = 50000.0  # nT
ANGLE_A, ANGLE_B, ANGLE_C = 45.0, 10.0, 30.0

# THE TEST BODY, IN METRES. ITS X EDGES FALL EXACTLY ON OBSERVATION POINTS, WHICH IS WHAT MAKES A NODE PLACED AT
# THE CALCULATION ELEVATION COINCIDE WITH ONE OF THEM
BODY_X1, BODY_X2 = 1500.0, 2500.0
PROFILE_START, PROFILE_END, PROFILE_NUM = -2000.0, 6000.0, 81

# DEPTH RANGES RELATIVE TO THE CALCULATION ELEVATION, WHICH IS ZERO THROUGHOUT
DEPTH_RANGES = [(1000.0, 2500.0),     # WHOLLY BELOW, THE ORDINARY CASE
                (-2500.0, -1000.0),   # WHOLLY ABOVE, E.G. TOPOGRAPHY
                (0.0, 2500.0),        # TOP EDGE EXACTLY ON IT, E.G. A LAYER DRAWN FROM THE SURFACE
                (-500.0, 2500.0),     # CROSSING IT, E.G. AN OUTCROPPING BODY
                (-1500.0, 1500.0)]    # STRADDLING IT SYMMETRICALLY

RTOL = 1.0e-3


@pytest.fixture
def elevation_profile():
    """X and Z coordinates of the observation points, at the default calculation elevation of zero."""
    xp = np.linspace(PROFILE_START, PROFILE_END, num=PROFILE_NUM)
    return xp, np.zeros_like(xp)


def body(z1, z2):
    """The test body's polygon nodes for a given depth range."""
    return get_rectangle_points(BODY_X1, BODY_X2, z1, z2)


def gravity_polygon(z1, z2):
    return Polygon(body(z1, z2), {'density': DENSITY_CONTRAST})


def magnetic_polygon(z1, z2, y1=None, y2=None):
    return Polygon(body(z1, z2), {'susceptibility': SUSCEPTIBILITY,
                                  'angle_a': ANGLE_A, 'angle_b': ANGLE_B, 'angle_c': ANGLE_C,
                                  'f': FIELD_INTENSITY, 'y1': y1, 'y2': y2})


def analytical_line_mass_gz(xp, z1, z2, n_cells=400):
    """
    Gravity anomaly of the test body, by summing infinite horizontal line masses.

    The attraction of a line mass of linear density lambda at horizontal offset dx and depth dz is
    ``2*G*lambda*dz/(dx^2 + dz^2)``, so the strike direction is integrated exactly and only the cross section is
    discretised. Unlike a closed form for a prism it is valid for either sign of dz, which is the whole point
    here.

    :param xp: 1D array (float) : X coordinates of the observation points
    :param z1: float : Depth of the top of the body, in metres
    :param z2: float : Depth of the base of the body, in metres
    :param n_cells: int : Number of cells along each axis of the cross section

    :returns: 1D array (float) : The calculated anomaly (mGal)
    """
    xs = np.linspace(BODY_X1, BODY_X2, n_cells, endpoint=False) + (0.5 * (BODY_X2 - BODY_X1) / n_cells)
    zs = np.linspace(z1, z2, n_cells, endpoint=False) + (0.5 * (z2 - z1) / n_cells)
    cell_area = ((BODY_X2 - BODY_X1) / n_cells) * ((z2 - z1) / n_cells)
    cx, cz = [axis.ravel() for axis in np.meshgrid(xs, zs, indexing='ij')]

    anomaly = np.array([np.sum(cz / (((cx - xo) ** 2) + (cz ** 2))) for xo in xp])

    return anomaly * 2.0 * bott.G * DENSITY_CONTRAST * cell_area * bott.SI2MGAL


def finite_difference_vgg(xp, z1, z2, step=0.5):
    """The vertical gravity gradient as a central difference of the gravity kernel, in Eotvos."""
    above = bott.gz(xp, np.full_like(xp, -step), [gravity_polygon(z1, z2)])
    below = bott.gz(xp, np.full_like(xp, step), [gravity_polygon(z1, z2)])

    return ((below - above) / (2.0 * step)) * 1.0e4  # mGal/m -> Eotvos


def outside_body(xp):
    """Mask of observation points that do not lie within the body's horizontal span."""
    return ~((xp >= BODY_X1) & (xp <= BODY_X2))


@pytest.mark.parametrize('z1, z2', DEPTH_RANGES)
def test_no_kernel_returns_a_non_finite_value(elevation_profile, z1, z2):
    """
    Test every kernel returns finite values wherever the body sits relative to the calculation elevation

    A single nan blanks a whole anomaly curve, since the polygons are summed by superposition.
    """
    xp, zp = elevation_profile

    with warnings.catch_warnings():
        warnings.simplefilter('error')  # A DIVIDE BY ZERO OR INVALID VALUE MUST FAIL THE TEST, NOT PASS QUIETLY
        calculated = {'bott': bott.gz(xp, zp, [gravity_polygon(z1, z2)]),
                      'kim_and_wessel': kim_and_wessel.gz(xp, zp, [gravity_polygon(z1, z2)]),
                      'talwani_and_heirtzler': talwani_and_heirtzler.nt(xp, zp, [magnetic_polygon(z1, z2)]),
                      'rasmussen_and_pedersen': rasmussen_and_pedersen.nt(
                          xp, zp, [magnetic_polygon(z1, z2, -20000.0, 20000.0)])}

    for name, anomaly in calculated.items():
        assert np.all(np.isfinite(anomaly)), '%s returned %d non-finite values' % (
            name, np.sum(~np.isfinite(anomaly)))


# THE SYMMETRIC BODY IS EXCLUDED: ITS ANOMALY CANCELS TO ZERO, SO THERE IS NOTHING TO COMPARE RELATIVELY AGAINST.
# IT GETS ITS OWN TEST BELOW
@pytest.mark.parametrize('z1, z2', [pair for pair in DEPTH_RANGES if pair != (-1500.0, 1500.0)])
def test_gravity_matches_the_line_mass_solution(elevation_profile, z1, z2):
    """
    Test the gravity kernel is correct for a body above, below, on or crossing the calculation elevation
    """
    xp, zp = elevation_profile

    calculated = bott.gz(xp, zp, [gravity_polygon(z1, z2)])
    expected = analytical_line_mass_gz(xp, z1, z2)

    # THE REFERENCE HAS AN INTEGRABLE SINGULARITY AT POINTS INSIDE THE BODY, WHERE THE MIDPOINT RULE CONVERGES
    # ONLY SLOWLY, SO IT IS COMPARED WHERE IT IS TRUSTWORTHY
    keep = outside_body(xp)
    np.testing.assert_allclose(calculated[keep], expected[keep], rtol=RTOL,
                               atol=RTOL * np.max(np.abs(expected[keep])))


def test_a_body_straddling_the_elevation_symmetrically_cancels(elevation_profile):
    """
    Test a body centred on the calculation elevation produces no gravity anomaly

    The mass above the observation points attracts upwards exactly as much as the mass below attracts downwards,
    so the two cancel. This is worth asserting on its own because it is the one geometry whose correct answer is
    zero, which makes it useless for a relative comparison but a sharp check of the signs.
    """
    xp, zp = elevation_profile

    calculated = bott.gz(xp, zp, [gravity_polygon(-1500.0, 1500.0)])

    # SCALED AGAINST THE SAME BODY PLACED BELOW THE ELEVATION, SINCE THE CANCELLED TOTAL CARRIES NO SCALE OF ITS OWN
    reference_magnitude = np.max(np.abs(bott.gz(xp, zp, [gravity_polygon(1000.0, 2500.0)])))

    np.testing.assert_allclose(calculated, np.zeros_like(calculated),
                               atol=1.0e-9 * reference_magnitude)


@pytest.mark.parametrize('z1, z2', [(1000.0, 2500.0), (-2500.0, -1000.0), (-500.0, 2500.0)])
def test_vgg_matches_a_central_difference_of_the_gravity(elevation_profile, z1, z2):
    """
    Test the VGG kernel is correct for a body above, below or crossing the calculation elevation

    Taking the logarithm of a node's depth used to return nan for every one of these but the first.
    """
    xp, zp = elevation_profile

    calculated = kim_and_wessel.gz(xp, zp, [gravity_polygon(z1, z2)])
    expected = finite_difference_vgg(xp, z1, z2)

    np.testing.assert_allclose(calculated, expected, rtol=RTOL,
                               atol=RTOL * np.max(np.abs(expected)))


def test_vgg_on_the_body_surface_takes_the_exterior_limit(elevation_profile):
    """
    Test the VGG at points lying exactly on the body's top surface is the limit from outside it

    The vertical gradient is genuinely discontinuous across a mass surface, by 4*pi*G*rho, so a point on it has
    two answers. Returning the exterior one is the meaningful choice for a surface observation, and is what a
    central difference cannot give because it straddles the jump.
    """
    xp, zp = elevation_profile

    # THE BODY'S TWO TOP CORNERS ARE EXCLUDED. THE LIMIT THERE DEPENDS ON THE DIRECTION OF APPROACH, SO A CORNER
    # HAS NO SINGLE EXTERIOR VALUE TO COMPARE AGAINST
    on_surface = (xp > BODY_X1) & (xp < BODY_X2)

    calculated = kim_and_wessel.gz(xp, zp, [gravity_polygon(0.0, 2500.0)])

    # THE LIMIT APPROACHED FROM ABOVE THE BODY, WHICH IS OUTSIDE IT
    approached = kim_and_wessel.gz(xp, np.full_like(xp, -0.01), [gravity_polygon(0.0, 2500.0)])

    np.testing.assert_allclose(calculated[on_surface], approached[on_surface], rtol=1.0e-3)

    # AND THE JUMP ACROSS THE SURFACE IS THE EXPECTED 4*pi*G*rho
    below = kim_and_wessel.gz(xp, np.full_like(xp, 0.01), [gravity_polygon(0.0, 2500.0)])
    jump = calculated[on_surface] - below[on_surface]
    expected_jump = 4.0 * np.pi * bott.G * DENSITY_CONTRAST * kim_and_wessel.SI_TO_EOTVOS

    np.testing.assert_allclose(jump, expected_jump, rtol=1.0e-3)


def test_a_node_lying_on_an_observation_point_is_handled(elevation_profile):
    """
    Test a node placed exactly on an observation point does not produce a non-finite value

    The kernels all take the logarithm of the distance from a node to the observation point, which is singular
    when the two coincide. It happens whenever a node sits at the calculation elevation and its x falls on one
    of the observation points, which is easily done since models are drawn from the surface at round positions.
    """
    xp, zp = elevation_profile
    assert BODY_X1 in xp and BODY_X2 in xp, 'the body edges must land on observation points for this test'

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        gravity = bott.gz(xp, zp, [gravity_polygon(0.0, 2500.0)])
        vgg = kim_and_wessel.gz(xp, zp, [gravity_polygon(0.0, 2500.0)])
        magnetic = talwani_and_heirtzler.nt(xp, zp, [magnetic_polygon(0.0, 2500.0)])

    assert np.all(np.isfinite(gravity))
    assert np.all(np.isfinite(vgg))
    assert np.all(np.isfinite(magnetic))


@pytest.mark.parametrize('z1, z2', DEPTH_RANGES)
def test_magnetic_kernels_agree_outside_the_body(elevation_profile, z1, z2):
    """
    Test the 2D and 2.75D magnetic kernels agree outside the body, whatever its depth range

    They are known to disagree at observation points inside the polygon; gmg marks those rather than trusting
    them. Everywhere else the 2.75D kernel at a very large strike extent must reproduce the 2D one.
    """
    xp, zp = elevation_profile
    keep = outside_body(xp)

    two_d = talwani_and_heirtzler.nt(xp, zp, [magnetic_polygon(z1, z2)])
    two_seven_five_d = rasmussen_and_pedersen.nt(
        xp, zp, [magnetic_polygon(z1, z2, -1.0e9, 1.0e9)])

    np.testing.assert_allclose(two_seven_five_d[keep], two_d[keep], rtol=RTOL,
                               atol=RTOL * np.max(np.abs(two_d[keep])))
