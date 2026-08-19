"""
Benchmark the Talwani & Heirtzler (1964) magnetic algorithm against the
analytical solution for a 2D infinite horizontal cylinder.

**references**

Sleep & Fujita 1997. Pg. 223 eq. 6.45
"""
import numpy as np

from gmgpy import talwani_and_heirtzler
from gmgpy.polygon import Polygon

SUSCEPTIBILITY = 0.001
FIELD_INTENSITY = 50000.0  # nT
ANGLE_A = 0.0  # MAGNETISATION VECTOR ANGLE FROM HORIZONTAL, DEGREES
ANGLE_B = 0.0  # MAGNETISATION AZIMUTH FROM GEOGRAPHIC NORTH, DEGREES
ANGLE_C = 0.0  # PROFILE AZIMUTH FROM GEOGRAPHIC NORTH, DEGREES

# SEE test_01_bott_gravity.py FOR WHY BOTH A RELATIVE AND AN ABSOLUTE TERM ARE USED
RTOL = 1.0e-3


def analytical_nt(xp, cylinder, susceptibility, field_intensity):
    """
    Magnetic anomaly of a 2D infinite horizontal cylinder.

    :param xp: 1D array (float) : X coordinates of the observation points
    :param cylinder: the cylinder geometry (radius, centre_x, centre_z)
    :param susceptibility: float : Susceptibility contrast of cylinder (SI)
    :param field_intensity: float : Total regional field intensity (nT)

    :returns: 1D array (float) : The calculated anomaly (nT)
    """
    # INDUCED MAGNETISATION CONTRAST
    magnetisation = susceptibility * field_intensity

    dx = xp - cylinder.centre_x
    z2 = cylinder.centre_z ** 2
    return magnetisation * (cylinder.radius ** 2) * \
        (((dx ** 2) - z2) / (2. * (((dx ** 2) + z2) ** 2)))


def test_talwani_magnetic(cylinder, profile):
    """
    Test the talwani and heirtzler magnetic algorithm reproduces the analytical
    solution
    """
    xp, zp = profile
    polygons = [Polygon(cylinder.nodes, {'susceptibility': SUSCEPTIBILITY,
                                         'angle_a': ANGLE_A,
                                         'angle_b': ANGLE_B,
                                         'angle_c': ANGLE_C,
                                         'f': FIELD_INTENSITY})]

    calculated = talwani_and_heirtzler.nt(xp, zp, polygons)
    expected = analytical_nt(xp, cylinder, SUSCEPTIBILITY, FIELD_INTENSITY)

    np.testing.assert_allclose(calculated, expected, rtol=RTOL,
                               atol=RTOL * np.max(np.abs(expected)))


def magnetic_props(angle_a, angle_b, angle_c, **extra):
    """Build the polygon properties for a magnetic calculation."""
    props = {'susceptibility': SUSCEPTIBILITY,
             'angle_a': angle_a,
             'angle_b': angle_b,
             'angle_c': angle_c,
             'f': FIELD_INTENSITY}
    props.update(extra)
    return props


def test_layers_superpose(cylinder, prism, profile):
    """
    Test two layers modelled together give the same anomaly as the sum of the two modelled alone

    The two layers deliberately have *different* magnetisation directions. The projection onto the inducing
    field direction has to be done per layer for this to hold; doing it once against sums accumulated over
    every layer projects all of them onto whichever layer happened to be processed last.
    """
    xp, zp = profile
    a = Polygon(cylinder.nodes, magnetic_props(30.0, 0.0, 0.0))
    b = Polygon(prism.nodes, magnetic_props(-50.0, 40.0, 70.0))

    together = talwani_and_heirtzler.nt(xp, zp, [a, b])
    apart = (talwani_and_heirtzler.nt(xp, zp, [a]) +
             talwani_and_heirtzler.nt(xp, zp, [b]))

    np.testing.assert_allclose(together, apart, rtol=RTOL,
                               atol=RTOL * np.max(np.abs(apart)))


def test_layer_order_does_not_matter(cylinder, prism, profile):
    """
    Test the calculated anomaly does not depend on the order the layers are given in
    """
    xp, zp = profile
    a = Polygon(cylinder.nodes, magnetic_props(30.0, 0.0, 0.0))
    b = Polygon(prism.nodes, magnetic_props(-50.0, 40.0, 70.0))

    np.testing.assert_allclose(talwani_and_heirtzler.nt(xp, zp, [a, b]),
                               talwani_and_heirtzler.nt(xp, zp, [b, a]))


def test_layer_without_susceptibility_is_skipped(cylinder, profile):
    """
    Test a polygon carrying no susceptibility is ignored rather than raising

    This matches the behaviour of :func:`~gmgpy.bott.gz` for a polygon carrying no density.
    """
    xp, zp = profile
    gravity_only = Polygon(cylinder.nodes, {'density': 500.0})

    calculated = talwani_and_heirtzler.nt(xp, zp, [gravity_only])

    np.testing.assert_array_equal(calculated, np.zeros_like(xp))


def test_remanence_projects_onto_the_inducing_field(cylinder, profile):
    """
    Test the anomaly is projected onto the inducing field direction, not the magnetisation direction

    The two are the same only when magnetisation is induced, so supplying an inducing field direction equal
    to the magnetisation direction must reproduce the induced only default, and supplying a different one
    must not.
    """
    xp, zp = profile
    angle_a, angle_b, angle_c = 30.0, 10.0, 45.0

    induced = talwani_and_heirtzler.nt(
        xp, zp, [Polygon(cylinder.nodes, magnetic_props(angle_a, angle_b, angle_c))])
    stated = talwani_and_heirtzler.nt(
        xp, zp, [Polygon(cylinder.nodes, magnetic_props(angle_a, angle_b, angle_c,
                                                        inclination=angle_a, declination=angle_b))])
    remanent = talwani_and_heirtzler.nt(
        xp, zp, [Polygon(cylinder.nodes, magnetic_props(angle_a, angle_b, angle_c,
                                                        inclination=60.0, declination=-5.0))])

    np.testing.assert_allclose(stated, induced)
    assert not np.allclose(remanent, induced), 'remanence must change the projected anomaly'
