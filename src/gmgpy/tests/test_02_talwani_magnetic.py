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
