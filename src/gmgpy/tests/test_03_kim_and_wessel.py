"""
Benchmark the Kim & Wessel (2016) vertical gravity gradient (VGG) algorithm
against the analytical solution for a 2D infinite horizontal cylinder.

The VGG is the vertical derivative of the gravity anomaly used in
test_01_bott_gravity.py, in the z -> **DOWN** sense, so it is positive directly
above a body of positive density contrast.

**references**

Analytical solution for buried cylinder: Garland (1965) Pg. 70
"""
import numpy as np

from gmgpy import bott, kim_and_wessel
from gmgpy.polygon import Polygon

DENSITY_CONTRAST = 250.0  # kg/m^3

# SEE test_01_bott_gravity.py FOR WHY BOTH A RELATIVE AND AN ABSOLUTE TERM ARE USED
RTOL = 1.0e-3


def analytical_vgg(xp, cylinder, density_contrast):
    """
    Vertical gravity gradient of a 2D infinite horizontal cylinder.

    :param xp: 1D array (float) : X coordinates of the observation points
    :param cylinder: the cylinder geometry (radius, centre_x, centre_z)
    :param density_contrast: float : Density contrast of cylinder (kg/m^3)

    :returns: 1D array (float) : The calculated anomaly (Eotvos)
    """
    dx = xp - cylinder.centre_x
    z2 = cylinder.centre_z ** 2
    return ((2. * np.pi * bott.G * (cylinder.radius ** 2) * density_contrast) *
            ((z2 - (dx ** 2)) / (((dx ** 2) + z2) ** 2))) * kim_and_wessel.SI_TO_EOTVOS


def test_kim_vgg(cylinder, profile):
    """
    Test the kim and wessel vgg algorithm reproduces the analytical solution
    """
    xp, zp = profile
    polygons = [Polygon(cylinder.nodes, {'density': DENSITY_CONTRAST})]

    calculated = kim_and_wessel.gz(xp, zp, polygons)
    expected = analytical_vgg(xp, cylinder, DENSITY_CONTRAST)

    np.testing.assert_allclose(calculated, expected, rtol=RTOL,
                               atol=RTOL * np.max(np.abs(expected)))
