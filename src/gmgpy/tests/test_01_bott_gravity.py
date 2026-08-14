"""
Benchmark the Bott (1969) gravity algorithm against the analytical solution for
a 2D infinite horizontal cylinder.

**references**

Analytical solution for buried cylinder: Garland (1965) Pg. 70
"""
import numpy as np

from gmgpy import bott
from gmgpy.polygon import Polygon

DENSITY_CONTRAST = 250.0  # kg/m^3

# THE POLYGON IS A DISCRETISATION OF A CIRCLE, SO AGREEMENT IS CLOSE, NOT EXACT.
# THE ABSOLUTE TERM IS SCALED TO THE PEAK AMPLITUDE - RELATIVE ERROR ALONE IS
# UNBOUNDED WHERE AN ANOMALY CROSSES ZERO.
RTOL = 1.0e-3


def analytical_gz(xp, cylinder, density_contrast):
    """
    Gravity anomaly of a 2D infinite horizontal cylinder.

    :param xp: 1D array (float) : X coordinates of the observation points
    :param cylinder: the cylinder geometry (radius, centre_x, centre_z)
    :param density_contrast: float : Density contrast of cylinder (kg/m^3)

    :returns: 1D array (float) : The calculated anomaly (mGal)
    """
    dx = xp - cylinder.centre_x
    return ((2. * np.pi * bott.G * (cylinder.radius ** 2) * density_contrast) *
            (cylinder.centre_z / ((dx ** 2) + (cylinder.centre_z ** 2)))) * bott.SI2MGAL


def test_bott_gravity(cylinder, profile):
    """
    Test the bott gravity algorithm reproduces the analytical solution
    """
    xp, zp = profile
    polygons = [Polygon(cylinder.nodes, {'density': DENSITY_CONTRAST})]

    calculated = bott.gz(xp, zp, polygons)
    expected = analytical_gz(xp, cylinder, DENSITY_CONTRAST)

    np.testing.assert_allclose(calculated, expected, rtol=RTOL,
                               atol=RTOL * np.max(np.abs(expected)))
