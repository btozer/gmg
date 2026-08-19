"""
Shared geometry for the gmg analytical benchmark tests.

Each test compares one of the gmg forward modelling algorithms against the
closed form solution for a 2D infinite horizontal cylinder. The cylinder is
represented to the algorithms as a polygon of discrete nodes, so agreement is
close but not exact - the node count below sets how close.

.. note:: The coordinate system is z -> **DOWN**, matching the algorithms.

.. note:: The node ordering produced by :func:`get_circle_points` is the
          ordering the algorithms expect; reversing it inverts the sign of the
          calculated anomaly.
"""
from types import SimpleNamespace

import numpy as np
import pytest

# CYLINDER GEOMETRY (SI UNITS)
RADIUS = 1000.0
CENTRE_X = 2000.0
CENTRE_Z = 2000.0
N_NODES = 360

# PRISM GEOMETRY FOR THE 2.75D BENCHMARKS (SI UNITS). THE STRIKE DIRECTION EXTENTS ARE DELIBERATELY ASYMMETRIC ABOUT
# THE SECTION, SINCE A SYMMETRIC BODY CANNOT DETECT AN ERROR IN THE SIGN OF THE ALONG STRIKE DIRECTION.
PRISM_X1 = 1500.0
PRISM_X2 = 2500.0
PRISM_Z1 = 1000.0
PRISM_Z2 = 2500.0
PRISM_Y1 = -2000.0
PRISM_Y2 = 3000.0

# OBSERVATION PROFILE (SI UNITS)
PROFILE_START = 0.0
PROFILE_END = 4000.0
PROFILE_NUM = 4001

# A COARSER PROFILE FOR THE BENCHMARKS WHOSE REFERENCE SOLUTION COSTS PER OBSERVATION POINT
COARSE_PROFILE_NUM = 21


def get_circle_points(r, x, z, N):
    """
    :param r: float : Radius of circle (m)
    :param x: float : X position of circle center
    :param z: float : Z position of circle center
    :param N: int : Number of nodes defining the 'circle'

    :returns: 2D array (float): C1:X coordinate; C2:Z coordinate

    """
    angle = np.linspace(0., 2.*np.pi, N)
    output_points = np.zeros(shape=(len(angle), 2))

    for i in range(len(angle)):
        output_points[i, 0] = x + (r * np.cos(angle[i]))
        output_points[i, 1] = z + (r * np.sin(angle[i]))

    return output_points


@pytest.fixture
def cylinder():
    """The test cylinder: its defining parameters and its polygon nodes."""
    return SimpleNamespace(
        radius=RADIUS,
        centre_x=CENTRE_X,
        centre_z=CENTRE_Z,
        nodes=get_circle_points(RADIUS, CENTRE_X, CENTRE_Z, N_NODES),
    )


def get_rectangle_points(x1, x2, z1, z2):
    """
    :param x1: float : Minimum X coordinate of the rectangle
    :param x2: float : Maximum X coordinate of the rectangle
    :param z1: float : Minimum Z coordinate of the rectangle
    :param z2: float : Maximum Z coordinate of the rectangle

    :returns: 2D array (float): C1:X coordinate; C2:Z coordinate

    .. note:: The nodes are ordered in the same rotational sense as :func:`get_circle_points`, which is the ordering
              the algorithms expect.
    """
    return np.array([[x2, z1], [x2, z2], [x1, z2], [x1, z1]], dtype=float)


@pytest.fixture
def profile():
    """X and Z coordinates of the observation points."""
    xp = np.linspace(PROFILE_START, PROFILE_END, num=PROFILE_NUM)
    zp = np.zeros_like(xp)
    return xp, zp


@pytest.fixture
def prism():
    """
    The test prism: a right rectangular prism of finite extent in the strike direction.

    Unlike the cylinder, the polygon nodes describe this body *exactly*, so a 2.75D benchmark against it carries no
    discretisation error of its own.
    """
    return SimpleNamespace(
        x1=PRISM_X1,
        x2=PRISM_X2,
        z1=PRISM_Z1,
        z2=PRISM_Z2,
        y1=PRISM_Y1,
        y2=PRISM_Y2,
        nodes=get_rectangle_points(PRISM_X1, PRISM_X2, PRISM_Z1, PRISM_Z2),
    )


@pytest.fixture
def coarse_profile():
    """X and Z coordinates of a smaller set of observation points."""
    xp = np.linspace(PROFILE_START, PROFILE_END, num=COARSE_PROFILE_NUM)
    zp = np.zeros_like(xp)
    return xp, zp
