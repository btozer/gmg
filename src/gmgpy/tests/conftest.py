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

# OBSERVATION PROFILE (SI UNITS)
PROFILE_START = 0.0
PROFILE_END = 4000.0
PROFILE_NUM = 4001


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


@pytest.fixture
def profile():
    """X and Z coordinates of the observation points."""
    xp = np.linspace(PROFILE_START, PROFILE_END, num=PROFILE_NUM)
    zp = np.zeros_like(xp)
    return xp, zp
