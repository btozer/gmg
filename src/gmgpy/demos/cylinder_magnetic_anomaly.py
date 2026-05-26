"""
Analytical solution to calculate the magnetic anomaly 
caused by a 2D infinite cylinder.

Used to benchmark 2D forward modeling of magnetic anomalies
using the Talwani & Heirtzler (1964) algorithm.

**references**

Lowrie, W. Fundamentals of Geophysics pg. 330 eq. 5.54

Sleep & Fujita 1997. Pg. 223 eq. 6.45
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import argparse

def points(r, x1, y1, N):
    """
    Calculate the coordinates of points defining a circle.

    :param r: Radius of circle (m)
    :type r: float
    :param x1: X position of circle center
    :type x1: float
    :param y1: Y position of circle center
    :type y1: float
    :param N: Number of nodes defining the 'circle'
    :type N: int
    :return: 2D array (floats): C1:X coordinate; C2:Y coordinate
    :rtype: np.ndarray
    """
    angle = np.linspace(0., 2. * np.pi, N)
    output_points = np.zeros(shape=(len(angle), 2))

    for i in range(len(angle)):
        output_points[i, 0] = x1 + (r * np.cos(angle[i]))
        output_points[i, 1] = y1 + (r * np.sin(angle[i]))

    return output_points

def analytic_solution(x1, z, r, mc):
    """
    Calculate the magnetic anomaly along a profile.

    :param x1: X position of circle center
    :type x1: float
    :param z: Depth to center of Cylinder (m)
    :type z: float
    :param r: Radius of cylinder (m)
    :type r: float
    :param mc: Magnetization contrast of cylinder
    :type mc: float
    :return: The predicted anomaly from the infinite cylinder (nT)
    :rtype: tuple (np.ndarray, np.ndarray)
    """
    u = 1  # Magnetic permeability

    # Calculate points along profile
    x = np.linspace((x1 - r * 8), (x1 + r * 8), num=101)
    z2 = z**2

    # Calculate anomaly along the profile
    anomaly = np.zeros_like(x)
    for i in range(len(x)):
        # n_t = 0.5 * u * mc * (r**2) * ( (z2 - ((xp-x)**2)) / ((z2 + ((xp-x)**2))**2) )  ## Lowrie
        anomaly[i] = mc * (r**2) * ((((x[i] - x1)**2) - z2) / (2 * ((((x[i] - x1)**2) + z2)**2)))  ## Sleep & Fujita 1997

    return x, anomaly

def plot_results(x, anomaly, cylinder_points):
    """
    Plot the magnetic anomaly and cylinder points.

    :param x: X coordinates of the profile
    :type x: np.ndarray
    :param anomaly: Magnetic anomaly values
    :type anomaly: np.ndarray
    :param cylinder_points: Points defining the cylinder
    :type cylinder_points: np.ndarray
    """
    fig, axes = plt.subplots(nrows=2, ncols=1)

    # Set X limits
    plt.xlim(x.min() * 0.001, x.max() * 0.001)

    # Plot magnetic anomaly in the first axes
    axes[0].plot(x * 0.001, anomaly, c='r', zorder=1)
    axes[0].set_title('Magnetic Anomaly')
    axes[0].set_xlabel('Distance (km)')
    axes[0].set_ylabel('Anomaly (nT)')

    # Plot cylinder in the second axes
    axes[1].plot(cylinder_points[:, 0] * 0.001, cylinder_points[:, 1] * 0.001, c='b', zorder=1)
    axes[1].scatter(cylinder_points[:, 0] * 0.001, cylinder_points[:, 1] * 0.001, c='b', zorder=1)
    axes[1].invert_yaxis()
    axes[1].set_aspect('equal')
    axes[1].set_title('Cylinder Cross-Section')
    axes[1].set_xlabel('Distance (km)')
    axes[1].set_ylabel('Depth (km)')

    # Show figure
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Program parameters
    parser = argparse.ArgumentParser(
        description='Create a 2D circle cross-section through a cylinder and calculate the 2D magnetic anomaly',
        usage='#####################################################\n'
              '                %(prog)s -r -x -y -n -d                    \n'
              '                                                           \n'
              '                -r{Circle radius (m)}                      \n'
              '                -x{Circle center X location (m)}           \n'
              '                -y{Circle center Y location (m)}           \n'
              '                -n{Number of nodes defining circle}        \n'
              '                -k{Susceptibility contrast}                \n'
              '                -f{Earths magnetic field strength (nT)}    \n'
              '###########################################################\n',
        add_help=True)

    parser.add_argument('-r', help='Circle radius (m)', required=True, type=float)
    parser.add_argument('-x', help='Circle center X location (m)', required=True, type=float)
    parser.add_argument('-y', help='Circle center Y location (m)', required=True, type=float)
    parser.add_argument('-n', help='Number of nodes defining circle', required=True, type=int)
    parser.add_argument('-k', help='Susceptibility contrast', required=True, type=float)
    parser.add_argument('-f', help='Earths magnetic field strength (nT)', required=True, type=float)
    args = vars(parser.parse_args())

    # Calculate points
    cylinder_points = points(args['r'], args['x'], args['y'], args['n'])
    # Save points
    np.savetxt('mag_cylinder.txt', cylinder_points[::-1] * 0.001, delimiter=' ', fmt='%.5f %.5f')

    # Calculate magnetic anomaly
    mc = args['f'] * args['k']
    x, anomaly = analytic_solution(args['x'], args['y'], args['r'], mc)

    # Output calculated cylinder anomaly
    np.savetxt('mag_anom.xa', np.column_stack((x * 0.001, anomaly)), delimiter=' ', fmt='%.6f %.6f')

    # Plot results
    plot_results(x, anomaly, cylinder_points)
