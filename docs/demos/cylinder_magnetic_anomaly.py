"""
Calculate the magnetic anomaly cuased by an infinite cylinder.

**references**

Lowrie, W. Fundamental of Geophysics pg. 330 eq. 5.54
Sleep & Fujita 1997. Pg. 223 eq. 6.45

"""

import numpy as np
import math as m
import argparse

def points(r, x1, y1, N):
    """
    :param r: float : RADIUS (m)
    :param x1: float : X position of circle
    :param y1: float : Y position of circle
    :param N: int : Number of nodes defining the 'circle'

    :returns: 2D array (floats): C1:X-position; C2:Gravity_anom

    """
    angle = np.linspace(0., 2. * np.pi, N)
    output_points = np.zeros(shape=(len(angle), 2))

    for i in xrange(len(angle)):
        output_points[i, 0] = x1 + (r * np.cos(angle[i]))
        output_points[i, 1] = y1 + (r * np.sin(angle[i]))

    return output_points

def n_t(x, z, r, mc):
    """
    :param x: List of X points where the predicted anomaly will be calculated
    :param z: Depth to center of Cylinder (m)
    :param r: Radius of cylinder (m)
    :param mc: Magnetization contrast of cylinder
    :return: The prdicted anomly from the infinite cylinder (nT)
    """

    ## CALC POINTS ALONG PROFILE
    xp = np.linspace(-25000., 25000., (50000. / 100))

    ## MAGNETIC PERMEABILITY
    # u = (4*m.pi)*1E-7
    u = 1

    ## CALCULATE ANOMALY ALONG THE PROFILE
    z2 = z**2
    # n_t = 0.5 * u * mc * (r**2) * ( (z2 - ((xp-x)**2)) / ((z2 + ((xp-x)**2))**2) )  ## Lowrie
    n_t = mc * (r**2) * ( (((xp-x)**2) - z2) / (2 * ((((xp-x)**2) + z2 )**2)) )  ## S&F97

    ## REUTURN POINTS
    return np.column_stack((xp * 0.001, n_t))

if __name__ == "__main__":

    ## PROGRAM PARAMETERS
    parser = argparse.ArgumentParser(description='Create a 2D circle cross-section through a cylinder and calculate the '
                                                 '2D magnietc anomaly',
                                     usage='#####################################################\n'
                                           '                %(prog)s -r -x -y -n -d                    \n'
                                           '                                                           \n'
                                           '                -r{Circle radius (m)}                      \n'
                                           '                -x{Circle center X location (m)}           \n'
                                           '                -y{Circle center Y location (m)}           \n'
                                           '                -n{Number of nodes defining circle}        \n'
                                           '                -k{Susceptibility contras}                 \n'
                                           '                -f{Earths magnetic field strength (nT)}    \n'                                           
                                           '###########################################################\n'
                                           ' \n',
                                     add_help=True)

    parser.add_argument('-r', help='Circle radius (m)', required=True, type=float)
    parser.add_argument('-x', help='Circle center X location (m)', required=True, type=float)
    parser.add_argument('-y', help='Circle center Y location (m)', required=True, type=float)
    parser.add_argument('-n', help='Number of nodes defining circle', required=True, type=int)
    parser.add_argument('-k', help='Susceptibility contrast', required=True, type=float)
    parser.add_argument('-f', help='Earths magnetic field strength (nT)', required=True, type=float)
    args = vars(parser.parse_args())


    # SET INPUTS
    x = 5000
    z = 2000
    r = 1000

    ## CALCULATE POINTS
    points = points(args['r'], args['x'], args['y'], args['n'])
    ## SAVE POINTS
    np.savetxt('mag_cylinder.txt', points[::-1]*0.001, delimiter=' ', fmt='%.5f %.5f')


    # CALCUALTE ANOMALY
    k = args['k']
    f = args['f']
    mc = f * k
    mag_anomaly = n_t(args['x'], args['y'], args['r'], mc)

    ## OUTPUT CALCULATED CYLINDER ANOMALY
    np.savetxt('mag_anom.xa', mag_anomaly, fmt='%.6f %.6f', delimiter=' ')
