import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import argparse

## CONSTANTS
SI2MGAL = 100000.0  ## Conversion factor from SI units to mGal: :math:`1\ m/s^2 = 10^5\ mGal`
G = 0.00000000006673  ## The gravitational constant in :math:`m^3 kg^{-1} s^{-1}`

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


def analytic_solution(r, x1, y1, rho_D):
    """
    :param r: float : RADIUS (m)
    :param x1: float : X position of circle
    :param y1: float : Y position of circle
    :param rho_D: float : Density contrast of cylinder

    :returns:
                1D array : x : X position along the profile
                1D array : anomaly : Calulated anomly produced by the cylinder
    """
    x = np.linspace(-25000., 25000., (50000. / 100))
    h = y1
    anomaly = np.zeros_like(x)
    for i in xrange(len(x)):
        # anomaly[i] = ( (4. * G * rho_D * np.pi * (r**3) * h) / (3. * ((x[i]**2 + h**2)**(3./2.))) ) * SI2MGAL
        anomaly[i] = ((2. * np.pi * G * (r ** 2) * rho_D) * (h / ((x[i] ** 2) + (h ** 2)))) * SI2MGAL
    x = x + x1
    return x*0.001, anomaly


########################################################################################################################
if __name__ == '__main__':

    ## PROGRAM PARAMETERS
    parser = argparse.ArgumentParser(description='reate a 2D circle cross-section through a cylinder and calculate the '
                                                 '2D gravity anomaly',
                                     usage='#####################################################\n'
                                           '                %(prog)s -r -x -y -n -d                    \n'
                                           '                                                           \n'
                                           '                -r{Circle radius (m)}                      \n'
                                           '                -x{Circle center X location (m)}           \n'
                                           '                -y{Circle center Y location (m)}           \n'
                                           '                -n{Number of nodes defining circle}        \n'
                                           '                -d{Density contrast (kg/m^3)}              \n'
                                           '###########################################################\n'
                                           ' \n',
                                     add_help=True)

    parser.add_argument('-r', help='Circle radius (m)', required=True, type=float)
    parser.add_argument('-x', help='Circle center X location (m)', required=True, type=float)
    parser.add_argument('-y', help='Circle center Y location (m)', required=True, type=float)
    parser.add_argument('-n', help='Number of nodes defining circle', required=True, type=int)
    parser.add_argument('-d', help='Density contrast (kg/m^3)', required=True, type=float)
    args = vars(parser.parse_args())

    ## CALCULATE POINTS
    points = points(args['r'], args['x'], args['y'], args['n'])
    ## SAVE POINTS
    np.savetxt('grav_cylinder.txt', points[::-1]*0.001, delimiter=' ', fmt='%.5f %.5f')


    ### CALCULATE THE GRAVITY ANOMALY
    x, anomaly = analytic_solution(args['r'], args['x'], args['y'], args['d'])
    ## SAVE POINTS
    np.savetxt('grav_anomaly.xa', zip(x, anomaly), delimiter=' ', fmt='%.5f %.5f')

    ## PLOT
    # fig = plt.figure()
    # plt.plot(x, anomaly, c='r', zorder=1)
    # plt.show()
########################################################################################################################