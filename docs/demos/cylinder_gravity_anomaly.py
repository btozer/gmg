"""
Calculate the gravity anomaly caused by a 2D infinite cylinder.

**references**

Analytical solution for buried cylinder: Garland (1965) Pg. 70 

"""
import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import argparse

## CONSTANTS
SI2MGAL = 100000.0  ## Conversion factor from SI units to mGal: :math:`1\ m/s^2 = 10^5\ mGal`
G = 0.00000000006673  ## The gravitational constant in :math:`m^3 kg^{-1} s^{-1}`

def points(r, x, z, N):
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

def analytic_solution(r, x1, z, rho_D):
    """
    :param r: float : RADIUS (m)
    :param x1: float : X position of circle center
    :param z: float : Z position of circle center
    :param rho_D: float : Density contrast of cylinder

    :returns:
                1D array : x : X position along the profile
                1D array : anomaly : Calculated anomaly
    """
    x = np.linspace((x1-r*8), (x1+r*8), num=101)
    anomaly = np.zeros_like(x)
    for i in range(len(x)):
        anomaly[i] = ((2. * np.pi * G * (r ** 2) * rho_D) * \
                      (z / (((x[i]-x1) ** 2) + (z ** 2)))) * SI2MGAL
    return x, anomaly

###############################################################################
if __name__ == '__main__':

    ## PROGRAM PARAMETERS
    parser = argparse.ArgumentParser(
        description='create a 2D circle cross-section '
                    'through a cylinder and calculate the '
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

    parser.add_argument('-r', help='Circle radius (m)', 
                        required=True, type=float)
    parser.add_argument('-x', help='Circle center X location (m)', 
                        required=True, type=float)
    parser.add_argument('-z', help='Circle center Z location (m)', 
                        required=True, type=float)
    parser.add_argument('-n', help='Number of nodes defining circle', 
                        required=True, type=int)
    parser.add_argument('-d', help='Density contrast (kg/m^3)', 
                        required=True, type=float)
    args = vars(parser.parse_args())

    ## ------------------------------------------------------------------------
    ## GENERATE DISCRETE POINTS DEFINING THE CYLINDER
    cylinder_points = points(args['r'], args['x'], args['z'], args['n'])
    
    ## SAVE POINTS TO DISC
    np.savetxt('grav_cylinder.txt', cylinder_points[::-1]*0.001, 
               delimiter=' ', fmt='%.5f %.5f')
    ## ------------------------------------------------------------------------

    ## ------------------------------------------------------------------------
    ## CALCULATE THE GRAVITY ANOMALY
    x, anomaly = analytic_solution(args['r'], args['x'], args['z'], args['d'])
    
    ## SAVE GRAVITY ANOMALY POINTS
    np.savetxt('grav_anomaly.xa', np.column_stack((x*0.001, anomaly)), 
               delimiter=' ', fmt='%.5f %.5f')
    ## ------------------------------------------------------------------------

    ## ------------------------------------------------------------------------
    ## PLOT USING km AS SCALE
    ## CREATE FIG
    fig, axes = plt.subplots(nrows=2, ncols=1)
    
    ## SET X LIMITS
    plt.xlim(x.min()*0.001, x.max()*0.001)

    ## PLOT GRAVITY ANOMALY IN FIRST AXES
    axes[0].plot(x*0.001, anomaly, c='r', zorder=1)
    
    ## PLOT CYLINDER IN SECOND AXES
    axes[1].plot(cylinder_points[:, 0]*0.001, 
                 cylinder_points[:, 1]*0.001, c='b', zorder=1)
    axes[1].scatter(cylinder_points[:, 0]*0.001, 
                    cylinder_points[:, 1]*0.001, c='b', zorder=1)
    axes[1].invert_yaxis() 
    axes[1].set_aspect('equal')
    
    ## SHOW FIGURE
    plt.tight_layout()
    plt.show()
    ## ------------------------------------------------------------------------
###############################################################################