"""
Calculate the magnetic anomaly caused by a 2D infinite cylinder.

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
    :param r: float : Radius of circle (m)
    :param x1: float : X position of circle center
    :param y1: float : Y position of circle center
    :param N: int : Number of nodes defining the 'circle'

    :returns: 2D array (floats): C1:X coordinate; C2:Y coordinate

    """
    angle = np.linspace(0., 2. * np.pi, N)
    output_points = np.zeros(shape=(len(angle), 2))

    for i in range(len(angle)):
        output_points[i, 0] = x1 + (r * np.cos(angle[i]))
        output_points[i, 1] = y1 + (r * np.sin(angle[i]))

    return output_points

def analytic_solution(x1, z, r, mc):
    """
    :param x1: float : X position of circle center
    :param z: Depth to center of Cylinder (m)
    :param r: Radius of cylinder (m)
    :param mc: Magnetization contrast of cylinder
    :return: The predicted anomaly from the infinite cylinder (nT)
    """

    ## MAGNETIC PERMEABILITY
    # u = (4*m.pi)*1E-7
    u = 1

    ## CALC POINTS ALONG PROFILE
    x = np.linspace((x1-r*8), (x1+r*8), num=101)
    print(x)
    ## CALCULATE ANOMALY ALONG THE PROFILE
    z2 = z**2

    anomaly = np.zeros_like(x)
    for i in range(len(x)):
        # n_t = 0.5 * u * mc * (r**2) * ( (z2 - ((xp-x)**2)) / ((z2 + ((xp-x)**2))**2) )  ## Lowrie
        anomaly[i] = mc * (r**2) * ( (((x[i]-x1)**2) - z2) / (2 * ((((x[i]-x1)**2) + z2 )**2)) )  ## S&F97
    print(anomaly)
    ## RETURN POINTS
    return x, anomaly

if __name__ == "__main__":

    ## PROGRAM PARAMETERS
    parser = argparse.ArgumentParser(
        description='Create a 2D circle cross-section'
                    'through a cylinder and calculate the'
                    '2D magnetic anomaly',
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

    parser.add_argument('-r', help='Circle radius (m)', 
                        required=True, type=float)
    parser.add_argument('-x', help='Circle center X location (m)', 
                        required=True, type=float)
    parser.add_argument('-y', help='Circle center Y location (m)', 
                        required=True, type=float)
    parser.add_argument('-n', help='Number of nodes defining circle', 
                        required=True, type=int)
    parser.add_argument('-k', help='Susceptibility contrast', 
                        required=True, type=float)
    parser.add_argument('-f', help='Earths magnetic field strength (nT)', 
                        required=True, type=float)
    args = vars(parser.parse_args())

    ## CALCULATE POINTS
    cylinder_points = points(args['r'], args['x'], args['y'], args['n'])
    ## SAVE POINTS
    np.savetxt('mag_cylinder.txt', cylinder_points[::-1]*0.001, delimiter=' ', 
               fmt='%.5f %.5f')

    ## CALCULATE MAGNETIC ANOMALY
    mc = args['f'] * args['k']
    x, anomaly = analytic_solution(args['x'], args['y'], args['r'], mc)

    ## OUTPUT CALCULATED CYLINDER ANOMALY
    np.savetxt('mag_anom.xa', np.column_stack((x*0.001, anomaly)), 
               delimiter=' ', fmt='%.6f %.6f')
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
