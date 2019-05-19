#!/usr/bin/env python

"""
PYTHON SCRIPT FOR CALCULATING HORIZONTAL DERIVATIVE OF A 1D ARRAY USING FINITE CENTRAL DIFFERENCE METHOD
WRITTEN BY B. TOZER UNIV. OXFORD 06.2016

NOTES::
INPUT DATA SHOULD BE XY SPACE DELIMITED ASCII TXT FILE
ORDER = DERIVATIVE ORDER i.e. 1=1st, 2=2nd etc
"""

# IMPORT MODULES
import numpy as np
from scipy import interpolate as ip
from scipy import signal
import matplotlib.pyplot as plt
import argparse

def read_data(raw_data):
    """
    READ IN INPUT DATA FROM FILE
    """

    data = np.genfromtxt(raw_data, delimiter=' ', dtype=(float))

    # First interpolate the data onto a even x spacing as required by the finite difference derivative approximation
    interpolate_func = ip.interp1d(data[:, 0], data[:, 1], kind='slinear')  #  create interpolation function

    x_min             = np.round(data[:, 0].min(), 0)+float(x_inc)
    x_max             = np.round(data[:, 0].max(), 0)-float(x_inc)
    x_interp_values   = np.arange(int(x_min), int(x_max)+1, float(x_inc))     #  x values for interpolated array
    y_interp_values = np.array(interpolate_func(x_interp_values))           #  y values for interpolated array
    input_interpolated     = np.column_stack((x_interp_values, y_interp_values))

    # IF THE DATA IS TO BE MEDIAN FILTERED BEFORE CALCULATING THE HORIZONTAL DERIVATIVE
    if args['s']:
        median_filter_points   = args['s']
        filtered_output        = np.zeros_like(input_interpolated)
        filtered_output[:, 0]  = input_interpolated[:, 0]
        filtered_output[:, 1]  = signal.medfilt(input_interpolated[:, 1], median_filter_points)
        input_interpolated = filtered_output

    return input_interpolated

def horizontal_derivative(input_interpolated, order, x_inc):
    """
    # CALCULATING HORIZONTAL DERIVATIVE OF A 1D ARRAY USING FINITE CENTRAL DIFFERENCE METHOD
    """

    # NOW CALCULATE THE DERIVATIVE USING THE INTERPOLATED (EVENLY SPACED) INPUT DATA
    N = len(input_interpolated)
    deriv = np.zeros(shape=((len(input_interpolated)-1), 2))
    print(input_interpolated)
    print("")

    for i in range(1, len(input_interpolated)-1):
        deriv[i, 0] = input_interpolated[i, 0]


        # CALC USING FORWARD METHOD
        deriv[i, 1] = abs((float(input_interpolated[i+1, 1])-float(input_interpolated[i, 1])))/(float(x_inc))

        # CALC USING CENTRAL DIFFERENCE
        # deriv[i, 1] = abs((float(input_interpolated[i+1, 1]) - float(input_interpolated[i-1, 1])))/(2.*float(x_inc))

    # SET FIRST AND LAST VALUES (WHICH ARE NOT CALCULATED BY FD METHOD) EQUAL TO SECOND AND SECOND-TO-LAST VALUES
    deriv[0, 0], deriv[0, 1] = deriv[1, 0]-1, deriv[1, 1]
    deriv[-1, 0], deriv[-1, 1] = deriv[-2, 0]+1, deriv[-2, 1]

    print(len(input_interpolated))
    print(len(deriv))

    # REPEAT IF ORDER IS > 1; e.g. second deriv order=2
    # while order > 1:
    #     deriv = horizontal_derivative(deriv, order-1, x_inc)
    #     order = order-1

    # RETURN CALCULATED DERIV
    return deriv

def plot(raw_data, input_interpolated, deriv, plot_name):
    """
    PLOT THE DATA AND DERIV
    """

    input_data = np.genfromtxt(raw_data, delimiter=' ', dtype=(float))

    fig = plt.figure(figsize=(11.69, 8.27), dpi=600)
    all_fontsize   = 4
    title_fontsize = 4

    # PLOT ON TOP AXIS
    ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

    # LABEL AXIS
    plt.xlabel('Distance (km)', fontsize=all_fontsize)
    plt.ylabel('Anomaly', fontsize=all_fontsize)

    # PLOT INPUT DATA
    plt.plot(input_data[:, 0], input_data[:, 1], c='blue', linewidth=0.2, zorder=1)
    plt.scatter(input_data[:, 0], input_data[:, 1], c='blue', s=0.5, zorder=2)

    # PLOT FILTERED INPUT DATA
    plt.plot(input_interpolated[:, 0], input_interpolated[:, 1], c='green', linewidth=0.2, zorder=1)
    plt.scatter(input_interpolated[:, 0], input_interpolated[:, 1], c='green', s=0.5, zorder=2)

    # PLOT DERIV
    # PLOT ON BOOTOM AXIS
    # ax = plt.subplot2grid((3,1), (2,0), rowspan=1, colspan=1)
    ax2 = ax.twinx()

    # LABEL AXIS
    plt.xlabel('Distance (km)', fontsize=all_fontsize)
    plt.ylabel('Derivative', fontsize=all_fontsize)

    # PLOT DERIV
    plt.plot(deriv[:, 0], deriv[:, 1], c='red', linewidth=0.2, zorder=1)
    plt.scatter(deriv[:, 0], deriv[:, 1], c='red', s=1, zorder=2)

    plt.show()
    exit(0)

    # WRITE OUT PDF FIG
    plt.savefig(plot_name, bbox_inches='tight', dpi='720')

if __name__=='__main__':

    # PROGRAM PARAMETERS
    parser = argparse.ArgumentParser(description='CALCULATING HORIZONTAL DERIVATIVE OF A 1D ARRAY USING FINITE CENTRAL '
                                                 'DIFFERENCE METHOD',
                                     usage='\n'
                                           '                %(prog)s -n -o -x optional{-p -w} \n'
                                           '\n'
                                           '                -n{Input file name} \n'
                                           '                -o{Derivative Order e.g. 1=1st, 2=2nd etc} \n'
                                           '                -x{X increment for interpolated data} \n'
                                           '                -s{Median smoothing filter value e.g. 5 = 5 point filter} \n'
                                           '                -p{Plot file name} \n'
                                           '                -w{Write file name} \n'
                                           '                -h{Print help} \n',
                                     add_help=True)

    parser.add_argument('-n', help='NAME OF INPUT FILE', required=True)
    parser.add_argument('-o', help='Derivative order', default=1, choices=(1, 2, 3, 4), type=int, required=True)
    parser.add_argument('-x', help='X INCREMENT FOR INTERPOLATED DATA', type=float, required=True)
    parser.add_argument('-s', help='USE A MEDIAN FILTER ON THE DATA', type=int, required=False)
    parser.add_argument('-p', help='PLOT FIGURE TO PDF: <name_of_file> e.g.-pderiv_plot.pdf', default='derv_plot.pdf',
                        required=False, type=str)
    parser.add_argument('-w', help='OUTPUT DERIV TXT FILE: <name_of_file> default=deriv.txt', default='deriv.txt',
                        required=False, type=str)

    args = vars(parser.parse_args())

    # READ IN DATA
    raw_data = args['n']
    order = args['o']
    x_inc = args['x']

    # RUN ALGORITHM
    input_derv_data = read_data(raw_data)
    deriv = horizontal_derivative(input_derv_data, order, x_inc)

    # SAVE OUTPUT AS .txt FILE
    if args['w']:
        output_name = str(args['w'])
        np.savetxt(output_name, deriv, delimiter=" ", fmt="%f %f")

    # # SAVE PLOT
    if args['p']:
        plot_name = str(args['p'])
        np.set_printoptions(suppress=True)
        np.set_printoptions(precision=4)
        # print(deriv.astype(float))

        plot(raw_data, input_derv_data, deriv, plot_name)
    # else:
    #     pass
