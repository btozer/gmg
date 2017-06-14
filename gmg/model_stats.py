"""
Calculate rms and chi-squared values for observed vs calculated values.
x = the x points with observed data values
calc = array of x's vs calculated values
"""

from scipy import interpolate
import numpy as np
from numpy import mean, sqrt, square

def rms(obs_x, obs_y, calc_x, calc_y):
    """#% CALCULATE RMS MISFIT BETWEEN OBSERVED AND CALCULATED VALUES"""

    '# %INTERPOLATE THE VALUES OF THE CALCULATED ANOMALY AT THE OBSERVED SAMPLE POINTS'
    f = interpolate.interp1d(calc_x, calc_y)  # %create interpolation func for calculated data
    intp_calc_y = f(obs_x)  # %interpolate calc data onto obs X points

    '# %CALCULATE THE RESIDUALS'
    res = (obs_y - intp_calc_y)
    residuals = np.column_stack((obs_x, res))

    '# %CALCULATE THE RMS VALUE'
    rms = round(sqrt(mean(square(res))), 2)

    return rms, residuals

def chi_squared(x, calc):
    pass