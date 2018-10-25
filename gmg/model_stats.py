"""
CALCULATE RMS AND CHI-SQUARED VALUES FOR OBSERVED VS CALCULATED VALUES
"""

from scipy import interpolate
import numpy as np
from numpy import mean, sqrt, square


def rms(obs_x, obs_y, calc_x, calc_y):
    """CALCULATE RMS MISFIT BETWEEN OBSERVED AND CALCULATED VALUES"""

    # INTERPOLATE THE VALUES OF THE CALCULATED ANOMALY AT THE OBSERVED SAMPLE POINTS
    f = interpolate.interp1d(calc_x, calc_y)  # CREATE INTERPOLATION FUNC FOR CALCULATED DATA
    intp_calc_y = f(obs_x)  # INTERPOLATE CALC DATA ONTO OBS X POINTS

    # %CALCULATE THE RESIDUALS
    res = (obs_y - intp_calc_y)
    residuals = np.column_stack((obs_x, res))

    # CALCULATE THE RMS VALUE, ROUNDED TO 2 D.P.
    rms_misfit = round(sqrt(mean(square(res))), 2)

    return rms_misfit, residuals


def chi_squared(x, calc):
    """CALCULATE CHI-SQUARED MISFIT BETWEEN OBSERVED AND CALCULATED VALUES"""
    pass
