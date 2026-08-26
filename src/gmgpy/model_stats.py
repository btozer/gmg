"""
CALCULATE MODEL MISFIT STATISTICS  VALUES FOR OBSERVED VS CALCULATED VALUES
    1. RMS 
    2. CHI-SQUARE
"""

from scipy import interpolate
import numpy as np
from numpy import mean, sqrt, square


def rms(obs_x, obs_y, calc_x, calc_y):
    """CALCULATE RMS MISFIT BETWEEN OBSERVED AND CALCULATED VALUES"""

    # LIMIT OBSERVED VALUES TO THE X RANGE OF THE CALCULATED VALUES
    obs = np.column_stack((obs_x, obs_y))
    obs = obs[obs[:, 0] >= calc_x.min()]
    obs = obs[obs[:, 0] <= calc_x.max()]

    # INTERPOLATE THE VALUES OF THE CALCULATED ANOMALY AT THE OBSERVED SAMPLE POINTS
    f = interpolate.interp1d(calc_x, calc_y)  # CREATE INTERPOLATION FUNC FOR CALCULATED DATA
    intp_calc_y = f(obs[:, 0])  # INTERPOLATE CALC DATA ONTO OBS X POINTS

    # CALCULATE THE RESIDUALS
    res = (obs[:, 1] - intp_calc_y)
    residuals = np.column_stack((obs[:, 0], res))

    # CALCULATE THE RMS VALUE, ROUNDED TO 4 D.P.
    rms_misfit = round(sqrt(mean(square(res))), 4)

    return rms_misfit, residuals


def best_fit_dc(obs_x, obs_y, calc_x, calc_y):
    """
    CALCULATE THE CONSTANT (DC) SHIFT THAT ZEROES THE MEAN RESIDUAL WHEN ADDED TO THE CALCULATED
    CURVE, I.E. mean(obs - calc) AT THE OBSERVED SAMPLE POINTS (GM-SYS STYLE AUTO DC SHIFT).
    USES THE SAME CLIPPING AND INTERPOLATION AS rms()
    """

    # LIMIT OBSERVED VALUES TO THE X RANGE OF THE CALCULATED VALUES
    obs = np.column_stack((obs_x, obs_y))
    obs = obs[obs[:, 0] >= calc_x.min()]
    obs = obs[obs[:, 0] <= calc_x.max()]

    # NO OBSERVED POINTS IN RANGE MEANS THERE IS NOTHING TO FIT TO (AND mean([]) WOULD BE nan)
    if len(obs) == 0:
        return 0.0

    # INTERPOLATE THE VALUES OF THE CALCULATED ANOMALY AT THE OBSERVED SAMPLE POINTS
    f = interpolate.interp1d(calc_x, calc_y)

    return float(mean(obs[:, 1] - f(obs[:, 0])))


def chi_square(x, calc):
    """CALCULATE CHI-SQUARE MISFIT BETWEEN OBSERVED AND CALCULATED VALUES"""
    pass
