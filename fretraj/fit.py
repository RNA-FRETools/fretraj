#!/usr/bin/env python3

import numpy as np
import scipy as sp


def gauss(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


def gauss_mixt(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)


def do_fit(fun, x, y, startParam, bounds=None, sigma=None):
    if bounds is None:
        params, cov = sp.optimize.curve_fit(fun, x, y, startParam)
    else:
        if sigma is not None:
            # e.g. sigma=1/(x+1)
            params, cov = sp.optimize.curve_fit(fun, x, y, startParam, bounds=bounds, sigma=sigma)
        else:
            params, cov = sp.optimize.curve_fit(fun, x, y, startParam, bounds=bounds)

    sigma_paramError = np.sqrt(np.diag(cov))
    var = fun.__code__.co_varnames[1:]
    return params, sigma_paramError, var


def binMidpoint(bins):
    binMid = (bins[1:] + bins[:-1]) / 2  # for making len(x)==len(y)
    halfbinWidth = (bins[1] - bins[0]) / 2
    return binMid, halfbinWidth
