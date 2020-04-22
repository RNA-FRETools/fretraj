#!/usr/bin/env python3

import numpy as np
try:
    import LabelLib as ll
except ModuleNotFoundError:
    _LabelLib_found = False
else:
    _LabelLib_found = True


def dist_mp(acv1, acv2):
    """
    Compute distance between mean dye positions of two accessible contact volumes

    Parameters
    ----------
    acv1 : fretraj.cloud.ACV
    acv2 : fretraj.cloud.ACV

    Returns
    -------
    R_mp : float
           distance between mean dye positions of acv1 and acv2

    Examples
    --------

    >>> avobj.dist_mp(acv1, acv2)

    """
    R_mp = np.sqrt(sum((acv1.mp - acv2.mp)**2))
    return R_mp


def mean_dist_DA_ll(acv1, acv2, n_dist=10**6):
    """
    Mean donor-acceptor distance between acv1 and acv2

    Parameters
    ----------
    acv1 : fretraj.cloud.ACV
    acv2 : fretraj.cloud.ACV
    n_dist : int, optional
             integer number of distances to calculate

    Returns
    -------
    mean_R_DA : float
    """
    mean_R_DA = ll.meanDistance(acv1.ll_Grid3D, acv2.ll_Grid3D, n_dist)
    return mean_R_DA


def dists_DA(acv1, acv2, n_dist=10**6, return_weights=True):
    """
    Compute a randomly subsampled donor-acceptor distance distribution for the two accessible volume clouds

    Parameters
    ----------
    acv1 : fretraj.cloud.ACV
    acv2 : fretraj.cloud.ACV
    n_dist : int, optional
             integer number of distances to calculate

    Returns
    -------
    mean_R_DA : ndarray
                1D array of length n_dist containing distances between acv1 and acv2
    """
    if return_weights:
        idx1 = np.random.choice(len(acv1.cloud_xyzqt), n_dist)
        idx2 = np.random.choice(len(acv2.cloud_xyzqt), n_dist)
        R_DA = np.sqrt(np.sum((acv1.cloud_xyzqt[idx1, 0:3] - acv2.cloud_xyzqt[idx2, 0:3])**2, axis=1))
        w = acv1.cloud_xyzqt[idx1, 3] * acv2.cloud_xyzqt[idx2, 3]
        w = w / w.sum()
        R_DA = np.vstack((R_DA, w)).T
    else:
        w1 = acv1.cloud_xyzqt[:, 3] / sum(acv1.cloud_xyzqt[:, 3])
        w2 = acv2.cloud_xyzqt[:, 3] / sum(acv2.cloud_xyzqt[:, 3])
        idx1 = np.random.choice(len(acv1.cloud_xyzqt), n_dist, p=w1)
        idx2 = np.random.choice(len(acv2.cloud_xyzqt), n_dist, p=w2)
        R_DA = np.sqrt(np.sum((acv1.cloud_xyzqt[idx1, 0:3] - acv2.cloud_xyzqt[idx2, 0:3])**2, axis=1))
    return R_DA


def mean_dist_DA(acv1, acv2, R_DA=None, n_dist=10**6, verbose=False):
    """
    Parameters
    ----------
    structure
    n_dist : int, optional

    Returns
    -------

    Examples
    --------

    >>> acv1 =
    >>> mean_dist_DA(acv1, acv2, n_dist=10**6)

    """
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        if verbose:
            print('> calculating R_DAs')
    if R_DA.ndim > 1:
        mean_R_DA = np.dot(R_DA[:, 0], R_DA[:, 1]) / R_DA[:, 1].sum()
    else:
        np.mean(R_DA)
    return mean_R_DA


def std_dist_DA(acv1, acv2, R_DA=None, n_dist=10**6, verbose=False):
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        if verbose:
            print('> calculating R_DAs')
    if R_DA.ndim > 1:
        std_R_DA = np.sqrt(np.dot(R_DA[:, 0]**2, R_DA[:, 1]) - np.dot(R_DA[:, 0], R_DA[:, 1])**2)
    else:
        std_R_DA = np.std(R_DA)
    return std_R_DA


def FRET_DA(acv1, acv2, R_DA=None, R0=54, n_dist=10**6, verbose=False):
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        if verbose:
            print('> calculating R_DAs')
    if R_DA.ndim > 1:
        E_DA = np.vstack((R0**6 / (R0**6 + R_DA[:, 0]**6), R_DA[:, 1])).T
    else:
        E_DA = R0**6 / (R0**6 + R_DA**6)
    return E_DA


def mean_FRET_DA(acv1, acv2, E_DA=None, R_DA=None, R0=54, n_dist=10**6, verbose=False):
    if E_DA is None:
        E_DA = FRET_DA(acv1, acv2, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating E_DAs')
    if E_DA.ndim > 1:
        mean_E_DA = np.dot(E_DA[:, 0], E_DA[:, 1]) / E_DA[:, 1].sum()
    else:
        mean_E_DA = np.mean(E_DA)
    return mean_E_DA

def std_FRET_DA(acv1, acv2, E_DA=None, R_DA=None, R0=54, n_dist=10**6, verbose=False):
    if E_DA is None:
        E_DA = FRET_DA(acv1, acv2, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating E_DAs')
    if E_DA.ndim > 1:
        std_E_DA = np.sqrt(np.dot(E_DA[:, 0]**2, E_DA[:, 1]) - np.dot(E_DA[:, 0], E_DA[:, 1])**2)
    else:
        std_E_DA = np.std(E_DA)
    return std_E_DA


def mean_FRET_DA_ll(acv1, acv2, R0=54, n_dist=10**6):
    return ll.meanEfficiency(acv1.ll_Grid3D, acv2.ll_Grid3D, R0, n_dist)


def mean_dist_DA_fromFRET(acv1, acv2, mean_E_DA=None, E_DA=None, R_DA=None, R0=54, n_dist=10**6, verbose=False):
    """

    Examples
    --------

    # with known mean_E_DA (fastest)
    >>> mean_dist_DA_fromFRET(acv1, acv2, mean_E_DA=mymean_E_DA)

    # with known E_DA (mean_E_DA is calculated on-the-fly)
    >>> mean_dist_DA_fromFRET(acv1, acv2, E_DA=myE_DA)

    # with known R_DA (E_DA and mean E_DA are calculated on-the-fly)
    >>> mean_dist_DA_fromFRET(acv1, acv2, R_DA=myR_DA)


    """
    if mean_E_DA is None:
        mean_E_DA = mean_FRET_DA(acv1, acv2, E_DA, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating mean_E_DA')
    mean_R_DA_E = R0 * (1 / mean_E_DA - 1)**(1 / 6)
    return mean_R_DA_E

def std_dist_DA_fromFRET(acv1, acv2, mean_E_DA=None, sigma_E_DA=None, R_DA=None, R0=54, n_dist=10**6, verbose=False):
    if (mean_E_DA is None) or (sigma_E_DA is None):
        mean_E_DA = mean_FRET_DA(acv1, acv2, E_DA, R_DA, R0, n_dist, verbose)
        sigma_E_DA = std_FRET_DA(acv1, acv2, E_DA, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating mean_E_DA and sigma_E_DA')
    std_R_DA_E = np.sqrt((R0 / (6 * ( 1 / mean_E_DA - 1)**(5 / 6) * mean_E_DA**2 ))**2 * sigma_E_DA**2)  # error propagation: sqrt( (df/dx)^2 * sigma_x^2 ) 
    return std_R_DA_E


def dists_DA_ll(acv1, acv2, n_dist=10**6, return_weights=True):
    R_DA = np.array(ll.sampleDistanceDistInv(acv1.ll_Grid3D, acv2.ll_Grid3D, n_dist))
    if return_weights:
        w = np.full((1, n_dist), 1 / n_dist)
        R_DA = np.vstack((R_DA, w)).T
    return R_DA


def dist_attach(attach1_xyz, attach2_xyz):
    R_attach = np.sqrt(sum((attach1_xyz - attach2_xyz)**2))
    return R_attach
