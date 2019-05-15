#!/usr/bin/env python3

import numpy as np
import LabelLib as ll


def dist_mp(acv1, acv2):
    """
    Parameters
    ----------

    Returns
    -------

    Examples
    --------

    >>> Rmp()

    """
    R_mp = np.sqrt(sum((acv1.mp - acv2.mp)**2))
    return R_mp


def mean_dist_DA_ll(acv1, acv2, n_dist=10**6):
    """
    Compute a randomly subsampled donor-acceptor distance distribution for the two accessible volume clouds

    Parameters
    ----------
    av1 : n x 3 array of xyz coordinates
    av1 : m x 3 array of xyz coordinates
    nof_dist : integer number of distances to calculate

    Returns
    -------
    1D array of distances between av1 and av2
    """
    R_DA = ll.meanDistance(acv1.ll_Grid3D, acv2.ll_Grid3D, n_dist)
    return R_DA


def dists_DA(acv1, acv2, n_dist=10**6):
    """
    Compute a randomly subsampled donor-acceptor distance distribution for the two accessible volume clouds

    Parameters
    ----------
    av1 : n x 3 array of xyz coordinates
    av1 : m x 3 array of xyz coordinates
    nof_dist : integer number of distances to calculate

    Returns
    -------
    1D array of distances between av1 and av2
    """
    idx1 = np.random.randint(0, len(acv1.cloud_xyzqt), n_dist)
    idx2 = np.random.randint(0, len(acv2.cloud_xyzqt), n_dist)
    r = np.sqrt(np.sum((acv1.cloud_xyzqt[idx1, 0:3] - acv2.cloud_xyzqt[idx2, 0:3])**2, axis=1))
    w = acv1.cloud_xyzqt[idx1, 3] * acv2.cloud_xyzqt[idx2, 3]
    w = w / w.sum()
    R_DA = np.vstack((r, w)).T
    return R_DA


def mean_dist_DA(acv1, acv2, R_DA=None, n_dist=10**6):
    """
    Parameters
    ----------
    structure

    Returns
    -------

    Examples
    --------

    >>> acv1 =
    >>> mean_dist_DA(acv1, acv2, n_dist=10**6)

    """
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        print('> calculating R_DAs')
    mean_R_DA = np.dot(R_DA[:, 0], R_DA[:, 1]) / R_DA[:, 1].sum()
    return mean_R_DA


def std_dist_DA(acv1, acv2, R_DA=None, n_dist=10**6):
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        print('> calculating R_DAs')
    std_R_DA = np.sqrt(np.dot(R_DA[:, 0]**2, R_DA[:, 1]) - np.dot(R_DA[:, 0], R_DA[:, 1])**2)
    return std_R_DA


def FRET_DA(acv1, acv2, R_DA=None, R0=54, n_dist=10**6):
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        print('> calculating R_DAs')
    E_DA = np.vstack((R0**6 / (R0**6 + R_DA[:, 0]**6), R_DA[:, 1])).T
    return E_DA


def mean_FRET_DA(acv1, acv2, E_DA=None, R_DA=None, R0=54, n_dist=10**6):
    if E_DA is None:
        E_DA = FRET_DA(acv1, acv2, R_DA, R0, n_dist)
        print('> calculating E_DAs')
    mean_E_DA = np.dot(E_DA[:, 0], E_DA[:, 1]) / E_DA[:, 1].sum()
    return mean_E_DA


def mean_FRET_DA_ll(acv1, acv2, R0=54, n_dist=10**6):
    return ll.meanEfficiency(acv1.ll_Grid3D, acv2.ll_Grid3D, R0, n_dist)


def mean_dist_DA_fromFRET(acv1, acv2, mean_E_DA=None, E_DA=None, R_DA=None, R0=54, n_dist=10**6):
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
        mean_E_DA = mean_FRET_DA(acv1, acv2, E_DA, R_DA, R0, n_dist)
        print('> calculating mean_E_DA')
    mean_R_DA_E = R0 * (1 / mean_E_DA - 1)**(1 / 6)
    return mean_R_DA_E


def dists_DA_ll(acv1, acv2, n_dist=10**6):
    # check weight calculation in LabelLib
    dists_DA = ll.sampleDistanceDistInv(acv1.ll_Grid3D, acv2.ll_Grid3D, n_dist)
    return dists_DA


def dist_attach(attach1_xyz, attach2_xyz):
    R_attach = np.sqrt(sum((attach1_xyz - attach2_xyz)**2))
    return R_attach
