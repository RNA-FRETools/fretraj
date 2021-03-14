#!/usr/bin/env python3

import numpy as np
try:
    import LabelLib as ll
except ModuleNotFoundError:
    _LabelLib_found = False
else:
    _LabelLib_found = True


def dist_mp(acv1, acv2):
    """Compute the distance between mean dye positions R_MP of two accessible
    contact volumes

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs

    Returns
    -------
    float
        distance between the mean dye positions of acv1 and acv2
    """
    R_mp = np.sqrt(sum((acv1.mp - acv2.mp)**2))
    return R_mp


def dists_DA(acv1, acv2, n_dist=10**6, return_weights=True):
    """Compute a randomly subsampled donor-acceptor distance distribution
    for the two accessible volume clouds

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    n_dist : int, default=10E6
        number of distances to calculate
    return_weights : bool, default=True
        additional vector with weights

    Returns
    -------
    ndarray
        ndarray of shape 1 x n_dist or 2 x n_dist (if weights are included)
        containing distances between acv1 and acv2
    """
    if return_weights:
        idx1 = np.random.choice(len(acv1.cloud_xyzqt), n_dist)
        idx2 = np.random.choice(len(acv2.cloud_xyzqt), n_dist)
        R_DA = np.sqrt(np.sum((acv1.cloud_xyzqt[idx1, 0:3]
                               - acv2.cloud_xyzqt[idx2, 0:3])**2, axis=1))
        w = acv1.cloud_xyzqt[idx1, 3] * acv2.cloud_xyzqt[idx2, 3]
        w = w / w.sum()
        R_DA = np.vstack((R_DA, w)).T
    else:
        w1 = acv1.cloud_xyzqt[:, 3] / sum(acv1.cloud_xyzqt[:, 3])
        w2 = acv2.cloud_xyzqt[:, 3] / sum(acv2.cloud_xyzqt[:, 3])
        idx1 = np.random.choice(len(acv1.cloud_xyzqt), n_dist, p=w1)
        idx2 = np.random.choice(len(acv2.cloud_xyzqt), n_dist, p=w2)
        R_DA = np.sqrt(np.sum((acv1.cloud_xyzqt[idx1, 0:3]
                               - acv2.cloud_xyzqt[idx2, 0:3])**2, axis=1))
    return R_DA


def dists_DA_ll(acv1, acv2, n_dist=10**6, return_weights=True):
    """Compute a randomly subsampled donor-acceptor distance distribution
    for the two accessible volume clouds

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    n_dist : int, default=10E6
        number of distances to calculate
    include_weights : bool, default=True
        additional vector with weights

    Returns
    -------
    ndarray
        ndarray of shape 1 x n_dist or 2 x n_dist (if weights are included)
        containing distances between acv1 and acv2
    """
    R_DA = np.array(ll.sampleDistanceDistInv(acv1.ll_Grid3D, acv2.ll_Grid3D,
                    n_dist))
    if return_weights:
        w = np.full((1, n_dist), 1 / n_dist)
        R_DA = np.vstack((R_DA, w)).T
    return R_DA


def mean_dist_DA(acv1, acv2, R_DA=None, n_dist=10**6, verbose=False):
    """Calculate the average donor-acceptor distance <R_DA>
    between donor and acceptor

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    R_DA : ndarray
        ndarray of length n_dist containing distances between acv1 and acv2
    n_dist : int, default=10E6
        number of distances to calculate
    verbose: bool

    Returns
    -------
    float
        average distance between donor and acceptor
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


def mean_dist_DA_ll(acv1, acv2, n_dist=10**6):
    """Calculate the average donor-acceptor distance <R_DA> between
    acv1 and acv2 using LabelLib

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    n_dist : int, default=10E6
        number of distances to calculate

    Returns
    -------
    float
        average distance between donor and acceptor
    """
    mean_R_DA = ll.meanDistance(acv1.ll_Grid3D, acv2.ll_Grid3D, n_dist)
    return mean_R_DA


def std_dist_DA(acv1, acv2, R_DA=None, n_dist=10**6, verbose=False):
    """Calculate the standard deviation of the donor acceptor distances sigma_R_DA

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    R_DA : ndarray
        ndarray of length n_dist containing distances between acv1 and acv2
    n_dist : int, default=10E6
        number of distances to calculate
    verbose : bool, default=False
        be verbose about calculation steps

    Returns
    -------
    float
        standard deviation of donor acceptor distances
    """
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        if verbose:
            print('> calculating R_DAs')
    if R_DA.ndim > 1:
        sigma_R_DA = np.sqrt(np.dot(R_DA[:, 0]**2, R_DA[:, 1])
                             - np.dot(R_DA[:, 0], R_DA[:, 1])**2)
    else:
        sigma_R_DA = np.std(R_DA)
    return sigma_R_DA


def FRET_DA(acv1, acv2, R_DA=None, R0=54, n_dist=10**6, verbose=False):
    """Calculate the FRET efficiencies E for each donor acceptor distance

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    R_DA : ndarray, default=None
        ndarray of length n_dist containing distances between acv1 and acv2
        (if None calculate on the fly)
    R0 : float, default=54
        Förster radius in Angstrom
    n_dist : int, default=10E6
        number of distances to calculate
    verbose : bool, default=False
        be verbose about calculation steps

    Returns
    -------
    ndarray
        ndarray of shape 1 x n_dist or 2 x n_dist (if weights are included)
        containing FRET efficiencies between acv1 and acv2

    Notes
    -----
    The transfer efficiency *E* is given by

    .. math::  E = \\frac{R_0^6}{R_0^6+R_\\text{DA}^6}

    """
    if R_DA is None:
        R_DA = dists_DA(acv1, acv2, n_dist)
        if verbose:
            print('> calculating R_DAs')
    if R_DA.ndim > 1:
        E_DA = np.vstack((R0**6 / (R0**6 + R_DA[:, 0]**6), R_DA[:, 1])).T
    else:
        E_DA = R0**6 / (R0**6 + R_DA**6)
    return E_DA


def mean_FRET_DA(acv1, acv2, E_DA=None, R_DA=None, R0=54, n_dist=10**6,
                 verbose=False):
    """Calculate the average FRET efficiency <E> between donor and acceptor cloud

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    E_DA : ndarray, default=None
        ndarray of length n_dist containing FRET efficiencies
        between acv1 and acv2 (if None calculate on the fly)
    R_DA : ndarray, default=None
        ndarray of length n_dist containing distances
        between acv1 and acv2 (if None calculate on the fly)
    R0 : float, default=54
        Förster radius in Angstrom
    n_dist : int, default=10E6
        number of distances to calculate
    verbose : bool, default=False
        be verbose about calculation steps

    Returns
    -------
    float
        average FRET efficiency between donor and acceptor cloud
    """
    if E_DA is None:
        E_DA = FRET_DA(acv1, acv2, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating E_DAs')
    if E_DA.ndim > 1:
        mean_E_DA = np.dot(E_DA[:, 0], E_DA[:, 1]) / E_DA[:, 1].sum()
    else:
        mean_E_DA = np.mean(E_DA)
    return mean_E_DA


def std_FRET_DA(acv1, acv2, E_DA=None, R_DA=None, R0=54, n_dist=10**6,
                verbose=False):
    """Calculate the standard deviation of the FRET efficiencies sigma_E

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    E_DA : ndarray, default=None
        ndarray of length n_dist containing FRET efficiencies
        between acv1 and acv2 (if None calculate on the fly)
    R_DA : ndarray, default=None
        ndarray of length n_dist containing distances
        between acv1 and acv2 (if None calculate on the fly)
    R0 : float, default=54
        Förster radius in Angstrom
    n_dist : int, default=10E6
        number of distances to calculate
    verbose : bool, default=False
        be verbose about calculation steps

    Returns
    -------
    float
        standard deviation of the FRET efficiencies
    """
    if E_DA is None:
        E_DA = FRET_DA(acv1, acv2, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating E_DAs')
    if E_DA.ndim > 1:
        sigma_E_DA = np.sqrt(np.dot(E_DA[:, 0]**2, E_DA[:, 1])
                             - np.dot(E_DA[:, 0], E_DA[:, 1])**2)
    else:
        sigma_E_DA = np.std(E_DA)
    return sigma_E_DA


def mean_FRET_DA_ll(acv1, acv2, R0=54, n_dist=10**6):
    """Calculate the average FRET efficiency <E> between donor and acceptor cloud
    using LabelLib

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    R0 : float, default=54
        Förster radius in Angstrom
    n_dist : int, default=10E6
        number of distances to calculate

    Returns
    -------
    float
        average FRET efficiency between donor and acceptor cloud
    """
    return ll.meanEfficiency(acv1.ll_Grid3D, acv2.ll_Grid3D, R0, n_dist)


def mean_dist_DA_fromFRET(acv1, acv2, mean_E_DA=None, E_DA=None, R_DA=None,
                          R0=54, n_dist=10**6, verbose=False):
    """Calculate the FRET averaged donor acceptor distance <R_DA>_E

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    mean_E_DA : float, default=None
        average FRET efficiency between donor and acceptor cloud
        (if None calculate on the fly)
    E_DA : ndarray, default=None
        ndarray of length n_dist containing FRET efficiencies
        between acv1 and acv2 (if None calculate on the fly)
    R_DA : ndarray, default=None
        ndarray of length n_dist containing distances
        between acv1 and acv2 (if None calculate on the fly)
    R0 : float, default=54
        Förster radius in Angstrom
    n_dist : int, default=10E6
        number of distances to calculate
    verbose : bool, default=False
        be verbose about calculation steps

    Returns
    -------
    float
        FRET averaged donor acceptor distance

    Examples
    --------

    with known :math:`\\langle E_\\text{DA}\\rangle` (fastest)

    >>> mean_dist_DA_fromFRET(acv1, acv2, mean_E_DA=mymean_E_DA)

    with known :math:`E_\\text{DA}`, i.e. :math:`\\langle E_\\text{DA}\\rangle` is calculated on-the-fly

    >>> mean_dist_DA_fromFRET(acv1, acv2, E_DA=myE_DA)

    with known :math:`R_\\text{DA}`, i.e. :math:`E_\\text{DA}` and
    :math:`\\langle E_\\text{DA}\\rangle` are calculated on-the-fly

    >>> mean_dist_DA_fromFRET(acv1, acv2, R_DA=myR_DA)


    """
    if mean_E_DA is None:
        mean_E_DA = mean_FRET_DA(acv1, acv2, E_DA, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating mean_E_DA')
    mean_R_DA_E = R0 * (1 / mean_E_DA - 1)**(1 / 6)
    return mean_R_DA_E


def std_dist_DA_fromFRET(acv1, acv2, mean_E_DA=None, sigma_E_DA=None,
                         E_DA=None, R_DA=None, R0=54, n_dist=10**6,
                         verbose=False):
    """Calculate the standard deviation of the FRET averaged donor acceptor
    distances sigma_R_DA_E by error propagation

    Parameters
    ----------
    acv1, acv2 : fretraj.cloud.ACV
        donor and acceptor ACVs
    mean_E_DA : float, default=None
        average FRET efficiency between donor and acceptor cloud
        (if None calculate on the fly)
    sigma_E_DA : float, default=None
        standard deviation of the FRET efficiencies
    E_DA : ndarray, default=None
        ndarray of length n_dist containing FRET efficiencies
        between acv1 and acv2 (if None calculate on the fly)
    R_DA : ndarray, default=None
        ndarray of length n_dist containing distances
        between acv1 and acv2 (if None calculate on the fly)
    R0 : float, default=54
         Förster radius in Angstrom
    n_dist : int, default=10E6
        number of distances to calculate
    verbose : bool, default=False
        be verbose about calculation steps

    Returns
    -------
    float
        standard deviation of the FRET averaged donor acceptor distances
    """
    if (mean_E_DA is None) or (sigma_E_DA is None):
        mean_E_DA = mean_FRET_DA(acv1, acv2, E_DA, R_DA, R0, n_dist, verbose)
        sigma_E_DA = std_FRET_DA(acv1, acv2, E_DA, R_DA, R0, n_dist, verbose)
        if verbose:
            print('> calculating mean_E_DA and sigma_E_DA')
    # error propagation: sqrt( (df/dx)^2 * sigma_x^2 )
    sigma_R_DA_E = np.sqrt((R0/(6 * (1/mean_E_DA-1)**(5/6) * mean_E_DA**2))**2
                           * sigma_E_DA**2)
    return sigma_R_DA_E


def dist_attach(attach1_xyz, attach2_xyz):
    """Calculate the distance between the attachment sites

    Parameters
    ----------
    attach1_xyz, attach2_xyz : ndarray
        one-dimensional array of x-,y-,z-coordinates of the attachment points

    Returns
    -------
    float
        distance between the attachment sites
    """
    R_attach = np.sqrt(sum((attach1_xyz - attach2_xyz)**2))
    return R_attach


def R_DAE_to_Rmp(mean_R_DA_E):
    """Conversion function for <R_DA,E> to R_mp derived in Kalinin, *Nat. Methods* (2012)

    Parameters
    ----------
    mean_R_DA_E : float
        FRET averaged donor acceptor distance

    Returns
    -------
    float
        distance between the mean dye positions

    Notes
    -----
    The conversion function between :math:`\\langle R_\\text{DA}\\rangle_E` and :math:`R_\\text{mp}` is a third-order
    polynomial:

    .. math::  R_\\text{mp} = 1.109\\cdot 10^{-5}\\cdot\\langle R_\\text{DA}\\rangle_E^3 - 7.286\\cdot 10^{-3}\\cdot
        \\langle R_\\text{DA}\\rangle_E^2 + 1.979\\cdot\\langle R_\\text{DA}\\rangle_E - 34.345

    References
    ----------

    .. [3] Kalinin, S. et al. "A toolkit and benchmark study for FRET-restrained high-precision \
        structural modeling", *Nat. Methods* **9**, 1218–1225 (2012).
    """
    R_mp = 1.109*10**-5*mean_R_DA_E**3 - 7.286*10**-3*mean_R_DA_E**2
    + 1.979*mean_R_DA_E - 34.345
    return R_mp
