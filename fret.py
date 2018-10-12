#!/usr/bin/env python3
# compute FRET from two accessible-contact volume clouds
# F.Steffen, University of Zurich

import numpy as np
import re, sys
import argparse
import configparser
from scipy import spatial

def parseCmd():
    """
    Parse the command line to get the input ACV clouds

    Returns
    -------
    acv_don : string with input donor ACV filename
    acv_acc : string with input donor ACV filename
    R0 : float specifying the Förster radius of the FRET pair
    """
    version = 1.0

    parser = argparse.ArgumentParser(description='calculate FRET from accessible-contact clouds')
    parser.add_argument('--version', action='version', version='%(prog)s ' + str(version))
    parser.add_argument('-p', help='ACV parameter file (.dat)', required=True)
    args = parser.parse_args()
    paramFile = args.p
    return paramFile

def getConfig(paramFile):
    """
    Parse the parameter file to get the configuration settings for acvCloud

    Parameters
    ----------
    paramFile : string with parameter filename

    Returns
    -------
    config : SafeConfigParser object
    """
    # default parameters
    config = configparser.SafeConfigParser({'R0':'54'})
    config.readfp(open(paramFile))
    return config

class ACV:
    def __init__(self, config, id):
        self.id = id
        self.loadCoords(config)
        self.mean_pos()

    def loadCoords(self, config):
        """
        Load coordinates of accessible volume cloud

        Parameters
        ----------
        config : SafeConfigParser object
        """
        self.filename = config.get('clouds', self.id)
        try:
            data = np.loadtxt(self.filename, skiprows=1, usecols=[1,2,3,4])
        except:
            data = np.loadtxt(self.filename, skiprows=1, usecols=[1,2,3])
        self.coords = data[:,0:3]
        if data.shape[1]>3:
            self.weights = data[:,3]
        else:
            self.weights = np.array([1.0]*self.coords.shape[0])

    def mean_pos(self):
        """
        Compute mean dye position
        """
        #self.mp = np.mean(self.coords, 0)
        self.mp = np.dot(self.weights,self.coords)/self.weights.sum()


def dist_mp(av1, av2):
    """
    Compute distance between mean dye positions of the two accessible volume clouds

    Parameters
    ----------
    av1 : n x 3 array of xyz coordinates
    av1 : m x 3 array of xyz coordinates

    Returns
    -------
    float of distance between av1.mp and av2.mp
    """
    return np.sqrt(sum((av1.mp-av2.mp)**2));


def dists_DA(av1, av2, nof_dist=10**6):
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
    idx_don = np.random.randint(0,av1.coords.shape[0], nof_dist)
    idx_acc = np.random.randint(0,av2.coords.shape[0], nof_dist)
    r = np.sqrt(np.sum((av1.coords[idx_don]-av2.coords[idx_acc])**2, axis=1))
    w = av1.weights[idx_don]*av2.weights[idx_acc]
    w = w/w.sum()
    return np.vstack((r,w)).T


def avg_dist_DA(R_DA):
    """
    Compute the (weighted) average donor-acceptor distance between the two accessible volume clouds

    Parameters
    ----------
    R_DA : k x 2 array of donor-acceptor distances and associated weights

    Returns
    -------
    float of average distance between av1 and av2
    """
    r = R_DA[:,0]
    w = R_DA[:,1]
    return np.dot(w,r)/w.sum()


def std_dist_DA(R_DA):
    """
    Compute the (weighted) width (= standard deviation) of the donor-acceptor distance distribution for the two accessible volume clouds

    Parameters
    ----------
    R_DA : k x 2 array of donor-acceptor distances and associated weights

    Returns
    -------
    float of width of donor-acceptor distribution for av1 and av2
    """
    r = R_DA[:,0]
    w = R_DA[:,1]
    return np.sqrt(np.dot(w,r**2) - np.dot(w,r)**2)


def mean_fret(R_DA, R0):
    """
    Compute the (weighted) mean FRET efficiency for the two accessible volume clouds

    Parameters
    ----------
    R_DA : k x 2 array of donor-acceptor distances and associated weights
    R0 : float of the Förster radius

    Returns
    -------
    float of mean FRET efficiency for av1 and av2
    """
    r = R_DA[:,0]
    w = R_DA[:,1]
    E = R0**6/(R0**6+r**6)
    return np.dot(w,E)/w.sum()


def mean_fret_dist(E_avg, R0):
    """
    Compute FRET-averaged distance between the two accessible volume clouds

    Parameters
    ----------
    E_avg : float of mean FRET efficiency for av1 and av2
    R0 : float of the Förster radius

    Returns
    -------
    float of FRET-averaged distance
    """
    return R0*(1/E_avg-1)**(1/6)


def writeOutput(Rmp, R_DA_avg, R_DA_std, E_avg, R_DA_E):
    print('ACV FRET calculation')
    print('----------------')
    print('Rmp = {:.1f}'.format(Rmp))
    print('<R_DA> = {:.1f}'.format(R_DA_avg))
    print('<sigma_DA> = {:.1f}'.format(R_DA_std))
    print('<E> = {:.2f}'.format(E_avg))
    print('<R_DA>_E = {:.1f}'.format(R_DA_E))
    print('----------------')


def fretObs(av1, av2, R0):
    Rmp = dist_mp(av1, av2)
    R_DA = dists_DA(av1, av2, 10**6)
    R_DA_avg = avg_dist_DA(R_DA)
    R_DA_std = std_dist_DA(R_DA)
    E_avg = mean_fret(R_DA, R0)
    R_DA_E = mean_fret_dist(E_avg, R0)
    writeOutput(Rmp, R_DA_avg, R_DA_std, E_avg, R_DA_E)


if __name__ == "__main__":
    paramFile = parseCmd()
    config = getConfig(paramFile)
    R0 = config.getfloat('förster', 'R0')
    av1 = ACV(config, 'av_don')
    av2 = ACV(config, 'av_acc')
    fretObs(av1, av2, R0)
