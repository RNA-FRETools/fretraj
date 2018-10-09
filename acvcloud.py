#!/usr/bin/env python3
# compute accessible-contact volume clouds for a given PDB get_structure
# F.Steffen, University of Zurich

import numpy as np
import re, sys
import heapq
import argparse
import configparser
from Bio import PDB
from scipy import spatial
from itertools import groupby
from operator import itemgetter

def parseCmd():
    """
    Parse the command line to get the input PDB file and the dye parameters

    Returns
    -------
    pdbFile : string with input PDB filename
    paramFile : string with parameter filename
    """
    version = 1.0

    parser = argparse.ArgumentParser(description='compute accessible-contact clouds for a given PDB structure')
    parser.add_argument('--version', action='version', version='%(prog)s ' + str(version))
    parser.add_argument('-i', help='Input PDB structure (.pdb)', required=True)
    parser.add_argument('-p', help='Parameter file (.dat)', required=True)
    args = parser.parse_args()
    pdbFile = args.i
    paramFile = args.p
    return pdbFile, paramFile


def getConfig(paramFile):
    """
    Parse the parameter file to get the configuration settings for acvCloud

    Parameters
    ----------
    paramFile : string with parameter filename

    Returns
    -------
    config : SafeConfigParser object with the configuration settings
    """
    # default parameters
    config = configparser.SafeConfigParser({"gammaCorr":"yes", "trajTimeInterval":"1", "deltat":"10"})
    config.readfp(open(paramFile))
    return config


def getPDBcoord(pdbFile, serialID):
    """
    Get coordinates of the input pdb file and the dye attachement position

    Parameters
    ----------
    pdbFile : string with input PDB filename
    paramFile : string with parameter filename

    Returns
    -------
    coords_target : n x 3 array of atomic coordinates of the target structure
    element_target : list of elements in the target structure
    coords_AttachPos : list of atomic coordinates of the attachment position
    """
    # read pdb files
    p = PDB.PDBParser(PERMISSIVE=1, QUIET=1)
    structure = p.get_structure("rna", pdbFile)

    # get coordinates of structure
    coords_target = []
    element_target = []
    atoms = structure.get_atoms()
    for atom in atoms:
        coords_target.append(atom.get_coord())
        element_target.append(re.findall("[A-Z]", atom.get_name())[0]) # suppose that the first uppercase alphabetic character is the element (HO is H whereas OH is O)
    coords_target = np.array(coords_target, dtype=np.dtype('float64'))
    element_target = np.array(element_target)

    # get coordinates of attachment position (specified by serialID)
    selection = PDB.Selection.unfold_entities(structure, target_level='A') # unfold to the atom level ("A")
    serial_numbers = [atom.serial_number for atom in selection]
    selection_dict = dict(zip(serial_numbers, selection))
    coords_AttachPos = list(selection_dict[serialID].get_coord().astype(np.dtype('float64')))
    return coords_target, element_target, coords_AttachPos


def generate_grid(linker, spacing, coords_AttachPos):
    """
    Calculate regular grid around attachment position

    Parameters
    ----------
    linker : dye linker length in Angstrom
    spacing : spacing between grid points
    coords_AttachPos : list of atomic coordinates of the attachment position

    Returns
    -------
    grid : 2*linker+1 x 3 array of grid points

    """
    row = 0
    grid = np.zeros(((2*linker+1)**3,3))

    for i in range(2*linker+1):
        for j in range(2*linker+1):
            for k in range(2*linker+1):
                grid[row,0] = i*spacing+coords_AttachPos[0]-linker
                grid[row,1] = j*spacing+coords_AttachPos[1]-linker
                grid[row,2] = k*spacing+coords_AttachPos[2]-linker
                row += 1
    return grid


def getRNAvolume(linker, element_target, coords_target, vdWRadius, grid, Kd_grid, CVthickness):
    """
    Compute a spacing fill model of the RNA based on atomic coordinates and VdW radii of the involved elements

    Parameters
    ----------
    linker : dye linker length in Angstrom
    element_target : list of elements in the target structure
    coords_AttachPos : list of atomic coordinates of the attachment position
    vdWRadius : dictionary of the VdW radius for each element
    grid : 2*linker+1 x 3 array of grid points

    Returns
    -------
    grid_RNAvol : m x 3 array of coordinates making up the RNA volume
    grid_notRNAvol : k x 3 array of grid points without those belonging to the spacing fill model
    """
    index_RNAvol = []
    for i in range(len(element_target)):
        #print(element_target[i])

        index_RNAvol.append(Kd_grid.query_ball_point(coords_target[i], vdWRadius[element_target[i]]+CVthickness))
    index_RNAvol = list(set([item for sublist in index_RNAvol for item in sublist if item]))
    grid_RNAvol = grid[index_RNAvol]
    index_notRNAvol = list(set(range(len(grid))) - set(index_RNAvol)) # "-" is the difference operator
    grid_notRNAvol = grid[index_notRNAvol]

    tuple_gridindex = np.vstack(np.unravel_index(range(len(grid)), (2*linker,2*linker,2*linker))).T
    #tuple_RNAnotvol_gridindex = tuple_gridindex[index_notRNAvol,:]
    return grid_RNAvol, grid_notRNAvol


def dijkstra(distanceDict,adjDict, src, N, linker):
    """
    Find all grid points which have a path length from the attachment point that is shorter than or equal to the linker length using Dijkstra's algorithm

    Parameters
    ----------
    distanceDict : dictionary of distance lists for every grid point in grid_notRNAvol : k x 3 array of grid points without those belonging to the
    adjDict : dictionary of adjacency list
    src : source node
    N : number of nodes
    linker : dye linker length in Angstrom

    Returns
    -------
    nodeswithin : list of grid points having a path length shorter than or equal to the linker length
    """
    # initalize
    dist = [float('inf') for i in range(N+1)]

    pq = [(0, src)]

    while pq:
        current_weight, u = heapq.heappop(pq)


        i = 0
        for v in adjDict[u]:
            weight = current_weight + distanceDict[u][i]

            if dist[v] > weight:
                dist[v] = weight
                heapq.heappush(pq, (weight, v))
            i += 1

    # select nodes with dist smaller than n
    nodeswithin = []
    for i in range(N):
        if dist[i] <= linker:
            nodeswithin.append(i)

    return nodeswithin


def runACV(config):
    """
    Run the ACV calculation

    Parameters
    -----------
    config : SafeConfigParser object with the configuration settings

    Returns
    -------
    grid_final : Accessible volume of dye
    grid_notRNAvol : k x 3 array of grid points without those belonging to the
    grid : 2*linker+1 x 3 array of grid points
    coords_AttachPos : list of atomic coordinates of the attachment position
    grid_CV : Contact volume of dye
    grid_FV : Free volume of dye (FV = AV - CV)
    """

    # parameters
    serialID = config.getint('dye parameters', 'serialID')
    linker = config.getint('dye parameters', 'linker')
    CVthickness = config.getfloat('dye parameters', 'CVthick')
    spacing = config.getfloat('grid', 'spacing')

    coords_target, element_target, coords_AttachPos = getPDBcoord(pdbFile, serialID)

    grid = generate_grid(n, spacing, coords_AttachPos)
    Kd_grid = spatial.cKDTree(grid)

    # kick out the edges of the grid (make it a sphere with a radius of the maximal linker length)
    index = Kd_grid.query_ball_point(coords_AttachPos, n)
    grid = grid[index]

    # construct KD-Trees
    Kd_coords = spatial.cKDTree(coords_target)
    Kd_grid = spatial.cKDTree(grid)

    vdWRadius = {"H":1.1, "C":1.7, "O": 1.52, "N":1.55, "P":1.8}
    grid_RNAvol, grid_notRNAvol = getRNAvolume(n, element_target, coords_target, vdWRadius, grid, Kd_grid, CVthickness)


    #grid_notRNAvol = np.vstack((grid_notRNAvol,coords_AttachPos))
    Kd_gridnotRNAvol = spatial.cKDTree(grid_notRNAvol)


    N = len(grid_notRNAvol)

    distances_src, ind_src = Kd_gridnotRNAvol.query(coords_AttachPos, k=13**3, distance_upper_bound=6)

    distances, ind = Kd_gridnotRNAvol.query(grid_notRNAvol, k=7**3, distance_upper_bound=3.01)

    distanceDict = {i: list(distances[i,np.isfinite(distances[i,:])]) for i in range(N)}
    adjDict = {i: list(ind[i,np.isfinite(distances[i,:])]) for i in range(N)}

    adjDict[N] = ind_src
    distanceDict[N] = distances_src

    #print distances[distances[:,0]==10]
    #print [(k, list(list(zip(*g))[1])) for k, g in groupby(distances.keys(), itemgetter(0))]

    #edges = [key[0][1] for key in distances[5000,:].items()]
    #N = 6
    #distances = {(1,3):5, (1,3):2, (1,4):2, (2,3):2, (2,4):4, (3,4):1}

    nodeswithin = dijkstra(distanceDict, adjDict, N, N, n)

    grid_final = grid_notRNAvol[nodeswithin]

    # CV
     # construct KD-Trees
    Kd_gridRNAvol = spatial.cKDTree(grid_RNAvol)
    Kd_grid_final = spatial.cKDTree(grid_final)

    # compute distance matrix
    distMat_CV = spatial.cKDTree.sparse_distance_matrix(Kd_gridRNAvol, Kd_grid_final, CVthickness)

    # get grid indices within CV thickness from the RNA volume
    index_CV = list(set([item[1] for item in distMat_CV.keys()]))

    # get the coordinates of the CV
    grid_CV = grid_final[index_CV]


    # get coordinates of the FV (difference of coordinates of AV that are not in CV)
    nrows, ncols = grid_final.shape
    dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [grid_final.dtype]}
    grid_FV = np.setdiff1d(grid_final.view(dtype), grid_CV.view(dtype))
    grid_FV = grid_FV.view(grid_final.dtype).reshape(-1, ncols)


    return grid_final, grid_notRNAvol, grid, coords_AttachPos, grid_CV, grid_FV


def writeXYZ(grid_final, grid_notRNAvol, grid, coords_AttachPos, grid_CV, grid_FV):
    """
    Write XYZ coordinates of accessible, contact and free volume

    Parameters
    ----------
    grid_final : Accessible volume of dye
    grid_notRNAvol : k x 3 array of grid points without those belonging to the
    grid : 2*linker+1 x 3 array of grid points
    coords_AttachPos : list of atomic coordinates of the attachment position
    grid_CV : Contact volume of dye
    grid_FV : Free volume of dye (FV = AV - CV)
    """
    filenameRNA = 'RNA.xyz'
    f = open(filenameRNA,'w')
    f.write('attachmentCoords '+" ".join(map(str, coords_AttachPos))+'\n')
    for i in range(len(grid)):
        f.write('D '+" ".join(map(str, grid[i]))+'\n')
    f.close()

    filenameRNA = 'RNA_not.xyz'
    f = open(filenameRNA,'w')
    f.write('attachmentCoords '+" ".join(map(str, coords_AttachPos))+'\n')
    for i in range(len(grid_notRNAvol)):
        f.write('D '+" ".join(map(str, grid_notRNAvol[i]))+'\n')
    f.close()

    filenameRNA = 'RNA_final.xyz'
    f = open(filenameRNA,'w')
    f.write('attachmentCoords '+" ".join(map(str, coords_AttachPos))+'\n')
    for i in range(len(grid_final)):
        f.write('D '+" ".join(map(str, grid_final[i]))+'\n')
    f.close()


    # contact volume
    filenameCV = 'CV.xyz'
    f = open(filenameCV,'w')
    f.write('attachmentCoords '+" ".join(map(str, coords_AttachPos))+'\n')
    for i in range(len(grid_CV)):
        f.write('D '+" ".join(map(str, grid_CV[i]))+'\n')
    f.close()

    # Free volume
    filenameFV = 'FV.xyz'
    f = open(filenameFV,'w')
    f.write('attachmentCoords '+" ".join(map(str, coords_AttachPos))+'\n')
    for i in range(grid_FV.shape[0]):
        f.write('D '+" ".join(map(str, grid_FV[i]))+'\n')
    f.close()



if __name__ == "__main__":
    pdbFile, paramFile = parseCmd()
    config = getConfig(paramFile)
    grid_final, grid_notRNAvol, grid, coords_AttachPos, grid_CV, grid_FV = runACV(config)
    writeXYZ(grid_final, grid_notRNAvol, grid, coords_AttachPos, grid_CV, grid_FV)
