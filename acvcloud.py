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
    config = configparser.SafeConfigParser({"spacing":"1"})
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
    pdbParser = PDB.PDBParser(PERMISSIVE=1, QUIET=1)
    structure = pdbParser.get_structure("bio", pdbFile)

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

    # remove the edges of the grid (make it a sphere with a radius of the maximal linker length)
    Kd_grid = spatial.cKDTree(grid)
    sphereIdx = Kd_grid.query_ball_point(coords_AttachPos, linker)
    grid = grid[sphereIdx]
    return grid


def getPDBvolume(linker, element_target, coords_target, vdWRadius, grid, Kd_grid, spacing):
    """
    Compute a spacing fill model of the PDB structure based on atomic coordinates and VdW radii of the involved elements

    Parameters
    ----------
    linker : dye linker length in Angstrom
    element_target : list of elements in the target structure
    coords_AttachPos : list of atomic coordinates of the attachment position
    vdWRadius : dictionary of the VdW radius for each element
    grid : 2*linker+1 x 3 array of grid points

    Returns
    -------
    grid_PDBvol : m x 3 array of coordinates making up the PDB volume
    grid_notPDBvol : k x 3 array of grid points without those belonging to the spacing fill model
    """
    PDBvolIdx = []
    for i,elem in enumerate(element_target):
        PDBvolIdx.append(Kd_grid.query_ball_point(coords_target[i], vdWRadius[elem]+spacing/2))  # find all grid points within the VdW radii (+ half of the grid spacing) of the atoms belonging to the target
        # Note: adding half the grid spacing helps finding missing points due to the finite spacing of the grid

    PDBvolIdx = list(set([item for sublist in PDBvolIdx for item in sublist if item])) # get indices of grid points inside the PDB volume
    grid_PDBvol = grid[PDBvolIdx]
    notPDBvolIdx = list(set(range(len(grid))) - set(PDBvolIdx)) # "-" is the difference operator
    grid_notPDBvol = grid[notPDBvolIdx]
    return grid_PDBvol, grid_notPDBvol

def buildNeighborList(Kd_gridnotPDBvol, coords_AttachPos, grid_notPDBvol, nof_nodes):
    """
    Build list of neighbors of all grid points not belonging to the PDB volume

    Parameters
    ----------
    Kd_gridnotPDBvol :
    coords_AttachPos :
    grid_notPDBvol :
    nof_nodes :

    Returns:
    --------
    distDict :
    adjDict :
    """
    # compute neighbor list for grid points (nodes)
    distances_src, idx_src = Kd_gridnotPDBvol.query(coords_AttachPos, k=13**3, distance_upper_bound=6)  # return k neighbors from attachement point with a max. distance of 6 Angstrom
    distances, idx = Kd_gridnotPDBvol.query(grid_notPDBvol, k=7**3, distance_upper_bound=3.01) # return k neighbors from attachement point with a max. distance of 3.01 Angstrom
    distDict = {i: list(distances[i,np.isfinite(distances[i,:])]) for i in range(nof_nodes)}
    adjDict = {i: list(idx[i,np.isfinite(distances[i,:])]) for i in range(nof_nodes)}
    adjDict[nof_nodes] = idx_src
    distDict[nof_nodes] = distances_src
    return distDict, adjDict


def dijkstra(distDict, adjDict, src, nof_nodes, linker):
    """
    Find all grid points which have a path length from the attachment point that is shorter than or equal to the linker length using Dijkstra's algorithm

    Parameters
    ----------
    distDict : dictionary of distance lists for every grid point in grid_notPDBvol
    adjDict : dictionary of adjacency list
    src : source node
    N : number of nodes
    linker : dye linker length in Angstrom

    Returns
    -------
    nodeswithin : list of grid points having a path length shorter than or equal to the linker length
    """
    # initalize
    dist = [float('inf') for i in range(nof_nodes+1)]
    pq = [(0, src)]
    while pq:
        current_weight, u = heapq.heappop(pq)
        i = 0
        for v in adjDict[u]:
            weight = current_weight + distDict[u][i]
            if dist[v] > weight:
                dist[v] = weight
                heapq.heappush(pq, (weight, v))
            i += 1

    # select nodes with dist smaller than linker
    nodeswithin = []
    for i in range(nof_nodes):
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
    grid_AV : Accessible volume of dye
    grid_notPDBvol : k x 3 array of grid points without those belonging to the
    coords_AttachPos : list of atomic coordinates of the attachment position
    grid_CV : Contact volume of dye
    grid_FV : Free volume of dye (FV = AV - CV)
    """

    # parameters
    serialID = config.getint('dye parameters', 'serialID')
    linker = config.getint('dye parameters', 'linker')
    CVthickness = config.getfloat('dye parameters', 'CVthick')
    spacing = config.getfloat('grid', 'spacing')

    # PDB coordinates
    coords_target, element_target, coords_AttachPos = getPDBcoord(pdbFile, serialID)
    Kd_coords = spatial.cKDTree(coords_target)

    # grid construction
    grid = generate_grid(linker, spacing, coords_AttachPos)
    Kd_grid = spatial.cKDTree(grid)

    # compute PDB volume
    vdWRadius = {'H':1.1, 'C':1.7, 'O': 1.52, 'N':1.55, 'P':1.8}
    grid_PDBvol, grid_notPDBvol = getPDBvolume(linker, element_target, coords_target, vdWRadius, grid, Kd_grid, spacing)
    Kd_gridPDBvol = spatial.cKDTree(grid_PDBvol)
    Kd_gridnotPDBvol = spatial.cKDTree(grid_notPDBvol)
    nof_nodes = len(grid_notPDBvol)

    # compute nodes with path lengths shorter than linker length
    distDict, adjDict = buildNeighborList(Kd_gridnotPDBvol, coords_AttachPos, grid_notPDBvol, nof_nodes)
    nodeswithin = dijkstra(distDict, adjDict, nof_nodes, nof_nodes, linker)
    grid_AV = grid_notPDBvol[nodeswithin]
    Kd_grid_AV = spatial.cKDTree(grid_AV)

    # compute CV
    distMat_CV = spatial.cKDTree.sparse_distance_matrix(Kd_gridPDBvol, Kd_grid_AV, CVthickness)
    idx_CV = list(set([item[1] for item in distMat_CV.keys()])) # get grid indices within CV thickness from the PDB volume
    grid_CV = grid_AV[idx_CV]

    #from scipy.spatial import ConvexHull
    #hull = ConvexHull(grid_PDBvol)
    #surf = grid_PDBvol[hull.vertices]
    #writeXYZ('surf', surf, coords_AttachPos)


    # compute FV (difference of coordinates of AV that are not in CV)
    nrows, ncols = grid_AV.shape
    dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [grid_AV.dtype]}
    grid_FV = np.setdiff1d(grid_AV.view(dtype), grid_CV.view(dtype))
    grid_FV = grid_FV.view(grid_AV.dtype).reshape(-1, ncols)

    return coords_AttachPos, grid_PDBvol, grid_notPDBvol, grid_AV, grid_CV, grid_FV


def writeXYZ(name, coords, coords_AttachPos):
    """
    Write XYZ coordinates of accessible, contact and free volume

    Parameters
    ----------
    name : name of the coordinate file (without suffix)
    coords : n x 3 array of grid points
    coords_AttachPos : list of atomic coordinates of the attachment position
    """
    filename = '{}.xyz'.format(name)
    f = open(filename, 'w')
    f.write('attachmentCoords '+" ".join(map(str, coords_AttachPos))+'\n')
    for c in coords:
        f.write('D '+" ".join(map(str, c))+'\n')
    f.close()



if __name__ == "__main__":
    pdbFile, paramFile = parseCmd()
    config = getConfig(paramFile)
    coords_AttachPos, grid_PDBvol, grid_notPDBvol, grid_AV, grid_CV, grid_FV = runACV(config)
    writeXYZ('PDB', grid_PDBvol, coords_AttachPos)
    writeXYZ('AV', grid_AV, coords_AttachPos)
    writeXYZ('CV', grid_CV, coords_AttachPos)
    writeXYZ('FV', grid_FV, coords_AttachPos)
