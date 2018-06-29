#!/usr/bin/python3
# compute accessible-contact volume clouds for a given PDB get_structure
# F.Steffen, University of Zurich

import numpy as np
import re
import heapq
import argparse
import configparser
from Bio import PDB
from scipy import spatial
from itertools import groupby
from operator import itemgetter

def parseCmd():
    """
    parses the command line to get the input PDB file and the dye parameters
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
    # default parameters
    config = configparser.SafeConfigParser({"gammaCorr":"yes", "trajTimeInterval":"1", "deltat":"10"})
    config.readfp(open(paramFile))
    return config


def getPDBcoord(pdbFile, serialID):
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


def generate_grid(n, spacing, coords_AttachPos):
    row = 0
    grid = np.zeros(((2*n+1)**3,3))

    for i in range(2*n+1):
        for j in range(2*n+1):
            for k in range(2*n+1):
                grid[row,0] = i*spacing+coords_AttachPos[0]-n
                grid[row,1] = j*spacing+coords_AttachPos[1]-n
                grid[row,2] = k*spacing+coords_AttachPos[2]-n
                row += 1
    return grid


def getRNAvolume(n, element_target, coords_target, vdWRadius, grid, Kd_grid, CVthickness):
    index_RNAvol = []
    for i in range(len(element_target)):
        #print(element_target[i])

        index_RNAvol.append(Kd_grid.query_ball_point(coords_target[i], vdWRadius[element_target[i]]+CVthickness))
    index_RNAvol = list(set([item for sublist in index_RNAvol for item in sublist if item]))
    grid_RNAvol = grid[index_RNAvol]
    index_notRNAvol = list(set(range(len(grid))) - set(index_RNAvol)) # "-" is the difference operator
    grid_notRNAvol = grid[index_notRNAvol]

    tuple_gridindex = np.vstack(np.unravel_index(range(len(grid)), (2*n,2*n,2*n))).T
    tuple_RNAnotvol_gridindex = tuple_gridindex[index_notRNAvol,:]
    return grid_RNAvol, grid_notRNAvol, tuple_RNAnotvol_gridindex


def dijkstra(distanceDict,adjDict, src, N, n):
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
        if dist[i] < n:
            nodeswithin.append(i)

    return nodeswithin


def runACV(config):
    # parameters
    serialID = config.getint('dye parameters', 'serialID')
    n = config.getfloat('dye parameters', 'n')
    CVthickness = config.getfloat('dye parameters', 'CVthick')
    spacing = config.getfloat('grid', 'spacing')

    coords_target, element_target, coords_AttachPos = getPDBcoord(pdbFile, serialID)

    grid = generate_grid(n, spacing, coords_AttachPos)
    Kd_grid = spatial.cKDTree(grid)

    # kick out the edges of the grid (make it a sphere with a radius of the masimal linker length)
    index = Kd_grid.query_ball_point(coords_AttachPos, n)
    grid = grid[index]

    # construct KD-Trees
    Kd_coords = spatial.cKDTree(coords_target)
    Kd_grid = spatial.cKDTree(grid)

    vdWRadius = {"H":1.1, "C":1.7, "O": 1.52, "N":1.55, "P":1.8}
    grid_RNAvol, grid_notRNAvol, tuple_RNAnotvol_gridindex = getRNAvolume(n, element_target, coords_target, vdWRadius, grid, Kd_grid, CVthickness)


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

    return grid_final, grid_notRNAvol, grid


def writeXYZ(grid_final, grid_notRNAvol, grid):
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




if __name__ == "__main__":
    pdbFile, paramFile = parseCmd()
    config = getConfig(paramFile)
    runACV(config)
