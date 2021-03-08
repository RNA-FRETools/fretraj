#!/usr/bin/env python3
import numpy as np
import numba as nb
import heapq
import packaging.version

# Note: Reflected lists are deprecated in numba >=0.51. Typed lists are available from numba >=0.45.
# The heapq used for the priority queue in the Dijkstra algorithm does not accept typed list though (numba issue #4926).
# Thus, the Numba version needs to be pinned. Numba <=0.50 issues scheduled deprecation warnings
# but does run the Dijkstra implementation in nopython mode with Python 3.7 and 3.8 (not 3.9, numba issue #6345).
# Numba 0.50 requires numpy<1.20 (numba issue #6041)

_hasTypedList = packaging.version.parse(nb.__version__) >= packaging.version.parse('0.45.0')

_dist_list = np.sqrt(np.array([1, 2, 3, 5, 6]))  # reduction to 74 neighbors


class Grid3D:
    """Class object holding a 3D grid

    Parameters
    ----------
    mol_xyzr : ndarray
               array of x-,y-,z-coordinates and VdW radii with a shape [n_atoms, 4]
    attach_xyz : ndarray
                one-dimensional array of x-,y-,z-coordinates of the attachment point
                (corresponds to the center of the grid)
    linker_length : float
                    length of the dye linker in Angstrom
    linker_width : float
                   diameter of the dye linker in Angstrom
    dye_radii : ndarray([3,1])
                array of dye radii in Angstrom with shape [3,1]
    grid_spacing : float
        spacing between grid points (in A)

    Notes
    -----
    Attributes of the LabelLib class
        - discStep : float  (the grid spacing)
        - originXYZ : numpy.array  (x-/y-/z-coordinate of the grid origin)
        - shape : numpy.array      (number of grid points in x-/y-/z-direction)
        - grid : numpy.array       (flattened list of grid point values)
    """

    def __init__(self, mol_xyzr, attach_xyz, linker_length, linker_width, dye_radii, grid_spacing, simulation_type):
        self.discStep = grid_spacing
        self.attach_xyz = attach_xyz
        self.grid_3d, self.originXYZ, self.originAdj = self.make_grid(attach_xyz, linker_length, grid_spacing)
        self.shape = np.array(self.grid_3d.shape)
        self.halfCubeLength = min(self.shape) * self.discStep * 0.5
        self.grid_3d = self.block_molecule(mol_xyzr, 0.5 * linker_width)
        maxR = 3 * grid_spacing
        maxR_source = linker_width + grid_spacing
        self.grid_3d = self.dijkstra_init(maxR, maxR_source)
        self.grid_3d = self.setAboveTreshold(self.grid_3d, linker_length, -4)
        self.grid_3d = self.setAboveTreshold(self.grid_3d, 0, 1)
        if simulation_type == 'AV1':
            dye_radii = [dye_radii[0]]
        self.grid_3d = self.excludeConcentricSpheres(mol_xyzr, dye_radii, 1)
        self.grid = self.grid_3d.flatten(order='F')

    @staticmethod
    @nb.jit(forceobj=True)
    def make_grid(attach_xyz, linker_length, grid_spacing):
        """
        Build a 3D grid around the attachment point

        Parameters
        ----------
        attach_xyz : ndarray
                     one-dimensional array of x-,y-,z-coordinates of the attachment point
                     (corresponds to the center of the grid)
        linker_length : float
                        length of the dye linker in Angstrom
        grid_spacing : float
            spacing between grid points (in A)

        Returns
        -------
        grid_3d : ndarray
                  3-dimensional array of grid points with a shape of 2*ll_padRound+1
        xyz_min : ndarray
        """
        ll_pad = linker_length + 3 * grid_spacing  # pad the linker_length
        ll_padRound = np.ceil(ll_pad / grid_spacing) * grid_spacing  # round to next higher grid value
        xyz_min = attach_xyz - ll_padRound
        originAdj = xyz_min - 0.5 * grid_spacing
        gridptsPerEdge = 2 * int(ll_padRound / grid_spacing + 0.5)  # + 1
        grid_3d = np.full([gridptsPerEdge] * 3, np.inf)
        for i in range(gridptsPerEdge):
            for j in range(gridptsPerEdge):
                for k in range(gridptsPerEdge):
                    ijk = np.array([i, j, k])
                    maxRSq = int(linker_length**2 / grid_spacing**2 + 0.5)
                    ijk0 = (attach_xyz - xyz_min) / grid_spacing
                    dSq = np.sum((ijk - ijk0)**2)
                    if dSq > maxRSq:
                        grid_3d[(i, j, k)] = -1
        return grid_3d, xyz_min, originAdj

    @staticmethod
    def _xyz2idx(xyz, originAdj, grid_spacing, decimals=6):
        """Get the ijk grid indices for a set of xyz values

        Parameters
        ----------
        xyz :
        originAdj :
        grid_spacing : float
            spacing between grid points (in A)

        Note
        ----
        Round to n-decimal places (here: 4) before casting to integer
        """
        return np.round(((xyz - originAdj) / grid_spacing), decimals).astype(np.int64)

    def block_molecule(self, mol_xyzr, extraClash):
        """Block the grid points which are within the VdW radius of any atom of the biomolecule plus an extra clash radius

        Parameters
        ----------
        mol_xyzr : ndarray
            array of x-,y-,z-coordinates and VdW radii with a shape [n_atoms, 4]
        extraClash : float
            clash radius added to the VdW radius

        Returns
        -------
        ndarray
            3-dimensional array of grid points with a shape of 2*adjL+1
        """
        maxVdW_extraClash = np.max(mol_xyzr[:, 3]) + extraClash
        if _hasTypedList:
            neighbor_list = nb.typed.List()
            [neighbor_list.append(n) for n in self.sortedNeighborIdx(maxVdW_extraClash)]
        else:
            neighbor_list = self.sortedNeighborIdx(maxVdW_extraClash)
        ijk_atom = self._xyz2idx(mol_xyzr[:, 0:3], self.originAdj, self.discStep)
        outDistSq = (self.halfCubeLength + maxVdW_extraClash)**2
        distSq = np.sum((mol_xyzr[:, 0:3] - self.attach_xyz)**2, 1)
        grid_3d = self._carve_VdWextraClash(self.grid_3d, mol_xyzr, neighbor_list, ijk_atom, extraClash, distSq,
                                            outDistSq, self.shape)
        return grid_3d

    @staticmethod
    @nb.jit(nopython=True)
    def _carve_VdWextraClash(grid_3d, mol_xyzr, neighbor_list, ijk_atom, extraClash, distSq, outDistSq, grid_shape):
        """
        Loop through the atoms and assign -1 to all grid values that are within the atoms VdW radius
        plus an extraClash value
        """
        for m in range(mol_xyzr.shape[0]):
            if distSq[m] > outDistSq:
                continue
            for n in neighbor_list:
                if n[0] > mol_xyzr[m, 3] + extraClash:
                    break
                ijk_n = ijk_atom[m] + n[1]
                i, j, k = ijk_n
                if not np.any(ijk_n < 0) and not np.any(grid_shape - ijk_n <= 0):
                    grid_3d[i, j, k] = -1
        return grid_3d

    def dijkstra_init(self, maxR, maxR_source):
        """
        Initialize the Dijkstra search algorithm

        Parameters
        ----------
        maxR : float
               radius within to search for neighbors
        maxRsource : float
                     radius within to search for neighbors in the first round of the Dijkstra algorithm

        Returns
        -------
        grid_3d : ndarray
                  3-dimensional array of grid points with a shape of 2*adjL+1
        """
        ai, aj, ak = self._xyz2idx(self.attach_xyz, self.originAdj, self.discStep)
        self.grid_3d[ai, aj, ak] = 0
        priority_queue = [(0.0, (ai, aj, ak))]
        if _hasTypedList:
            # priority_queue = nb.typed.List()
            # [priority_queue.append(x) for x in [(0.0, (ai, aj, ak))]]
            edges_ess = nb.typed.List()
            [edges_ess.append(x) for x in self.essential_neighbors(self.sortedNeighborIdx(maxR), _dist_list)]
            edges_src = nb.typed.List()
            [edges_src.append(x) for x in self.sortedNeighborIdx(maxR_source)]
        else:
            edges_ess = self.essential_neighbors(self.sortedNeighborIdx(maxR), _dist_list)
            edges_src = self.sortedNeighborIdx(maxR_source)
        grid_3d = self.dijkstra(self.grid_3d, edges_ess, edges_src, priority_queue, (ai, aj, ak))
        return grid_3d

    def sortedNeighborIdx(self, maxR):
        """
        Build a neighbor list with a maximal extent given by maxR and sorted by increasing distance
        from the origin

        Parameters
        ----------
        maxR : float
               radius within to search for neighbors

        Returns
        -------
        idxs : list of 2-tuples of numpy.ndarray and float
               the tuples contain the ijk indices and the distance from the origin (0,0,0)
        """
        idxs = self._neighborIdx(maxR, self.discStep)
        idxs.sort(key=lambda x: x[0])
        return idxs

    @staticmethod
    @nb.jit(forceobj=True)
    def _neighborIdx(maxR, grid_spacing):
        """Build a list of neighboring indices to the origin (0,0,0)

        Parameters
        ----------
        maxR : float
            radius within which to search for neighbors
        grid_spacing : float
            spacing between grid points (in A)

        Returns
        -------
        idxs : list of 2-tuples of numpy.ndarray and float
            the tuples contain the ijk indices and the distance from the origin (0,0,0), e.g. [(1.0,array([-1,0,0]))]

        Notes
        -----
        The neighbor list is built from an origin set at (0,0,0). imaxR defines the number of indices
        along one axis (i,j or k) to build the index cube.
        For imaxR=3, first a cube of 7*7*7=343 indices is built (1).
        This is reduced in a second step to 123 indices that have a distance < imaxR from the origin (2).
        Finally, they are reduced 74 essential neighbors by the distances in _dist_list (3).
        (this last reduction speeds up the Dikstra algorithm due to the shorter neighbor list without loosing much
        accuracy)
        """
        idxs = []
        maxRSq = maxR**2
        imaxR = int(maxR / grid_spacing + 0.5)  # rounds to next integer
        # (1) build a cube with (2*imaxR+1)**3 indices
        for k in range(-imaxR, imaxR + 1, 1):
            for j in range(-imaxR, imaxR + 1, 1):
                for i in range(-imaxR, imaxR + 1, 1):
                    ijk = np.array([i, j, k], dtype=int)
                    # (2) only accept indices with 0 < dS < imaxRSq
                    dSq = np.sum(ijk**2) * grid_spacing**2
                    if dSq <= maxRSq:  # and dSq > 0:
                        d = np.sqrt(dSq)
                        idxs.append((d, ijk))
        return idxs

    def essential_neighbors(self, idxs, dist_list):
        """
        Reduce the neighbor list to those indices with a distance from the origin featured in dist_list

        Parameters
        ----------
        idxs : list of 2-tuples of numpy.ndarray and float
               the tuples contain the ijk indices and the distance from the origin (0,0,0)

        idxs_ess : list of 2-tuples of numpy.ndarray and float
                   list n essential neighbors from the origin (0,0,0)

        Notes
        -----
        If dist_list=sqrt([1, 2, 3, 5, 6]) then the the neighbor list is reduced to 74 neighbors
        which make up a spherical shape and are a good compromise between neighbor space and time efficiency
        in the Dijkstra algorithm. The smaller this list the faster the Dijkstra algorithm will run
        at the expense of the accuracy of the resulting accessible volume (i.e. its "sphericalness")
        """
        idxs_ess = []
        for e in idxs:
            minDiff = min(abs(e[0] - dist_list * self.discStep))
            if minDiff < (0.01 * self.discStep):
                idxs_ess.append(e)
        return idxs_ess

    @staticmethod
    @nb.jit(nopython=True)
    def dijkstra(grid_3d, edges_ess, edges_src, priority_queue, start_ijk):
        """
        Djikstra algorithm with a priority queue

        Parameters
        ----------
        grid_3d : numpy.ndarray
                  3-dimensional array of grid points with a shape given by n_xyz
        edges_ess : nb.typed.List of 2-tuples of numpy.ndarray and float
                    list n essential neighbors from the origin (0,0,0)
        edges_src : list of 2-tuples of numpy.ndarray and float
                    nb.typed.List of neighbors in the initialization round of the Dijkstra algorithm
        priority_queue : list of 2-tuples of numpy.ndarray and float
                         list with distance and index of origin

        Returns
        -------
        grid_3d : numpy.ndarray([nx,ny,nz])
                  3-dimensional array of grid points with a shape given by n_xyz
        """
        while priority_queue:
            r, idx = heapq.heappop(priority_queue)
            if r > grid_3d[idx]:
                continue
            if idx == start_ijk:
                edges = edges_src
            else:
                edges = edges_ess
            for e in edges:
                i = idx[0] + e[1][0]
                j = idx[1] + e[1][1]
                k = idx[2] + e[1][2]
                val = grid_3d[i, j, k]
                if val < 0:
                    continue
                t_r = r + e[0]
                if t_r < val:
                    grid_3d[i, j, k] = t_r
                    heapq.heappush(priority_queue, (t_r, (i, j, k)))
        return grid_3d

    @staticmethod
    @nb.jit(nopython=True)
    def setAboveTreshold(grid_3d, treshold, new_value):
        """
        Reassign grid values which are above a specified treshold

        Parameters
        ----------
        treshold : float
        new_value : float

        Returns
        -------
        grid_3d : numpy.ndarray([nx,ny,nz])
                  3-dimensional array of grid points with a shape given by n_xyz
        """
        for i in range(grid_3d.shape[0]):
            for j in range(grid_3d.shape[1]):
                for k in range(grid_3d.shape[2]):
                    if grid_3d[(i, j, k)] > treshold:
                        grid_3d[(i, j, k)] = new_value
        return grid_3d

    def excludeConcentricSpheres(self, mol_xyzr, dye_radii, maxRho):
        """
        Parameters
        ----------
        mol_xyzr : numpy.ndarray
                   array of x-,y-,z-coordinates and VdW radii with a shape [n_atoms, 4]
        dye_radii : ndarray([3,1])
                    array of dye radii in Angstrom with shape [3,1]
        maxRho : float
                 maximum grid value in the free volume

        Returns
        -------
        grid_3d : numpy.ndarray([nx,ny,nz])
                  3-dimensional array of grid points with a shape given by n_xyz
        """
        maxVdW_extraClash = np.max(mol_xyzr[:, 3]) + np.max(dye_radii)
        if _hasTypedList:
            neighbor_list = nb.typed.List()
            [neighbor_list.append(n) for n in self.sortedNeighborIdx(maxVdW_extraClash)]
            dye_radii_sorted = nb.typed.List()
            [dye_radii_sorted.append(x) for x in sorted(dye_radii)]
        else:
            neighbor_list = self.sortedNeighborIdx(maxVdW_extraClash)
            dye_radii_sorted = sorted(dye_radii)
        rhos = np.linspace(0, maxRho, len(dye_radii) + 1)
        ijk_atom = self._xyz2idx(mol_xyzr[:, 0:3], self.originAdj, self.discStep)
        outdistSq = (self.halfCubeLength + maxVdW_extraClash)**2
        distSq = np.sum((mol_xyzr[:, 0:3] - self.attach_xyz)**2, 1)
        grid_3d = self._assignRho(self.grid_3d, mol_xyzr, neighbor_list, ijk_atom, dye_radii_sorted, rhos, distSq,
                                  outdistSq, self.shape)
        return grid_3d

    @staticmethod
    @nb.jit(nopython=True)
    def _assignRho(grid_3d, mol_xyzr, neighbor_list, ijk_atom, dye_radii_sorted, rhos, distSq, outdistSq, grid_shape):
        """
        Loop through the atoms and radii and reassign the grid values with a number 0 < rho < 1
        """
        n = len(neighbor_list)
        n_radii = len(dye_radii_sorted)
        for m in range(mol_xyzr.shape[0]):
            if distSq[m] > outdistSq:
                continue
            else:
                s = 0
                for r in range(n_radii):
                    effR = dye_radii_sorted[r] + mol_xyzr[m, 3]
                    while neighbor_list[s][0] <= effR:
                        ijk_n = ijk_atom[m] + neighbor_list[s][1]
                        i, j, k = ijk_n
                        if not np.any(ijk_n < 0) and not np.any(grid_shape - ijk_n <= 0):
                            gridval = grid_3d[i, j, k]
                            grid_3d[i, j, k] = min(rhos[r], gridval)
                        s += 1
                        if s >= n:
                            break
        return grid_3d

    def addWeights(self, das_xyzrm):
        """
        Assign a weight to the grid values which are part of the contact volume

        Parameters
        ----------
        das_xyzrm : numpy.ndarray([5,n_atoms])
                    array of marked coordinates and padded vdW radii (n_atoms = number of atoms in mdtraj.Trajectory)

        Returns
        -------
        grid_3d : numpy.ndarray([nx,ny,nz])
                  3-dimensional array of grid points with a shape given by n_xyz
        """
        maxClash = np.max(das_xyzrm[:, 3]) + self.discStep
        if _hasTypedList:
            neighbor_list = nb.typed.List()
            [neighbor_list.append(n) for n in self.sortedNeighborIdx(maxClash)]
        else:
            neighbor_list = self.sortedNeighborIdx(maxClash)
        ijk_atom = self._xyz2idx(das_xyzrm[:, 0:3], self.originAdj, self.discStep)
        outDistSq = (self.halfCubeLength + maxClash)**2
        distSq = np.sum((das_xyzrm[:, 0:3] - self.attach_xyz)**2, 1)
        grid_3d = self._assignDensity(self.grid_3d, das_xyzrm, neighbor_list, ijk_atom, distSq, outDistSq, self.shape)
        return grid_3d

    @staticmethod
    @nb.jit(nopython=True)
    def _assignDensity(grid_3d, das_xyzrm, neighbor_list, ijk_atom, distSq, outDistSq, grid_shape):
        """
        Loop through the atoms and assign the density
        """
        n = len(neighbor_list)
        for m in range(das_xyzrm.shape[0]):
            if distSq[m] > outDistSq:
                continue
            else:
                s = 0
                while neighbor_list[s][0] <= das_xyzrm[m, 3]:
                    ijk_n = ijk_atom[m] + neighbor_list[s][1]
                    i, j, k = ijk_n
                    if not np.any(ijk_n < 0) and not np.any(grid_shape - ijk_n <= 0):
                        gridval = grid_3d[i, j, k]
                        if gridval > 0:
                            grid_3d[i, j, k] += das_xyzrm[m, 4]
                    s += 1
                    if s >= n:
                        break
        return grid_3d
