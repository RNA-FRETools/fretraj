#!/usr/bin/env python3

import numpy as np
import mdtraj as md
import copy

try:
    from fretraj import export
except ImportError: # for PyMOL plugin
    from . import export

class Plumed:

    def __init__(self, structure, volumes, selection='name P', cutoff=1, pseudo_element=md.element.bromine):
        self.struct = copy.deepcopy(structure)
        chain = self.struct.top.add_chain()
        pseudo_serials = []
        self.dist_restraints = []
        self.xyz = self.struct.xyz * 10
        for volume in volumes:
            res = self.struct.top.add_residue(pseudo_element.symbol, chain, max([res.resSeq for res in self.struct.top.residues])+1)
            atom = self.struct.top.add_atom(pseudo_element.symbol, pseudo_element, res, self.struct.top.n_atoms+1)
            pseudo_serials.append(atom.serial)
            self.dist_restraints.append(self.distance_restraint(atom.serial, volume, selection, cutoff))
            self.xyz = np.hstack((self.xyz, np.array([[volume.acv.mp]])))
        self.n_atoms = self.struct.top.n_atoms


    def distance_restraint(self, pseudo_serial, volume, selection='name P', cutoff=15):
        """
        cutoff : float
                 radius within which to search for atoms in the selection (in Angstrom)
        """
        attach_id_mdtraj = volume.attach_id-1
        neighbor_idx = md.compute_neighbors(self.struct, cutoff / 10, [attach_id_mdtraj])[0]
        atomsele_idx = volume.structure.top.select(selection)
        atoms_idx = np.array([idx for idx in neighbor_idx if idx in atomsele_idx])
        atoms_ser = atoms_idx+1
        atoms_xyz = volume.structure.xyz[0,atoms_idx,:] * 10
        distances = np.sqrt(np.sum((volume.acv.mp - atoms_xyz)**2, axis=1)) / 10
        return {'pseudo_atom': {'xyz': volume.acv.mp / 10}, 
                'attach_atoms': {'serial': atoms_ser, 'xyz': volume.structure.xyz[0,atoms_idx,:]},
                'distances': distances}

    def write_plumed(self, filename, R_mp, kappa_dist=1000, Rmp_kappa=500):
        """
        """
        plumed_str = ''
        k = 1
        plumed_str += 'mp1: GROUP ATOMS=xxx\n'
        plumed_str += 'mp2: GROUP ATOMS=yyy\n\n'
        for dr in self.dist_restraints:
            for i,attach_atom in enumerate(dr['attach_atoms']['serial']):
                plumed_str += 'd{:d}{:d}: DISTANCE ATOMS=mp{:d},{}\n'.format(k,i,k,attach_atom)
            for i,dist in enumerate(dr['distances']):
                plumed_str += 'r{:d}{:d}: RESTRAINT ARG=d{:d}{:d} AT={:0.1f} KAPPA={}\n'.format(k,i,k,i,dist,kappa_dist)
            plumed_str += '\n'
            k += 1
        plumed_str += 'WHOLEMOLECULES ENTITY0=1-{:d}\n'.format(self.n_atoms)
        plumed_str += 'Rmp: DISTANCE ATOMS=mp1,mp2 NOPBC\n'
        plumed_str += 'RESTRAINT ARG=Rmp AT={:0.1f} KAPPA={}\n'.format(R_mp, Rmp_kappa)
        with open(filename, 'w') as f:
            f.write(plumed_str)

    def write_pseudo(self, pdb_filename, itp_filename, pseudo_element='MP'):
        with open(pdb_filename, 'w') as f:
            pdbstr = export._pdb_format.format('ATOM', 1, pseudo_element, ' ', pseudo_element, ' ', 1, ' ', 0.000, 0.000, 0.000, 1, 0, ' ', ' ')
            f.write(pdbstr)
        with open(itp_filename, 'w') as f:
            itpstr = '''[ moleculetype ]
; molname       nrexcl
MP             1       ; pseudoatom representing the dye mean position

[ atoms ]
; id    at type         res nr  residue name     at name  cg nr  charge
1       MP              1       MP             MP       1      0.00000
'''
            f.write(itpstr)
        print('(1) Copy the \"{}\" file into your force field directory\n'.format(itp_filename.split('/')[-1]))
        print('(2) Add the following lines to your \"topology.top\" file:\n; Include topology for dye mean position\n#include "<path/to/forcefield>/{}"\n\n#ifdef POSRES\n[ position_restraints ]\n1     1  1000  1000  1000\n#endif\n'.format(itp_filename.split('/')[-1]))
        print('(3) Add the following line to your \"atomtypes.atp\" file in the force field directory:\n{}        0.00000 ;\n'.format(pseudo_element))
        print('(4) Add the following lines to your \"ffnonbonded.itp\" file in the force field directory:\n; Dummy mass for dye mean position\n{}          35      79.90    0.0000  A   4.64693e-01  2.45414e-01'.format(pseudo_element))


    # def write_pseudoGRO(self, filename):
    #     pdbfile = md.formats.GroTrajectoryFile(filename, mode='w')
    #     pdbfile.write(self.xyz / 10, self.struct.top)


    def write_vmd(self, filename='vis.vmd'):
        """
        Generate a visualization file for VMD containing Tcl API commands
        """
        vmd_vis = '''color change rgb blue    0.20  0.33  0.72
color change rgb red     0.76  0.33  0.29
color change rgb gray    0.71  0.74  0.77
color change rgb orange  0.89  0.50  0.32
color change rgb green   0.42  0.70  0.51
color change rgb violet  0.55  0.39  0.70
axes location off

mol modselect 0 0 all not solvent
mol modstyle 0 0 Licorice 0.300000 12.000000 12.000000
mol modcolor 0 0 ColorID 2
mol modmaterial 0 0 AOChalky

mol addrep 0
mol modselect 1 0 {resname MP}
mol modcolor 1 0 ColorID 0
mol modstyle 1 0 VDW 1.500000 12.000000

set sel [atomselect 0 "resname MP"]
'''
        for i in range(2):
            vmd_vis += 'set pseudoid{:d} [lindex [$sel get index] {:d}]\n'.format(i,i)
            for serial in self.dist_restraints[i]['attach_atoms']['serial']:
                vmd_vis += 'label add Bonds 0/$pseudoid{:d} 0/{:d}\n'.format(i,serial-1)
        vmd_vis += 'label add Bonds 0/$pseudoid0 0/$pseudoid1\n'
        vmd_vis += 'label textsize 0.5\n'
        with open(filename, 'w') as f:
            f.write(vmd_vis)


    def write_pymol(self, filename='vis.py'):
        """
        Generate a visualization file for PyMOL containing Python API commands
        """
        pymol_vis = '''cmd.hide(\'lines\')
cmd.show(\'sticks\', \'all and not hydrogen\')
cmd.remove(\'solvent or (inorganic and not resn MP)\')
cmd.set(\'cartoon_ring_mode\',3)
cmd.set(\'cartoon_ring_finder\',0)
cmd.color(\'gray80\')
cmd.show(\'spheres\', \'resn MP\')
cmd.set(\'sphere_scale\', 1.5, \'resn MP\')
cmd.color('skyblue', \'resn MP\')
cmd.set(\'dash_radius\', 0.3)
cmd.set(\'dash_gap\', 1)
'''
        for i in range(2):
            pymol_vis +='cmd.pseudoatom(\'p{:d}\', pos=[{:0.2f},{:0.2f},{:0.2f}])\n'.format(i, *self.dist_restraints[i]['pseudo_atom']['xyz'] * 10)
            for k,serial in enumerate(self.dist_restraints[i]['attach_atoms']['serial']):    
                pymol_vis += 'cmd.distance(\'r{:d}{:d}\', \'p{:d}\', \'id {:d}\')\n'.format(i+1,k+1,i,serial)
        pymol_vis += 'cmd.color(\'gray30\', \'r*\')\n'
        pymol_vis += 'cmd.hide(\'labels\', \'r*\')\n'
        with open(filename, 'w') as f:
            f.write(pymol_vis)