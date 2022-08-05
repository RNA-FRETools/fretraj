#!/usr/bin/env python3

import subprocess
import os
import argparse
import platform

package_directory = os.path.dirname(os.path.abspath(__file__))

from fretraj import metadata
from fretraj import _LabelLib_found
from fretraj import _nglview_found
from fretraj import _burst_module_available
from fretraj import _jupyter_module_available
from fretraj import cloud
import mdtraj as md


def parseCmd():
    """Parse the command line to get the input PDB file and the dye parameters

    Returns
    -------
    tuple of str
        tuple containing the filenames of the input PDB and the parameter file, the labeling site identifier
        and the filename of the output ACV, i.e. (pdb_input, param_file, site, pdb_output)
    """
    parser = argparse.ArgumentParser(
        description="compute accessible-contact clouds for an MD trajectory or a given PDB structure"
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + str(metadata["Version"]))
    parser.add_argument(
        "--path", action="version", version=f"package directory: {package_directory}", help="Show package directory"
    )
    parser.add_argument("-i", "--input", help="Input PDB structure (.pdb)", required=True)
    parser.add_argument("-p", "--parameters", help="Parameter file for labels (.json)", required=True)
    parser.add_argument("-s", "--site", help="reference identifier for the labeling position (string)", required=True)
    parser.add_argument(
        "-o", "--output", help="Output file of accessible contact volume (.pdb, .xyz, .dx)", required=False
    )
    parser.add_argument(
        "--show-config",
        action="version",
        version=f"labellib: {_LabelLib_found} | nglview: {_nglview_found} | burst submodule: {_burst_module_available} | jupyter submodule: {_jupyter_module_available}",
        help="Show installed libraries and submodules",
    )
    args = parser.parse_args()
    in_filePDB = args.input
    param_fileJSON = args.parameters
    out_fileACV = args.output
    site = args.site
    return (in_filePDB, param_fileJSON, site, out_fileACV)


def main():
    """Calcualting an ACV from the command line

    Notes
    -----
    Run with:
    >>> fretraj -i <input pdb> -p <label parameter file (.json)> -s <label position identifier>
    """
    in_filePDB, param_fileJSON, site, out_fileACV = parseCmd()
    labels = cloud.labeling_params(param_fileJSON)
    struct = md.load_pdb(in_filePDB)
    print("Calculating accessible-contact volume...")
    av = cloud.Volume(struct, site, labels)
    if not out_fileACV:
        out_fileACV = "ACV.pdb"
    _, out_fileextACV = os.path.splitext(out_fileACV)
    out_path = os.path.dirname(os.path.abspath(out_fileACV))
    av.save_acv(out_fileACV, format=out_fileextACV[1:])
    print(f"{out_fileACV} saved to {out_path}")


def pymol_vis():
    """Visualize a PDB file or trajectory in PyMOL

    Notes
    -----
    Run with:
    >>> pymol_vis -c <structure file (.gro/.pdb)> -x <xtc trajectory> (optional) -v <pymol visualization file> (optional) -s <start,stop,stride (default 1,-1,1)>
    """
    if platform.system() != "Windows":
        subprocess.call(os.path.join(package_directory, "skripts", "pymol_vis.sh"))
    else:
        print("pymol_vis is only available for Unix operating systems")


def vmd_vis():
    """Visualize a PDB file or trajectory in VMD

    Notes
    -----
    Run with:
    >>> vmd_vis -c <structure file (.gro/.pdb)> -x <xtc trajectory> (optional) -v <vmd visualization file> (optional) -s <start,stop,stride (default 1,-1,1)>
    """
    if platform.system() != "Windows":
        subprocess.call(os.path.join(package_directory, "skripts", "vmd_vis.sh"))
    else:
        print("vmd_vis is only available for Unix operating systems")
