#!/usr/bin/env python3

import subprocess
import os

package_directory = os.path.dirname(os.path.abspath(__file__))


def pymol_vis():
    subprocess.call(os.path.join(package_directory, 'skripts', 'pymol_vis.sh'))


def vmd_vis():
    subprocess.call(os.path.join(package_directory, 'skripts', 'vmd_vis.sh'))
