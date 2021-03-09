#!/usr/bin/env python3

import subprocess


def pymol_vis():
    subprocess.call('../skripts/pymol_vis.sh')


def vmd_vis():
    subprocess.call('../skripts/vmd_vis.sh')
