#!/usr/bin/env python3

import subprocess


def pymol_vis():
    subprocess.call('pymol_vis.sh')


def vmd_vis():
    subprocess.call('vmd_vis.sh')
