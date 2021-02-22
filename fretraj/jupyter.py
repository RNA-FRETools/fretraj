#!/usr/bin/env python3

# start a PyMOL server session from a terminal:
#   pymol -R

# On Windows you may create a shortcut that executes the following command:
#   C:\path\to\PyMOLWin.exe -R

import os
import re


def connect2pymol():
    import xmlrpc.client as xmlrpclib
    cmd = xmlrpclib.ServerProxy('http://localhost:9123')
    curr_wd = os.getcwd()
    try:
        cmd.cd(curr_wd)
    except:
        cmd.cd(re.sub(r'/mnt/([a-z])', r'\1:', curr_wd))
    return cmd
