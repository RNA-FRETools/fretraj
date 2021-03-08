#!/usr/bin/env python3
import pytest
import numpy as np
from fretraj import export
from fretraj import cloud


CLOUD_XYZQT = np.array([[0.5, 1.0, 0.5, 1.0, 1],
                        [0.5, 1.0, 0.5, 8.0, 2]])
MP = CLOUD_XYZQT[0, 0:3]
ELEMENT = export.ELEMENT
ELEMENT_MP = export.ELEMENT_MP
NAME = export.NAME
NAME_MP = export.ELEMENT_MP


def test_xyz():
    s = export.xyz(CLOUD_XYZQT, MP)
    assert s.split('\n')[2] == f'{ELEMENT}\t0.500\t1.000\t0.500\t1.000'
    assert s.split('\n')[4] == f'{ELEMENT_MP}\t0.500\t1.000\t0.500\t-1.000'


def test_xyz_encodeElem():
    s = export.xyz(CLOUD_XYZQT, MP, encode_element=True)
    assert s.split('\n')[2] == 'F\t0.500\t1.000\t0.500\t1.000'
    assert s.split('\n')[3] == 'C\t0.500\t1.000\t0.500\t8.000'


def test_xyz_noweights():
    s = export.xyz(CLOUD_XYZQT, MP, write_weights=False)
    assert s.split('\n')[2] == f'{ELEMENT}\t0.500\t1.000\t0.500'
    assert s.split('\n')[4] == f'{ELEMENT_MP}\t0.500\t1.000\t0.500'


def test_xyz_noweights_encodeElem():
    s = export.xyz(CLOUD_XYZQT, MP, encode_element=True, write_weights=False)
    assert s.split('\n')[2] == 'F\t0.500\t1.000\t0.500'
    assert s.split('\n')[3] == 'C\t0.500\t1.000\t0.500'


def test_pdb():
    s = export.pdb(CLOUD_XYZQT, MP)
    assert s.split('\n')[0] == f'ATOM      1    {NAME}  FV     1       0.500   1.000   0.500 99.00  1.00           {ELEMENT}  '
    assert s.split('\n')[2] == f'ATOM      3   {NAME_MP}  MP     0       0.500   1.000   0.500 99.00 -1.00          {ELEMENT_MP}  '


def test_openDX():
    grid_3d = np.zeros((3, 3, 3))
    grid_3d[(0, 0, 0)] = 1
    grid_3d[(0, 0, 1)] = 10
    s = export.open_dx(grid_3d, (0, 0, 0), (1, 1, 1))
    assert s.split('\n')[1] == 'origin 0.000 0.000 0.000'
    assert s.split('\n')[7] == '1.0 10.0 0.0'


if __name__ == '__main__':
    pytest.main()
