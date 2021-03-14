#!/usr/bin/env python3
import pytest
import os
from fretraj import burst
import jsonschema

_TEST_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

BURST_PARAMETERS = {
        "dyes": {
            "tauD": 0.75,
            "tauA": 1.5,
            "QD": 0.2,
            "QA": 0.3
        },
        "sampling": {
            "nbursts": 2000,
            "skipframesatstart": 0,
            "skipframesatend": 1000,
            "multiprocessing": True
        },
        "fret": {
            "R0": 5.4,
            "kappasquare": 0.666666,
            "no_gamma": True,
            "quenching_radius": 1.0
        },
        "species": {
            "name": ["all"],
            "unix_pattern_rkappa": ["*.dat"],
            "unix_pattern_don_coords": ["*donor*.xvg"],
            "unix_pattern_acc_coords": ["*acceptor*.xvg"],
            "probability": [1]
        },
        "bursts": {
            "lower_limit": 20,
            "upper_limit": 150,
            "lambda": -2.3,
            "QY_correction": False,
            "averaging": "all"
        }
    }


@pytest.fixture
def patch_open(mocker):
    mocker.patch('fretraj.burst.open')
    mocker.patch('json.load', return_value=BURST_PARAMETERS)


@pytest.fixture
def validate_parameters(patch_open):
    return burst.readParameters('anyname.json')


def test_donor_lifetime(validate_parameters):
    assert validate_parameters['dyes']['tauD'] == 0.75


def test_missing_key(validate_parameters, monkeypatch):
    monkeypatch.delitem(BURST_PARAMETERS['dyes'], 'tauD')
    with pytest.raises(jsonschema.exceptions.ValidationError):
        burst.readParameters('anyname.json')


def test_wrong_type(validate_parameters, monkeypatch):
    monkeypatch.setitem(BURST_PARAMETERS['dyes'], 'tauD', '0.75')
    with pytest.raises(jsonschema.exceptions.ValidationError):
        burst.readParameters('anyname.json')


def testTrajectory_rkappa():
    traj = burst.Trajectory.from_file(rkappa_filename=os.path.join(_TEST_DIR, 'data', 'test_rkappa1.dat'))
    assert traj.time[0] == 0


def testTrajectory_rkappa_dyecoords():
    traj = burst.Trajectory.from_file(rkappa_filename=os.path.join(_TEST_DIR, 'data', 'test_rkappa1.dat'),
                                      don_coords_filename=os.path.join(_TEST_DIR, 'data', 'test_donorcoords1.xvg'),
                                      acc_coords_filename=os.path.join(_TEST_DIR, 'data', 'test_acceptorcoords1.xvg'))
    assert traj.time[0] == 0


def testSpecies():
    species = burst.Species(name='all', probability=1,
                            filelist_rkappa=[os.path.join(_TEST_DIR, 'data', 'test_rkappa1.dat'), 
                                             os.path.join(_TEST_DIR, 'data', 'test_rkappa2.dat')])
    assert species.name == 'all'
    assert species.trajectories[0].time[0] == 0
    assert pytest.approx(species.trajectories[0].weight, 0.01) == 0.4


def testEnsemble():
    ensemble = burst.Ensemble(os.path.join(_TEST_DIR, 'data'),
                              {'species': {'name': ['all'], 'unix_pattern_rkappa': ['*.dat'],
                                           'probability': [1], 'n_trajectory_splits': None}})
    assert len(ensemble.species[0].trajectories) == 2


if __name__ == '__main__':
    pytest.main()
