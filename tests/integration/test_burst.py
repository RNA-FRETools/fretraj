#!/usr/bin/env python3

import pytest
from fretraj import burst
import os

_TEST_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))


@pytest.fixture
def setup_parameters(mocker):
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
    mocker.patch('fretraj.burst.open')
    mocker.patch('json.load', return_value=BURST_PARAMETERS)
    return burst.readParameters('anyname.json')


def test_singleCore_withoutAnisotropy(setup_parameters):
    setup_parameters['sampling']['multiprocessing'] = False
    burst.Experiment(os.path.join(_TEST_DIR, 'data'), setup_parameters,
                     binwidth=0.025, verbose=False, compute_anisotropy=False)


def test_multiCore_withoutAnisotropy(setup_parameters):
    burst.Experiment(os.path.join(_TEST_DIR, 'data'), setup_parameters,
                     binwidth=0.025, verbose=False, compute_anisotropy=False)


def test_multiCore_withAnisotropy(setup_parameters):
    burst.Experiment(os.path.join(_TEST_DIR, 'data'), setup_parameters,
                     binwidth=0.025, verbose=False, compute_anisotropy=True)


if __name__ == "__main__":
    pytest.main()
