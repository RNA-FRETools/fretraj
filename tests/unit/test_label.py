#!/usr/bin/env python3
import pytest
import os
import mdtraj as md
from fretraj import cloud

_TEST_DIR = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture(scope='function')
def create_labels():
    minimal_labels = {"Position":
                            {"Cy3":
                                {"attach_id": 1,
                                 "linker_length": 20,
                                 "linker_width": 5,
                                 "dye_radius1": 8,
                                 "dye_radius2": 3,
                                 "dye_radius3": 3,
                                 },
                             "Cy5":
                                {"attach_id": 288,
                                 "linker_length": 20,
                                 "linker_width": 5,
                                 "dye_radius1": 9.5,
                                 "dye_radius2": 3,
                                 "dye_radius3": 3,
                                 },
                            },
                         "Distance": {"Cy3-Cy5":
                            {"R0": 54}
                            }
                        }
    return minimal_labels

def test_simulation_type(create_labels):
    cloud.check_labels(create_labels, verbose=False)
    assert create_labels['Position']['Cy3']['simulation_type'] == 'AV3'

def test_Position_missing(create_labels):
    create_labels.pop('Position', None)
    with pytest.raises(ValueError):
        cloud.check_labels(create_labels, verbose=False)

def test_linkerlength_missing(create_labels):
    create_labels['Position']['Cy3'].pop('linker_length', None)
    with pytest.raises(KeyError):
        cloud.check_labels(create_labels, verbose=False)

def test_attachID_wrongtype(create_labels):
    create_labels['Position']['Cy3']['attach_id'] = '1'
    with pytest.raises(TypeError):
        cloud.check_labels(create_labels, verbose=False)

def test_key_unrecognized(create_labels):
    create_labels['Position']['Cy3']['unknown_key'] = None
    with pytest.raises(KeyError):
        cloud.check_labels(create_labels, verbose=False)


if __name__ == '__main__':
    pytest.main()
