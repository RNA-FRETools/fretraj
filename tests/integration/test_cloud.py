#!/usr/bin/env python3

import pytest
import os
import mdtraj as md
from fretraj import cloud

_TEST_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))


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
                        "use_LabelLib": False,
                        },
                       "Cy5":
                       {"attach_id": 288,
                        "linker_length": 20,
                        "linker_width": 5,
                        "dye_radius1": 9.5,
                        "dye_radius2": 3,
                        "dye_radius3": 3,
                        "use_LabelLib": False,
                        },
                       },
                      "Distance": {"Cy3-Cy5":
                                   {"R0": 54}
                                   }
                      }
    cloud.check_labels(minimal_labels, verbose=False)
    return minimal_labels


@pytest.fixture(scope='module')
def load_structure():
    return md.load_pdb(os.path.join(_TEST_DIR, 'data/DNA.pdb'))


@pytest.fixture(scope='function')
def create_volumes(create_labels, load_structure):
    vol_Cy3 = cloud.Volume(load_structure, 'Cy3', create_labels)
    vol_Cy5 = cloud.Volume(load_structure, 'Cy5', create_labels)
    return vol_Cy3, vol_Cy5


class TestVolume:

    @pytest.mark.skipif(not cloud._LabelLib_found, reason="requires LabelLib")
    def test_volume_LabelLib(self, create_labels, load_structure):
        create_labels['Position']['Cy3']['use_LabelLib'] = True
        vol = cloud.Volume(load_structure, 'Cy3', create_labels)
        assert hasattr(vol.acv, 'mp')

    def test_volume_pythononly(self, create_labels, load_structure):
        vol = cloud.Volume(load_structure, 'Cy3', create_labels)
        assert hasattr(vol.acv, 'mp')

    def test_stateframe_mismatch(self, create_labels, load_structure):
        create_labels['Position']['Cy3']['state'] = create_labels['Position']['Cy3']['frame_mdtraj']
        vol = cloud.Volume(load_structure, 'Cy3', create_labels)
        assert vol.acv is None

    def test_stateframe_outofrange(self, create_labels, load_structure):
        create_labels['Position']['Cy3']['state'] = 2
        create_labels['Position']['Cy3']['frame_mdtraj'] = 1
        vol = cloud.Volume(load_structure, 'Cy3', create_labels)
        assert vol.acv is None

    def test_attachID_outofrange(self, create_labels, load_structure):
        create_labels['Position']['Cy3']['attach_id'] = 10**3
        vol = cloud.Volume(load_structure, 'Cy3', create_labels)
        assert vol.acv is None

    def test_attachID_buried(self, create_labels, load_structure):
        create_labels['Position']['Cy3']['attach_id'] = 158
        vol = cloud.Volume(load_structure, 'Cy3', create_labels)
        assert vol.acv is None

    def test_fromframes(self, create_labels, load_structure):
        vol = cloud.Volume.from_frames(load_structure, 'Cy3', create_labels, [0, 0])
        assert isinstance(vol[0], cloud.Volume)

    def test_fromattachID(self, create_labels, load_structure):
        vol = cloud.Volume.from_attachID(load_structure, 'Cy3', create_labels, [1, 288])
        assert isinstance(vol[0], cloud.Volume)


class TestFRET:

    def test_FRET_LabelLib(self, create_volumes, create_labels):
        fret = cloud.FRET(create_volumes[0], create_volumes[1], 'Cy3-Cy5', create_labels)
        assert hasattr(fret, 'R_mp')

    def test_FRET_pythononly(self, create_volumes, create_labels):
        create_volumes[0].use_LabelLib = False
        fret = cloud.FRET(create_volumes[0], create_volumes[1], 'Cy3-Cy5', create_labels)
        assert hasattr(fret, 'R_mp')


if __name__ == '__main__':
    pytest.main()
