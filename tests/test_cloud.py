#!/usr/bin/env python3

import unittest
import os
import mdtraj as md
from fretraj import cloud

_TEST_DIR = os.path.abspath(os.path.dirname(__file__))

class LabelCheck(unittest.TestCase):

    def setUp(self):
        self.minimal_labels = {"Position":
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

    def test_simulation_type(self):
        cloud.check_labels(self.minimal_labels, verbose=False)
        self.assertEqual(self.minimal_labels['Position']['Cy3']['simulation_type'], 'AV3')
        global _importChecked
        _importChecked = True

    def test_Position_missing(self):
        labels = self.minimal_labels.pop('Position', None)
        self.assertRaises(ValueError, cloud.check_labels, self.minimal_labels, verbose=False)

    def test_linkerlength_missing(self):
        self.minimal_labels['Position']['Cy3'].pop('linker_length', None)
        self.assertRaises(KeyError, cloud.check_labels, self.minimal_labels, verbose=False)

    def test_attachID_wrongtype(self):
        self.minimal_labels['Position']['Cy3']['attach_id'] = '1'
        self.assertRaises(TypeError, cloud.check_labels, self.minimal_labels, verbose=False)

    def test_key_unrecognized(self):
        self.minimal_labels['Position']['Cy3']['unknown_key'] = None
        self.assertRaises(KeyError, cloud.check_labels, self.minimal_labels, verbose=False)

class VolumeCheck(unittest.TestCase):

    def setUp(self):
        self.minimal_labels = {"Position":
                            {"Cy3":
                                {"attach_id": 1,
                                 "linker_length": 20,
                                 "linker_width": 5,
                                 "dye_radius1": 8,
                                 "dye_radius2": 3,
                                 "dye_radius3": 3,
                                 "use_LabelLib": True,
                                 },
                             "Cy5":
                                {"attach_id": 288,
                                 "linker_length": 20,
                                 "linker_width": 5,
                                 "dye_radius1": 9.5,
                                 "dye_radius2": 3,
                                 "dye_radius3": 3,
                                 "use_LabelLib": True,
                                 },
                            },
                         "Distance": {"Cy3-Cy5":
                            {"R0": 54}
                            }
                         }
        cloud.check_labels(self.minimal_labels, verbose=False)
        self.struct = md.load_pdb(os.path.join(_TEST_DIR, 'testdata/DNA.pdb'))

    def test_volume_LabelLib(self):
         vol = cloud.Volume(self.struct, 'Cy3', self.minimal_labels)
         self.assertTrue(hasattr(vol.acv, 'mp'))

    def test_volume_pythononly(self):
         self.minimal_labels['Position']['Cy3']['use_LabelLib'] = False
         vol = cloud.Volume(self.struct, 'Cy3', self.minimal_labels)
         self.assertTrue(hasattr(vol.acv, 'mp'))

    def test_stateframe_mismatch(self):
         self.minimal_labels['Position']['Cy3']['state'] = self.minimal_labels['Position']['Cy3']['frame_mdtraj'] 
         vol = cloud.Volume(self.struct, 'Cy3', self.minimal_labels)
         self.assertIsNone(vol.acv)

    def test_stateframe_outofrange(self):
         self.minimal_labels['Position']['Cy3']['state'] = 2
         self.minimal_labels['Position']['Cy3']['frame_mdtraj'] = 1
         vol = cloud.Volume(self.struct, 'Cy3', self.minimal_labels)
         self.assertIsNone(vol.acv)

    def test_attachID_outofrange(self):
         self.minimal_labels['Position']['Cy3']['attach_id'] = 10**3
         vol = cloud.Volume(self.struct, 'Cy3', self.minimal_labels)
         self.assertIsNone(vol.acv)

    def test_attachID_buried(self):
         self.minimal_labels['Position']['Cy3']['attach_id'] = 158
         vol = cloud.Volume(self.struct, 'Cy3', self.minimal_labels)
         self.assertIsNone(vol.acv)

    def test_fromframes(self):
        vol = cloud.Volume.from_frames(self.struct, 'Cy3', self.minimal_labels, [0,0])
        self.assertIsInstance(vol[0], cloud.Volume)

    def test_fromattachID(self):
        vol = cloud.Volume.from_attachID(self.struct, 'Cy3', self.minimal_labels, [1,288])
        self.assertIsInstance(vol[0], cloud.Volume)


class FRETCheck(unittest.TestCase):

    def setUp(self):
        self.minimal_labels = {"Position":
                                {"Cy3":
                                    {"attach_id": 1,
                                     "linker_length": 20,
                                     "linker_width": 5,
                                     "dye_radius1": 8,
                                     "dye_radius2": 3,
                                     "dye_radius3": 3,
                                     "use_LabelLib": True,
                                     },
                                 "Cy5":
                                    {"attach_id": 288,
                                     "linker_length": 20,
                                     "linker_width": 5,
                                     "dye_radius1": 9.5,
                                     "dye_radius2": 3,
                                     "dye_radius3": 3,
                                     "use_LabelLib": True,
                                     },
                                },
                             "Distance": {"Cy3-Cy5":
                                {"R0": 54}
                                }
                             }

        cloud.check_labels(self.minimal_labels, verbose=False)
        self.struct = md.load_pdb(os.path.join(_TEST_DIR, 'testdata/DNA.pdb'))
        self.vol_Cy3 = cloud.Volume(self.struct, 'Cy3', self.minimal_labels)
        self.vol_Cy5 = cloud.Volume(self.struct, 'Cy5', self.minimal_labels)

    def test_FRET_LabelLib(self):
         fret = cloud.FRET_Trajectory(self.vol_Cy3, self.vol_Cy5, 'Cy3-Cy5', self.minimal_labels)
         self.assertTrue(hasattr(fret, 'R_mp'))

    def test_FRET_pythononly(self):
         self.vol_Cy3.use_LabelLib = False
         fret = cloud.FRET_Trajectory(self.vol_Cy3, self.vol_Cy5, 'Cy3-Cy5', self.minimal_labels)
         self.assertTrue(hasattr(fret, 'R_mp'))


if __name__ == '__main__':
    unittest.main()
