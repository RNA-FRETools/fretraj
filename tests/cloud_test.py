#!/usr/bin/env python3
import sys

import unittest
import json
from fretraj import cloud
from fretraj import traj


class LabelCheck(unittest.TestCase):

    def setUp(self):
        with open('data/labels.json') as f:
            self.labels = json.load(f)
            self.pkey = '1'
            self.dkey = '1-59'

    def test_Position_missing(self):
        self.labels.pop('Position', None)
        labels = acv.check_labels(self.labels)
        self.assertEqual(labels['Position'], None)

    def test_SimulationType_missing(self):
        self.labels['Position'][self.pkey].pop('simulation_type', None)
        with self.assertRaises(KeyError):
            acv.check_labels(self.labels)

    def test_gridSpacing_default_value(self):
        self.labels['Position'][self.pkey].pop('grid_spacing', None)
        labels = acv.check_labels(self.labels)
        self.assertEqual(labels['Position'][self.pkey]['grid_spacing'], 0.5)

    def test_attachID_type(self):
        self.labels['Position'][self.pkey]['attach_id'] = "1"
        with self.assertRaises(TypeError):
            acv.check_labels(self.labels)

    def test_attachID_value(self):
        self.labels['Position'][self.pkey]['attach_id'] = 10
        labels = acv.check_labels(self.labels)
        self.assertEqual(labels['Position'][self.pkey]['attach_id'], 10)


if __name__ == '__main__':
    unittest.main()
