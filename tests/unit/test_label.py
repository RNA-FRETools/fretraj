#!/usr/bin/env python3
import pytest
import os
import argparse
import json
from fretraj import cloud


MINIMAL_LABELS = {"Position":
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


@pytest.fixture
def complete_labels():
    cloud.check_labels(MINIMAL_LABELS, verbose=False)
    return MINIMAL_LABELS


def test_simulation_type(complete_labels):
    assert complete_labels['Position']['Cy3']['simulation_type'] == 'AV3'


@pytest.mark.parametrize('label_dict, missing_item, error',
                         [(MINIMAL_LABELS, 'Position', ValueError),
                          (MINIMAL_LABELS['Position']['Cy3'], 'linker_length', KeyError)])
def test_missing_item(label_dict, missing_item, error, monkeypatch):
    monkeypatch.delitem(label_dict, missing_item)
    with pytest.raises(error):
        cloud.check_labels(MINIMAL_LABELS, verbose=False)


def test_attachID_wrongtype(monkeypatch):
    monkeypatch.setitem(MINIMAL_LABELS['Position']['Cy3'], 'linker_length', '1')
    with pytest.raises(TypeError):
        cloud.check_labels(MINIMAL_LABELS, verbose=False)


def test_unrecognized_key(monkeypatch):
    monkeypatch.setitem(MINIMAL_LABELS['Position']['Cy3'], 'unknown_key', 'some_value')
    with pytest.raises(KeyError):
        cloud.check_labels(MINIMAL_LABELS, verbose=False)


def test_AV1(complete_labels, monkeypatch):
    monkeypatch.setitem(MINIMAL_LABELS['Position']['Cy3'], 'simulation_type', 'AV1')
    cloud.check_labels(MINIMAL_LABELS, verbose=False)
    assert complete_labels['Position']['Cy3']['dye_radius2'] == 0


def test_parseCmd(mocker):
    mocker.patch('argparse.ArgumentParser.parse_args',
                 return_value=argparse.Namespace(input='test.pdb', parameters='param.json', output=None))
    out = cloud.parseCmd()
    assert out == ('test.pdb', 'param.json', None)


def test_labelingParams_valid(mocker):
    mocker.patch('fretraj.cloud.open')
    mocker.patch('json.load', return_value=MINIMAL_LABELS)
    labels_json = cloud.labeling_params('param.json')
    assert labels_json['Position']['Cy3']['attach_id'] == 1


def test_labelingParams_missing_key(mocker, monkeypatch):
    mocker.patch('fretraj.cloud.open')
    monkeypatch.delitem(MINIMAL_LABELS['Position']['Cy3'], 'linker_length')
    mocker.patch('json.load', return_value=MINIMAL_LABELS)
    labels_json = cloud.labeling_params('param.json')
    assert labels_json is None


def test_labelingParams_wrongtype(mocker, monkeypatch):
    mocker.patch('fretraj.cloud.open')
    monkeypatch.setitem(MINIMAL_LABELS['Position']['Cy3'], 'linker_length', '1')
    mocker.patch('json.load', return_value=MINIMAL_LABELS)
    labels_json = cloud.labeling_params('param.json')
    assert labels_json is None


if __name__ == '__main__':
    pytest.main()
