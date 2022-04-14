import json

from pathlib import Path

from fixtures import config, gene_configs, batch1

def test_make_directory_structure(batch1, config):
    test_directory_json = batch1.make_directory_structure(config)
    with open('tests/data/correct_ev2val_sp1_directory.json') as fin:
        correct_output = json.load(fin)
    with open(test_directory_json) as fin:
        test_output = json.load(fin)

    # remove version because this will change with each PR
    del correct_output['Batch']['pgxVersion']
    del test_output['Batch']['pgxVersion']

    assert sorted(test_output.items()) == sorted(correct_output.items())
    Path(test_directory_json).unlink()
