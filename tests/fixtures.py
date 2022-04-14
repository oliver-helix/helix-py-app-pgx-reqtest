import pytest
import json

from configparser import ConfigParser, ExtendedInterpolation

from helix_lib_core import fileops
from helix_lib_bioinfo.wrappers.config import load_tool_config

from helix_app_pgx.gene_config import construct_gene_configs
from helix_app_pgx.structures import Batch


@pytest.fixture
def config():
    config = ConfigParser(interpolation=ExtendedInterpolation())
    config_files = ['helix_app_pgx/config/reference.ini',
                    'helix_app_pgx/config/workflow.ini',
                    'helix_app_pgx/config/tools.ini',
                    'helix_app_pgx/config/platform/test.ini']
    config.read(config_files)
    load_tool_config(config_files)
    # make sure unit test scratch/output does not mix with production output
    config['output']['directory'] = '/tests/scratch'
    config['scratch']['directory'] = '/tests/scratch'
    fileops.makedirs(config['scratch']['directory'])
    fileops.makedirs(config['output']['samples_directory'])
    # we won't use actual s3 path below for testing
    config['output']['s3_dir'] = 's3_dir'
    return config


@pytest.fixture
def gene_configs(config):
    gene_configs = construct_gene_configs('CYP2C19,CYP2D6,VKORC1',
                                          config['reference']['prefix'],
                                          config['output']['directory'],
                                          config['scratch']['directory'])
    for gene_config in gene_configs.values():
        fileops.makedirs(gene_config.scratch_dir)
    return gene_configs


@pytest.fixture
def cyp2d6_gene_config(gene_configs):
    return gene_configs['CYP2D6']


@pytest.fixture
def cyp2c19_gene_config(gene_configs):
    return gene_configs['CYP2C19']


@pytest.fixture
def variant_gene_config(gene_configs):
    return gene_configs['VKORC1']


@pytest.fixture
def batch1(gene_configs):
    with open('tests/data/test_ev2val_sp1_batch.json') as json_file:
        input_json = json.load(json_file)
    batch1 = Batch(input_json, gene_configs)
    batch1.gene_batches['CYP2D6'].load_counts_df('tests/data/test_ev2val_sp1_counts_data.txt.gz')
    return batch1


@pytest.fixture
def batch2(gene_configs):
    """Batch with single sample for testing failure routine."""

    with open('tests/data/batch2.json') as json_file:
        input_json = json.load(json_file)

    return Batch(input_json, gene_configs)


@pytest.fixture
def batch3(gene_configs):
    """Batch with more than one assay version."""

    with open('tests/data/batch3.json') as json_file:
        input_json = json.load(json_file)

    return Batch(input_json, gene_configs)


@pytest.fixture
def cyp2c19_params(config, cyp2c19_gene_config):
    """Dictionary used during allele calling step."""
    return {'min_dp_per_copy': config['parameters'].getint('min_dp_per_copy'),
              'alt_threshold': config['parameters'].getfloat('alt_threshold'),
              'nomatch_pvalue': config['parameters'].getfloat('nomatch_pvalue'),
              'max_low_pvalue': config['parameters'].getint('max_low_pvalue'),
              'gene_copies': range(1, cyp2c19_gene_config.max_gene_copies + 1),
              'paralog_copies': range(0, cyp2c19_gene_config.max_gene_copies + 1),
              'assay_version': 'v3',   # CYP2C19 unit test cases come from v3
              'gene_config': cyp2c19_gene_config}


@pytest.fixture
def cyp2d6_params(config, cyp2d6_gene_config):
    """Dictionary used during allele calling step."""
    return {'min_dp_per_copy': config['parameters'].getint('min_dp_per_copy'),
              'alt_threshold': config['parameters'].getfloat('alt_threshold'),
              'nomatch_pvalue': config['parameters'].getfloat('nomatch_pvalue'),
              'max_low_pvalue': config['parameters'].getint('max_low_pvalue'),
              'gene_copies': range(1, cyp2d6_gene_config.max_gene_copies + 1),
              'paralog_copies': range(0, cyp2d6_gene_config.max_gene_copies + 1),
              'assay_version': 'v2',   # CYP2D6 unit test cases come from v2
              'gene_config': cyp2d6_gene_config}
