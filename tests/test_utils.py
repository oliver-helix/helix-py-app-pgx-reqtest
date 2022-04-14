import json
import shutil
import os
import filecmp
import subprocess

from collections import OrderedDict
from pathlib import Path
import pandas as pd
import pytest

from helix_app_pgx.utils import (create_empty_output_star, perform_failure_routine, get_sample_assay_versions,
                                 get_sample_flowcells, parse_bam_header, establish_assay_version, create_empty_output_variants,
                                 get_sample_gene_qc_metrics_variants, get_sample_gene_qc_metrics_star, compress_sample_files,
                                 print_jsons)
from helix_app_pgx.structures import Sample
from helix_app_pgx.call_alleles import calculate_likelihoods, create_sample_output_star, create_sample_output_variants
from helix_app_pgx.read_data import read_allele_depths, read_variants

from fixtures import (config, gene_configs, batch2, batch3, cyp2d6_gene_config, cyp2d6_params, cyp2c19_gene_config,
                      cyp2c19_params, variant_gene_config)


def test_get_sample_assay_versions(batch3):
    assert get_sample_assay_versions(batch3.get_gene_batches()[0]) == ['v1','v2']


def test_parse_bam_header():
    header_line = ('@RG\tID:3\tCN:HELIX_LAB_1\t'
                   'LB:PC0000032-DNA-B7-de88ce20-a324-43fa-8571-902f25151997\t'
                   'PG:2.8.1\tPL:ILLUMINA\tPM:HiSeq\tPU:H3JGMCFXX-4\t'
                   'SM:NA12878-PC0000032-B07')
    assert parse_bam_header(header_line) == {'CN': 'HELIX_LAB_1',
                                             'ID': '3',
                                             'LB': 'PC0000032-DNA-B7-de88ce20-a324-43fa-8571-902f25151997',
                                             'PG': '2.8.1',
                                             'PL': 'ILLUMINA',
                                             'PM': 'HiSeq',
                                             'PU': 'H3JGMCFXX-4',
                                             'SM': 'NA12878-PC0000032-B07'}

def test_get_sample_flowcells():
    assert get_sample_flowcells('tests/data/NA12878_slice.bam') == ['H3JGMCFXX']


def test_establish_assay_version(config, batch3):
    # batch3 originally has two samples, NA12878 and CL-13.
    # CL-13 was processed on a flowcell different to the batch flowcell.

    # Since these samples are duplicated across gene batches, we choose the first gene batch.
    assert sorted(batch3.get_gene_batches()[0].samples_for_analysis) == ['CL-13_133645816', 'NA12878_132964731']
    establish_assay_version(config, batch3)
    # After this function, CL-13 is identified as having an invalid assay version. It will be excluded
    # from further analysis and assay version will be inferred as v1 (NA12878's assay version)
    gene_batch_after = batch3.get_gene_batches()[0]
    assert gene_batch_after.samples_for_analysis == ['NA12878_132964731']
    assert gene_batch_after.samples['NA12878_132964731'].results['filter'] is None
    assert gene_batch_after.samples['CL-13_133645816'].results['filter'] == 'ASSAYINVALID'
    assert config['reference']['assay_version'] == 'v1'


def test_create_empty_output_star_lowcov():
    # The arguments in the next line are irrelevant to the test
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt', 'qctable_by_region.txt',
                    'v2')
    sample.results['filter'] = 'LOWCOV'
    gene_symbol = 'CYP2D6'
    create_empty_output_star(sample, gene_symbol)

    for k,v in sample.output.items():
        if k == 'filter':
            assert v[0] == 'LOWCOV'
        elif k != 'gene':
            assert v == None


def test_create_empty_output_star_lowq():
    # The arguments in the next line are irrelevant to the test
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt', 'qctable_by_region.txt',
                    'v2')
    sample.results['filter'] = 'LOWQ'
    gene_symbol = 'CYP2D6'
    create_empty_output_star(sample, gene_symbol)

    for k,v in sample.output.items():
        if k == 'filter':
            assert v[0] == 'LOWQ'
        elif k !='gene':
            assert v == None


def test_create_empty_output_star_default():
    # The arguments in the next line are irrelevant to the test
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v2')
    gene_symbol = 'CYP2D6'
    create_empty_output_star(sample, gene_symbol)

    for k, v in sample.output.items():
        if k != 'gene':
            assert v == None


def test_create_empty_output_variants():
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v2')
    gene_symbol = 'VKORC1'
    create_empty_output_variants(sample, gene_symbol)

    expected_result = OrderedDict([('gene', gene_symbol),
                                   ('genotypes', None),
                                   ('filter', None),
                                  ])
    assert sample.output == expected_result


def test_get_sample_gene_qc_metrics_variants():
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v2')
    gene_symbol = 'VKORC1'
    sample.output['filter'] = ['PASS']

    result = get_sample_gene_qc_metrics_variants(sample, gene_symbol)
    expected_result = OrderedDict([('gene', gene_symbol),
                                   ('noCandidates', 'False'),
                                   ('filter', ['PASS']),
                                   ('geneQc', 'PASS')
                                  ])

    assert result == expected_result

@pytest.mark.parametrize("filter,expected", [
    (['PASS'], 'PASS'),
    (['NOVEL'], 'PASS'),
    (['NOMATCH'], 'PASS'),
    (['NOMATCH', 'NOVEL'], 'PASS'),
    (['LOWQ'], 'FAIL'),
    (['LOWCOV'], 'FAIL'),
])
def test_get_sample_gene_qc_metrics_star(filter, expected):
    # Test that correct geneQc value is applied (when not in failure mode).
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v3')
    gene_symbol = 'CYP2D6'
    sample.output['filter'] = filter

    result = get_sample_gene_qc_metrics_star(sample, gene_symbol)

    assert result['filter'] == filter
    assert result['geneQc'] == expected


@pytest.mark.parametrize("filter", [
    ['PASS'],
    ['NOVEL'],
    ['NOMATCH'],
    ['NOMATCH', 'NOVEL'],
    ['LOWQ'],
    ['LOWCOV'],
])
def test_get_sample_gene_qc_metrics_star_fail_mode(filter):
    # Test that correct filter and geneQc value is applied (when in failure mode).
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v3')
    gene_symbol = 'CYP2D6'
    sample.output['filter'] = filter

    result = get_sample_gene_qc_metrics_star(sample, gene_symbol, failure_status=True)

    assert result['filter'] is None
    assert result['geneQc'] == 'FAIL'


def test_print_jsons(batch2, config, cyp2d6_params, cyp2c19_params, variant_gene_config):
    # Uses single-sample batch with input data cobbled together from various other tests.
    sample_id = 'CL-69-Plate_1-B01_133052682'

    # CYP2D6 and CYP2C19 calls use variants and copy number input data. VKORC1 only requires variant data.
    input_data = {'CYP2D6': {'variants': 'tests/data/novel_nomatch_vgraph.txt',
                             'copy_number': 'tests/data/novel_nomatch_exon_copy_emit.txt'},
                  'CYP2C19': {'variants': 'tests/data/CYP2C19_partial_deletion_vgraph.txt',
                              'copy_number': 'tests/data/CYP2C19_partial_deletion_exon_copy_emit.txt'},
                  'VKORC1': {'variants': 'tests/data/vkorc1_variants_nocall.txt',}
                  }

    for gene_symbol, params in zip(['CYP2D6', 'CYP2C19'],
                                   [cyp2d6_params, cyp2c19_params]):
        gene_config = params['gene_config']
        sample = batch2.get_gene_batch(gene_symbol).samples[sample_id]
        muts_all_df, novel_muts = read_allele_depths(input_data[gene_symbol]['variants'],
                                                    gene_config.mutations_lookup)
        sample.data['muts_all_df'] = muts_all_df
        sample.data['novel_muts'] = novel_muts
        sample.data['exon_emit_mat'] = pd.read_csv(input_data[gene_symbol]['copy_number'],
                                                    sep='\t')
        sample.results = calculate_likelihoods(sample, params, gene_config.scratch_dir)
        create_sample_output_star(sample, gene_config)

    # VKORC1 gene
    sample = batch2.get_gene_batch('VKORC1').samples[sample_id]
    sample.data['r2v_calls'] = read_variants(input_data['VKORC1']['variants'])
    create_sample_output_variants(sample, variant_gene_config)

    print_jsons(batch2, config)

    # Test failedGenes logic.
    qc_json = '{}/{}_{}'.format(config['output']['samples_directory'], sample_id, config['suffixes']['qc_json'])
    expected_output = [{'gene': 'CYP2C19',
                        'filter': ['NOMATCH'],
                        'geneQc': 'PASS'},
                       {'gene': 'CYP2D6',
                        'filter': ['NOMATCH', 'NOVEL'],
                        'geneQc': 'PASS'},
                       {'gene': 'VKORC1',
                        'filter': ['FAIL'],
                        'geneQc': 'FAIL'}
                      ]
    with open(qc_json) as fin:
        sample_qc = json.load(fin)

        # Check number of results are as expected and individual gene QC results agree.
        assert len(sample_qc['results']) == len(expected_output)
        for gene_qc_result, expected in zip(sample_qc['results'], expected_output):
            assert gene_qc_result['gene'] == expected['gene']
            assert gene_qc_result['filter'] == expected['filter']
            assert gene_qc_result['geneQc'] == expected['geneQc']

        # Check sample-level QC results.
        assert sample_qc['pgxQc'] == 'FAIL'
        assert sample_qc['failedGenes'] == ['VKORC1']

    # Spot-check star allele results in output JSON.
    output_json = '{}/{}_{}'.format(config['output']['samples_directory'], sample_id, config['suffixes']['final_json'])
    # The following expected output is a only a subset of the full output JSON contents.
    expected_output = [{'gene': 'CYP2C19',
                        'alleles': ['*1', '*37'],
                        'filter': ['NOMATCH']},
                       {'gene': 'CYP2D6',
                        'alleles': ['*29', '*29'],
                        'filter': ['NOMATCH', 'NOVEL']},
                      ]
    with open(output_json) as fin:
        sample_output = json.load(fin)
        # Check number of star allele results are as expected and that individual star allele results agree.
        assert len(sample_output['results']['starAlleles']) == len(expected_output)
        for gene_result, expected in zip(sample_output['results']['starAlleles'], expected_output):
            assert gene_result['gene'] == expected['gene']
            assert gene_result['alleles'] == expected['alleles']
            assert gene_result['filter'] == expected['filter']


def test_perform_failure_routine(batch2, config):
    perform_failure_routine(batch2, config)
    # Ensure the output and QC JSON files, which are ingested by Analysis-Workflow, is as expected.
    samples_dir = config['output']['samples_directory']
    star_gene_batches_unsorted = [v for v in batch2.gene_batches.values() if v.gene_config.result_type == 'star allele']
    star_gene_batches_sorted = sorted(star_gene_batches_unsorted, key=lambda x: x.gene_symbol)
    for sample_id in batch2.sample_ids:
        output_json = '{}/{}_{}'.format(samples_dir, sample_id, config['suffixes']['final_json'])
        with open(output_json) as out_json:
            sample_output = json.load(out_json)
        for i, star_gene_batch in enumerate(star_gene_batches_sorted):
            sample_output_results = sample_output['results']['starAlleles'][i]
            expected_output = {'gene': star_gene_batch.gene_symbol,
                               'alleles': None,
                               'callQuality': None,
                               'filter': None,
                               'degenerate': None,
                               'degenerateAlleles': None,
                               'copyNumber': None,
                               'variants': None}
            assert all(sample_output_results[k] == expected_output[k] for k in expected_output.keys())

    variant_gene_batches_unsorted = [v for v in batch2.gene_batches.values() if v.gene_config.result_type == 'variant']
    variant_gene_batches_sorted = sorted(variant_gene_batches_unsorted, key=lambda x: x.gene_symbol)
    for sample_id in batch2.sample_ids:
        output_json = '{}/{}_{}'.format(samples_dir, sample_id, config['suffixes']['final_json'])
        with open(output_json) as out_json:
            sample_output = json.load(out_json)
        for i, variant_gene_batch in enumerate(variant_gene_batches_sorted):
            sample_output_results = sample_output['results']['variantGenotypes'][i]
            expected_output = {'gene': variant_gene_batch.gene_symbol,
                               'genotypes': None,
                               'filter': None,
                              }
            assert all(sample_output_results[k] == expected_output[k] for k in expected_output.keys())


        qc_json = '{}/{}_{}'.format(samples_dir, sample_id, config['suffixes']['qc_json'])
        with open(qc_json) as qc_json:
            sample_qc = json.load(qc_json)
            assert sample_qc['pgxQc'] == 'FAIL'
            assert sample_qc['failedGenes'] == ['CYP2C19', 'CYP2D6', 'VKORC1']
            assert sample_qc['analysisStatus'] == 'Aborted'
            assert sample_qc['pgxReadCountQc'] == 'FAIL'

            for sample_qc_results in sample_qc['results']:
                expected_qc = {'gene': sample_qc_results['gene'],
                               'noCandidates': 'False',
                               'filter': None,
                               'geneQc': 'FAIL'}
                assert all(sample_qc_results[k] == expected_qc[k] for k in expected_qc.keys())

    # Ensure that missing assay version is handled appropriately
    batch_qc_json = '{}/{}_{}'.format(config['output']['directory'],
                                      batch2.batch_id,
                                      config['suffixes']['batch_qc_json'])
    with open(batch_qc_json) as json_file:
        batch_qc = json.load(json_file)

    assert batch_qc['assayVersion'] is None


def test_compress_sample_files(tmpdir):
    test_dir = str(tmpdir)
    source_dir = 'tests/data/compress'
    batch_id = 'batch'
    archive_suffix = 'sample_files.tar.gz'

    # Files that are expected to be compressed in the tarball.
    compressed_file_names = ['sample_round1_likelihoods.txt', 'sample_vgraph.txt']
    # Files that are expected not be compressed in the tarball.
    uncompressed_file_names = ['batch_calls_summary.txt', 'sample_output.json', 'sample_qc.json']

    # First, check that compress_sample_files will do nothing if directory is empty.
    compress_sample_files(batch_id, test_dir, archive_suffix)
    assert os.listdir(test_dir) == []

    # Next, copy files to be compressed to testing location.
    for fname in os.listdir(source_dir):
        shutil.copy('{}/{}'.format(source_dir, fname), test_dir)

    # Run compress_sample_files() on directory which now actually contains files.
    compress_sample_files(batch_id, test_dir, archive_suffix)
    tarball_path = os.path.join('{}/{}_{}'.format(test_dir, batch_id, archive_suffix))

    # Check that the expected files are included in the tarball with their original copies removed.
    assert (sorted(subprocess.check_output(['tar', '-tf', tarball_path]).decode().splitlines())
            == sorted(['{}/{}'.format(test_dir.lstrip('/'), name) for name in compressed_file_names]))
    for fname in compressed_file_names:
        assert not os.path.exists('{}/{}'.format(test_dir, fname))
    # Check that the batch and sample JSON files remain.
    for fname in uncompressed_file_names:
        assert filecmp.cmp('{}/{}'.format(test_dir, fname),
                           '{}/{}'.format(source_dir, fname))
