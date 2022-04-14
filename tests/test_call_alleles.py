import pytest

from collections import OrderedDict

import pandas as pd
import numpy as np

from helix_app_pgx.call_alleles import hits_within_qual, allele_combo_to_numbers
from helix_app_pgx.call_alleles import get_allele_structure, get_defining_variants, get_combo_variant_ploidy, \
                                       get_combo_structure, get_mutation_statistics, hybrids_are_complementary_cyp2d6, \
                                       loop_thru_combos, get_possible_alleles, create_sample_output_variants, get_quality, \
                                       get_filter, calculate_likelihoods, combo_violates_constraints_cyp2d6, \
                                       combo_violates_constraints_cyp2c19, create_sample_output_star
from helix_app_pgx.read_data import read_allele_depths, read_variants
from helix_app_pgx.structures import Sample
from fixtures import config, cyp2d6_gene_config, cyp2c19_gene_config, gene_configs, variant_gene_config, cyp2c19_params, \
                     cyp2d6_params


def test_get_allele_structure(cyp2d6_gene_config):
    assert np.array_equal(get_allele_structure('deletion', cyp2d6_gene_config, 2, 1), np.array([0]*9+[1]*9))
    assert np.array_equal(get_allele_structure('*36', cyp2d6_gene_config, 1, 1), cyp2d6_gene_config.struct_dict['*36'])


def test_get_defining_variants(cyp2d6_gene_config):
    defining_variants = cyp2d6_gene_config.defining_variants
    deletion_allele = cyp2d6_gene_config.deletion_allele

    # *14A does not exist in the current version of the reportable range, but *112 does
    with pytest.raises(KeyError):
        get_defining_variants('*14A', defining_variants, deletion_allele)

    assert defining_variants.loc['*112'].equals(
        get_defining_variants('*112', defining_variants, deletion_allele))


def test_get_defining_variants_ext(cyp2d6_gene_config):
    defining_variants = cyp2d6_gene_config.defining_variants
    defining_variants_ext = cyp2d6_gene_config.defining_variants_ext
    deletion_allele = cyp2d6_gene_config.deletion_allele

    # *1A is not in defining_variants but is in defining_variants_ext
    with pytest.raises(KeyError):
        get_defining_variants('*1A', defining_variants, deletion_allele)

    assert defining_variants_ext.loc['*1A'].equals(
        get_defining_variants('*1A', defining_variants_ext, deletion_allele))


def test_hits_within_qual1():
    # equal first elements
    likelihoods = pd.read_csv('tests/data/round1_likelihoods.txt', sep='\t', index_col=0)
    top_combos = hits_within_qual(likelihoods, 0)
    assert top_combos[0] == "('*10A', '*4M', 'null', 'null', 1, 1, 2)"
    assert top_combos[1] == "('*1B', '*4A', 'null', 'null', 1, 1, 2)"
    assert top_combos[2] == "('*39', '*4J', 'null', 'null', 1, 1, 2)"
    assert len(top_combos) == 3


def test_hits_within_qual2():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods_lowq.txt', sep='\t', index_col=0)
    top_combos = hits_within_qual(likelihoods)
    assert top_combos[0] == "('*10A', '*4M', 'null', 'null', 1, 1, 2)"
    assert top_combos[1] == "('*1B', '*4A', 'null', 'null', 1, 1, 2)"
    assert top_combos[2] == "('*39', '*4J', 'null', 'null', 1, 1, 2)"
    assert top_combos[3] == "('*4D', '*4M', 'null', 'null', 1, 1, 2)"
    assert len(top_combos) == 4


def test_hits_within_qual3():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods_neginf.txt', sep='\t', index_col=0)
    top_combos = hits_within_qual(likelihoods)
    assert top_combos[0] == "('*5', '*5', 'null', 'null', 1, 1, 4)"
    assert top_combos[1] == "('*5', '*5', '*16', 'null', 1, 1, 3)"
    assert len(top_combos) == 2


def test_get_quality_filter1():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods.txt', sep='\t', index_col=0)
    assert get_quality(likelihoods) == 99


def test_get_quality2():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods_lowq.txt', sep='\t', index_col=0)
    assert get_quality(likelihoods) == 3


def test_get_quality3():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods_neg_inf.txt', sep='\t', index_col=0)
    assert get_quality(likelihoods) == 0


def test_get_filter1():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods.txt', sep='\t', index_col=0)
    quality = get_quality(likelihoods)
    assert get_filter(likelihoods, quality) == 'PASS'


def test_get_filter2():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods_lowq.txt', sep='\t', index_col=0)
    quality = get_quality(likelihoods)
    assert get_filter(likelihoods, quality) == 'LOWQ'


def test_get_filter3():
    likelihoods = pd.read_csv('tests/data/round1_likelihoods_neg_inf.txt', sep='\t', index_col=0)
    quality = get_quality(likelihoods)
    assert get_filter(likelihoods, quality) == 'NOMATCH'


def test_allele_combo_to_numbers():
    assert allele_combo_to_numbers(('*5', '*4J', '*4N', '*68', 1, 1, 2)) == [4, 4, 5, 68]
    assert allele_combo_to_numbers(('*3A', '*4A', 'null', 'null', 1, 1, 3)) == [3, 4]


def test_get_mutation_statistics1(cyp2d6_params):
    params = cyp2d6_params
    cyp2d6_gene_config = params['gene_config']

    muts_all_df, _ = read_allele_depths('tests/data/na12877_vgraph.txt', cyp2d6_gene_config.mutations_lookup)
    muts_df = muts_all_df.loc[cyp2d6_gene_config.defining_variants.columns]

    putative_solution = ('*4D', '*4A', '*68', 'null', 1, 1, 2)
    dm = get_combo_variant_ploidy(putative_solution,
                                  candidate_structure=get_combo_structure(putative_solution, cyp2d6_gene_config),
                                  gene_config=cyp2d6_gene_config)
    muts_df['n_alts'] = dm['n_alt']
    muts_df['ploidies'] = dm['ploidy']

    mutation_likelihood, mutation_pvalue, low_coverage, mutation_clash, mutation_txids = get_mutation_statistics(muts_df, params)

    assert not any(low_coverage)
    assert any(mutation_clash)


def test_get_mutation_statistics2(cyp2d6_params):
    # correct solutions
    params = cyp2d6_params
    cyp2d6_gene_config = params['gene_config']

    muts_all_df, _ = read_allele_depths('tests/data/high_1846GtoA_ref_error_vgraph.txt', cyp2d6_gene_config.mutations_lookup)
    muts_df = muts_all_df.loc[cyp2d6_gene_config.defining_variants.columns]

    putative_solution = ('*4A', '*4A', '*68', '*68', 1, 1, 2)
    dm = get_combo_variant_ploidy(putative_solution,
                                  candidate_structure=get_combo_structure(putative_solution, cyp2d6_gene_config),
                                  gene_config=cyp2d6_gene_config)
    muts_df['n_alts'] = dm['n_alt']
    muts_df['ploidies'] = dm['ploidy']

    mutation_likelihood, mutation_pvalue, low_coverage, mutation_clash, mutation_txids = get_mutation_statistics(muts_df, params)

    assert not any(low_coverage)
    assert not any(mutation_clash)


def test_get_mutation_statistics3(cyp2d6_params):
    # require extremely high coverage and pvalue
    params = cyp2d6_params
    cyp2d6_gene_config = params['gene_config']

    # Override defaults for this test
    params['min_dp_per_copy'] = 50
    params['nomatch_pvalue'] = 1.e-2

    muts_all_df, _ = read_allele_depths('tests/data/high_1846GtoA_ref_error_vgraph.txt', cyp2d6_gene_config.mutations_lookup)
    muts_df = muts_all_df.loc[cyp2d6_gene_config.defining_variants.columns]

    putative_solution = ('*4A', '*4A', '*68', '*68', 1, 1, 2)
    dm = get_combo_variant_ploidy(putative_solution,
                                  candidate_structure=get_combo_structure(putative_solution, cyp2d6_gene_config),
                                  gene_config=cyp2d6_gene_config)
    muts_df['n_alts'] = dm['n_alt']
    muts_df['ploidies'] = dm['ploidy']

    mutation_likelihood, mutation_pvalue, low_coverage, mutation_clash, mutation_txids = get_mutation_statistics(muts_df, params)

    assert any(low_coverage)
    assert any(mutation_clash)


def test_get_mutation_statistics4(cyp2d6_params):
    # Homozygous deletion should not get NOMATCH filter
    params = cyp2d6_params
    cyp2d6_gene_config = params['gene_config']

    muts_all_df, _ = read_allele_depths('tests/data/homdel_vgraph.txt', cyp2d6_gene_config.mutations_lookup)
    muts_df = muts_all_df.loc[cyp2d6_gene_config.defining_variants.columns]

    putative_solution = ('*5', '*5', 'null', 'null', 1, 1, 2)
    dm = get_combo_variant_ploidy(putative_solution,
                                  candidate_structure=get_combo_structure(putative_solution, cyp2d6_gene_config),
                                  gene_config=cyp2d6_gene_config)
    muts_df['n_alts'] = dm['n_alt']
    muts_df['ploidies'] = dm['ploidy']

    mutation_likelihood, mutation_pvalue, low_coverage, mutation_clash, mutation_txids = get_mutation_statistics(muts_df, params)

    assert not any(mutation_clash)


def test_hybrids_are_complementary_cyp2d6(cyp2d6_gene_config):
    assert hybrids_are_complementary_cyp2d6('*78', '*68', cyp2d6_gene_config) is True
    assert hybrids_are_complementary_cyp2d6('*36', '*80', cyp2d6_gene_config) is True
    assert hybrids_are_complementary_cyp2d6('*13', 'null', cyp2d6_gene_config) is False
    assert hybrids_are_complementary_cyp2d6('*4N', '*57', cyp2d6_gene_config) is False


def test_combo_violates_constraints_cyp2d6(cyp2d6_gene_config):
    assert combo_violates_constraints_cyp2d6('*1A', '*4A', '*68', 'null', cyp2d6_gene_config) is False
    assert combo_violates_constraints_cyp2d6('*1A', '*4A', '*68', '*16', cyp2d6_gene_config) is True
    assert combo_violates_constraints_cyp2d6('*1A', '*4A', '*78', 'null', cyp2d6_gene_config) is True
    assert combo_violates_constraints_cyp2d6('*1A', '*4A', '*16', 'null', cyp2d6_gene_config) is True
    assert combo_violates_constraints_cyp2d6('*5', '*4A', '*16', 'null', cyp2d6_gene_config) is False
    assert combo_violates_constraints_cyp2d6('*1A', '*4A', '*13', 'null', cyp2d6_gene_config) is True


def test_combo_violates_constraints_cyp2c19(cyp2c19_gene_config):
    assert combo_violates_constraints_cyp2c19('CYP2C19*2.001', 'CYP2C19*11.001', 'CYP2C19*37.001', 'null',
                                              cyp2c19_gene_config) is True
    assert combo_violates_constraints_cyp2c19('CYP2C19*36.001', 'CYP2C19*11.001', 'CYP2C19*37.001', 'null',
                                              cyp2c19_gene_config) is False


def test_loop_thru_combos1(cyp2d6_params):
    params = cyp2d6_params
    cyp2d6_gene_config = params['gene_config']

    sample_emit_df = pd.read_csv('tests/data/star_13_exon_copy_emit.txt', sep='\t')
    sample_emit_df.set_index('exon', inplace=True)

    muts_all_df, _ = read_allele_depths('tests/data/star_13_vgraph.txt', cyp2d6_gene_config.mutations_lookup)
    muts_defining_df = muts_all_df.loc[cyp2d6_gene_config.defining_variants.columns]

    possible_nonhybrid, possible_hybrid = get_possible_alleles(muts_defining_df, params)
    mle_exon_ploidies = sample_emit_df.loc[cyp2d6_gene_config.all_exons, 'mle_n'].astype(int)
    round1 = loop_thru_combos(possible_hybrid, possible_nonhybrid, sample_emit_df, mle_exon_ploidies, muts_defining_df,
                              params)
    assert round1.iloc[0]['allele_combo'] == ('*2A', '*2A', '*77', 'null', 1, 1, 1)


def test_loop_thru_combos2(cyp2d6_params):
    # When the structural log likelihood is -inf, ensure results are sorted by mutation log likelihood.
    params = cyp2d6_params
    cyp2d6_gene_config = params['gene_config']

    sample_emit_df = pd.read_csv('tests/data/neginf_exon_copy_emit.txt', sep='\t')
    sample_emit_df.set_index('exon', inplace=True)

    muts_all_df, _ = read_allele_depths('tests/data/neginf_vgraph.txt', cyp2d6_gene_config.mutations_lookup)
    muts_defining_df = muts_all_df.loc[cyp2d6_gene_config.defining_variants.columns]

    possible_nonhybrid, possible_hybrid = get_possible_alleles(muts_defining_df, params)
    mle_exon_ploidies = sample_emit_df.loc[cyp2d6_gene_config.all_exons, 'mle_n'].astype(int)
    round1 = loop_thru_combos(possible_hybrid, possible_nonhybrid, sample_emit_df, mle_exon_ploidies, muts_defining_df,
                              params)

    assert round1.iloc[0]['allele_combo'] == ('*5', '*5', 'null', 'null', 1, 1, 4)


def test_loop_thru_combos3(cyp2c19_params):
    # Tests for the case where CYP2C19*37.001, the partial deletion allele, is expected.
    params = cyp2c19_params
    cyp2c19_gene_config = params['gene_config']
    sample_emit_df = pd.read_csv('tests/data/CYP2C19_partial_deletion_exon_copy_emit.txt', sep='\t')
    sample_emit_df.set_index('exon', inplace=True)

    muts_all_df, _ = read_allele_depths('tests/data/CYP2C19_partial_deletion_vgraph.txt', cyp2c19_gene_config.mutations_lookup)
    muts_defining_df = muts_all_df.loc[cyp2c19_gene_config.defining_variants.columns]

    # Technically, the partial deletion allele is not a hybrid, allele, but to be consistent with the variable names
    # in call_alleles.py, we'll use 'hybrid' below to refer to a partial deletion.
    possible_nonhybrid, possible_hybrid = get_possible_alleles(muts_defining_df, params)

    mle_exon_ploidies = sample_emit_df.loc[cyp2c19_gene_config.all_exons, 'mle_n'].astype(int)
    round1 = loop_thru_combos(possible_hybrid, possible_nonhybrid, sample_emit_df, mle_exon_ploidies, muts_defining_df,
                              params)

    # In the tuple representation below, there is a *36 allele to indicate an absence of a full gene allele.
    assert round1.iloc[0]['allele_combo'] == ('CYP2C19*36.001', 'CYP2C19*1.002', 'CYP2C19*37.001', 'null', 1, 1, 0)


def test_loop_thru_combos4(cyp2c19_params):
    # Tests for the case where CYP2C19*37.001, the partial deletion allele, is not expected.
    # This is a check on the combo_violates_constraints_cyp2c19 function. Without that function,
    # the data files below would support an allele combination with *2, *11, and *37,
    # a partial deletion in tandem with a full gene allele, which is not empirically observed.

    params = cyp2c19_params
    cyp2c19_gene_config = params['gene_config']
    sample_emit_df = pd.read_csv('tests/data/CYP2C19_no_partial_deletion_exon_copy_emit.txt', sep='\t')
    sample_emit_df.set_index('exon', inplace=True)

    muts_all_df, _ = read_allele_depths('tests/data/CYP2C19_no_partial_deletion_vgraph.txt', cyp2c19_gene_config.mutations_lookup)
    muts_defining_df = muts_all_df.loc[cyp2c19_gene_config.defining_variants.columns]

    # Technically, the partial deletion allele is not a hybrid, allele, but to be consistent with the variable names
    # in call_alleles.py, we'll use 'hybrid' below to refer to a partial deletion.
    possible_nonhybrid, possible_hybrid = get_possible_alleles(muts_defining_df, params)

    mle_exon_ploidies = sample_emit_df.loc[cyp2c19_gene_config.all_exons, 'mle_n'].astype(int)
    round1 = loop_thru_combos(possible_hybrid, possible_nonhybrid, sample_emit_df, mle_exon_ploidies, muts_defining_df,
                              params)
    assert round1.iloc[0]['allele_combo'] == ('CYP2C19*2.001', 'CYP2C19*11.001', 'null', 'null', 1, 1, 0)


def test_loop_thru_combos5(cyp2c19_params):
    # Regression test to ensure bug CB-253 does not occur.
    # This test evaluates a sample with three copies of CYP2C19. Before the bug fix, it was possible for 
    # an allele combination that had 5 copies of some exons (4 full gene copies plus 1 partial gene) to be
    # considered a candidate allele combination. Since the maximum gene copy for CYP2C19 for calculating
    # exon copy likelihoods is 2 per chromosome, there were no emission probabilities available for the 
    # structure statistics calculation, leading to an error.
    params = cyp2c19_params
    cyp2c19_gene_config = params['gene_config']
    sample_emit_df = pd.read_csv('tests/data/CYP2C19_duplication_exon_copy_emit.txt', sep='\t')
    sample_emit_df.set_index('exon', inplace=True)

    muts_all_df, _ = read_allele_depths('tests/data/CYP2C19_duplication_vgraph.txt', cyp2c19_gene_config.mutations_lookup)
    muts_defining_df = muts_all_df.loc[cyp2c19_gene_config.defining_variants.columns]

    possible_nonhybrid, possible_hybrid = get_possible_alleles(muts_defining_df, params)

    mle_exon_ploidies = sample_emit_df.loc[cyp2c19_gene_config.all_exons, 'mle_n'].astype(int)
    round1 = loop_thru_combos(possible_hybrid, possible_nonhybrid, sample_emit_df, mle_exon_ploidies, muts_defining_df,
                              params)
    assert round1.iloc[0]['allele_combo'] == ('CYP2C19*2.001', 'CYP2C19*11.001', 'null', 'null', 1, 2, 0)


def test_calculate_likelihoods(cyp2c19_params):
    params = cyp2c19_params
    cyp2c19_gene_config = params['gene_config']

    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v3')

    muts_all_df, novel_muts = read_allele_depths('tests/data/CYP2C19_partial_deletion_vgraph.txt',
                                                 cyp2c19_gene_config.mutations_lookup)
    sample.data['muts_all_df'] = muts_all_df
    sample.data['novel_muts'] = novel_muts
    sample.data['exon_emit_mat'] = pd.read_csv('tests/data/CYP2C19_partial_deletion_exon_copy_emit.txt',
                                                sep='\t')
    sample_results = calculate_likelihoods(sample, params, cyp2c19_gene_config.scratch_dir)

    # In the tuple representation below, there is a *36 allele to indicate an absence of a full gene allele.
    assert sample_results['allele_call'] == ('CYP2C19*36.001', 'CYP2C19*1.002', 'CYP2C19*37.001', 'null', 1, 1, 0)


def test_create_sample_output_star(cyp2d6_gene_config, cyp2d6_params):
    # First, populate sample.results using test data and by calling calculate_likelihoods.
    params = cyp2d6_params
    cyp2d6_gene_config = params['gene_config']

    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v3')

    muts_all_df, novel_muts = read_allele_depths('tests/data/novel_nomatch_vgraph.txt',
                                                 cyp2d6_gene_config.mutations_lookup)
    sample.data['muts_all_df'] = muts_all_df
    sample.data['novel_muts'] = novel_muts
    sample.data['exon_emit_mat'] = pd.read_csv('tests/data/novel_nomatch_exon_copy_emit.txt',
                                                sep='\t')
    sample.results = calculate_likelihoods(sample, params, cyp2d6_gene_config.scratch_dir)

    # Modify sample.output in place.
    create_sample_output_star(sample, cyp2d6_gene_config)

    expected_novel_variant = ('chr22\t42130638\t.\tG\tA\t2707.77\t.\tAC=1;AF=0.5;AN=2;BaseQRankSum=-2.134;ClippingRankSum=0;DP=177;'
                              'ExcessHet=3.0103;FS=1.19;MLEAC=1;MLEAF=0.5;MQ=57.94;MQRankSum=1.143;QD=15.3;ReadPosRankSum=-2.243;'
                              'SOR=0.792;ANN=A|stop_gained|HIGH|CYP2D6|ENSG00000100197.21|transcript|ENST00000360608.9|'
                              'protein_coding|1/9|c.154C>T|p.Gln52*|269/1684|154/1494|52/497||;LOF=(CYP2D6|ENSG00000100197.21|'
                              '1|1.00);NMD=(CYP2D6|ENSG00000100197.21|1|1.00)\tGT:AD:DP:GQ:PL\t0/1:96,81:177:99:2736,0,3663')

    # Check that single-line fields are as expected.
    assert sample.output['alleles'] == ['*29', '*29']
    assert sample.output['callQuality'] == 99
    assert sample.output['filter'] == ['NOMATCH', 'NOVEL']
    assert sample.output['degenerate'] is False
    assert sample.output['degenerateAlleles'] is None
    assert sample.output['novel'] == [expected_novel_variant]

    # Spot-check first element of copyNumber value.
    expected_first_copy_number = {"exonId": "CYP2D6.e1",
                                  "callPloidy": 2,
                                  "maxLikelihoodPloidy": 2,
                                  "ploidyPval": "3.77e-01"
                                 }
    for k,v in sample.output['copyNumber'][0].items():
        assert v == expected_first_copy_number[k]

    # Spot-check first element of variants value.
    expected_first_variant = {"variantId": "g.42126605G>A",
                              "rsId": "rs568495591",
                              "txId": "4187C>T",
                              "refPloidy": 2,
                              "altPloidy": 0,
                              "refAd": 49,
                              "altAd": 0,
                             "mutPval": "6.11e-01"
                             }
    for k,v in sample.output['variants'][0].items():
        assert v == expected_first_variant[k]


def test_create_sample_output_variants(variant_gene_config):
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v2')
    sample.data['r2v_calls'] = read_variants('tests/data/vkorc1_variants.txt')
    create_sample_output_variants(sample, variant_gene_config)

    expected_genotypes = [OrderedDict([('rsId', 'rs9923231'),
                                       ('genomicId', ['g.31096368C>A', 'g.31096368C>G', 'g.31096368C>T']),
                                       ('transcriptId', None),
                                       ('genomicGenotype', ['C', 'T']),
                                       ('transcriptGenotype', None),
                                       ('filter', ['PASS']),
                                       ('genotypeQuality', 99),
                                       ('genomicRefAd', 37),
                                       ('genomicAltAd', [0, 0, 16])]),
                          ]
    expected_result = OrderedDict([('gene', 'VKORC1'),
                                  ('genotypes', expected_genotypes),
                                  ('filter', ['PASS'])])
    assert sample.output == expected_result


def test_create_sample_output_variants_no_rsid(variant_gene_config):
    sample = Sample('test_sample', '123456', 'AN-1', 'test.bam', 'test.bam.bai',
                    'test.bcf', 'test.bcf.csi', 'qctable.txt',
                    'qctable_by_region.txt',
                    'v2')
    sample.data['r2v_calls'] = read_variants('tests/data/vkorc1_variants_no_rsid.txt')
    create_sample_output_variants(sample, variant_gene_config)

    expected_genotypes = [OrderedDict([('rsId', 'rs9923231'),
                                       ('genomicId', ['g.31096368C>A', 'g.31096368C>G', 'g.31096368C>T']),
                                       ('transcriptId', None),
                                       ('genomicGenotype', None),
                                       ('transcriptGenotype', None),
                                       ('filter', ['LOWQ']),
                                       ('genotypeQuality', None),
                                       ('genomicRefAd', None),
                                       ('genomicAltAd', None)]),
                          ]
    expected_result = OrderedDict([('gene', 'VKORC1'),
                                  ('genotypes', expected_genotypes),
                                  ('filter', ['FAIL'])])
    assert sample.output == expected_result
