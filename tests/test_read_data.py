from helix_app_pgx.read_data import read_allele_depths, read_qc, read_qc_by_region, read_variants

from fixtures import config, cyp2d6_gene_config, gene_configs


def test_read_allele_depths1(cyp2d6_gene_config):
    muts_all_df, _ = read_allele_depths('tests/data/na12877_vgraph.txt', cyp2d6_gene_config.mutations_lookup)
    muts_df = muts_all_df.loc[cyp2d6_gene_config.defining_variants.columns]

    assert muts_df.loc['4181G>C','alt'] == 122

def test_read_allele_depths2(cyp2d6_gene_config):
    muts_all_df, novel_muts = read_allele_depths('tests/data/sa_sample_novel_lof_vgraph.txt',
                                                 cyp2d6_gene_config.mutations_lookup)
    vcf_fields = novel_muts[0].split('\t')
    assert vcf_fields[0] == 'chr22'
    assert len(vcf_fields)== 10

def test_read_qc():
    r2v_qc = read_qc('tests/data/sample_qctable.txt')
    assert r2v_qc['AT_DROPOUT'] == '23.749874'

def test_read_qc_by_region():
    r2v_qc = read_qc_by_region('tests/data/sample_qctable-by-region.txt')
    assert r2v_qc['CODING_PCT_OFF_BAIT'] == '0.403059'
    assert r2v_qc['MENDELIOMECORE_PCT_TARGET_BASES_20X'] == '0.970846'

def test_read_qc_by_region2():
    r2v_qc = read_qc_by_region('tests/data/sample_qctable-by-region_oldr2v.txt')
    assert r2v_qc['MENDELIOMECORE_PCT_TARGET_BASES_20X'] == '0.977512'

def test_read_variants():
    r2v_calls = read_variants('tests/data/vkorc1_variants.txt')
    assert eval(r2v_calls['rs9923231']['filter']) == ['PASS']
    assert eval(r2v_calls['rs9923231']['genomic_genotype']) == ['C','T']
