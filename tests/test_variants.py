import os

import pandas as pd
from pysam import VariantFile

from helix_app_pgx.variants import annotate_with_snpeff, parse_vgraph_out_cyp2d6, merge_pgx_variant_records, annotate_with_vgraph
from helix_lib_core.convert import tryfloat
from fixtures import config, cyp2d6_gene_config, cyp2c19_gene_config, gene_configs


def test_annotate_with_snpeff_cyp2d6(config):
    # test annotate novel LOF for CYP2D6
    input_vcf = '/src/tests/data/CYP2D6_novel_lof.vcf.gz'
    snpeff_vcf = os.path.join(config['scratch']['directory'], 'CYP2D6_novel_lof_snpeff.vcf.gz')

    annotate_with_snpeff(input_vcf, snpeff_vcf, config)

    rec = next(VariantFile(snpeff_vcf))
    assert 'LOF' in rec.info


def test_annotate_with_snpeff_cyp4f2(config):
    # test annotate novel LOF for CYP4F2
    input_vcf = '/src/tests/data/CYP4F2_novel_lof.vcf.gz'
    snpeff_vcf = os.path.join(config['scratch']['directory'], 'CYP4F2_novel_lof_snpeff.vcf.gz')

    annotate_with_snpeff(input_vcf, snpeff_vcf, config)

    rec = next(VariantFile(snpeff_vcf))
    assert 'LOF' in rec.info


def test_annotate_with_snpeff_cyp2c9(config):
    # test annotate novel LOF for CYP2C9
    input_vcf = '/src/tests/data/CYP2C9_novel_lof.vcf.gz'
    snpeff_vcf = os.path.join(config['scratch']['directory'], 'CYP2C9_novel_lof_snpeff.vcf.gz')

    annotate_with_snpeff(input_vcf, snpeff_vcf, config)

    rec = next(VariantFile(snpeff_vcf))
    assert 'LOF' in rec.info


def test_annotate_with_snpeff_cyp2c19(config):
    # test annotate novel LOF for CYP2C19
    input_vcf = '/src/tests/data/CYP2C19_novel_lof.vcf.gz'
    snpeff_vcf = os.path.join(config['scratch']['directory'], 'CYP2C19_novel_lof_snpeff.vcf.gz')

    annotate_with_snpeff(input_vcf, snpeff_vcf, config)

    rec = next(VariantFile(snpeff_vcf))
    assert 'LOF' in rec.info


def test_parse_vgraph_out_cyp2d6(config, cyp2d6_gene_config):
    annotated_vcf = '/src/tests/data/na12878_lof_vgraph_annotated.vcf'
    variants_txt = os.path.join(config['scratch']['directory'], 'na12878_vgraph_out.txt')
    gene_ref_vcf = cyp2d6_gene_config.gene_ref_vcf
    het_range_normal = [tryfloat(config['parameters']['het_normal_lower']),
                        tryfloat(config['parameters']['het_normal_upper'])]
    het_range_exon2 = [tryfloat(config['parameters']['het_exon2_lower']),
                       tryfloat(config['parameters']['het_exon2_upper'])]
    parse_vgraph_out_cyp2d6(annotated_vcf, variants_txt, gene_ref_vcf, het_range_normal, het_range_exon2)

    df = pd.read_csv(variants_txt, sep='\t')
    assert df.shape[0] == 117 # number of records in gene_ref_vcf
    assert all(df.loc[df['mutation'].isin(['g.42128945C>T', 'g.42128242delT'])]['call']==1)
    assert sum(df['call']) == 2


def test_parse_novel_vgraph_out_cyp2d6(config, cyp2d6_gene_config):
    annotated_vcf = '/src/tests/data/example_novel_lof_vgraph_annotated.vcf'
    variants_txt = os.path.join(config['scratch']['directory'], 'novel_vgraph_out.txt')
    gene_ref_vcf = cyp2d6_gene_config.gene_ref_vcf
    het_range_normal = [tryfloat(config['parameters']['het_normal_lower']),
                        tryfloat(config['parameters']['het_normal_upper'])]
    het_range_exon2 = [tryfloat(config['parameters']['het_exon2_lower']),
                       tryfloat(config['parameters']['het_exon2_upper'])]
    parse_vgraph_out_cyp2d6(annotated_vcf, variants_txt, gene_ref_vcf, het_range_normal, het_range_exon2)

    df = pd.read_csv(variants_txt, sep='\t')
    assert df.shape[0] == 118 # number of records in gene_ref_vcf plus 1
    assert df.iloc[117]['mutation'].startswith('Novel:')


def test_parse_vgraph_out_multi_alt(config, cyp2d6_gene_config):
    annotated_vcf = 'tests/data/example_multi_alt_vgraph_annotated.vcf'
    variants_txt = os.path.join(config['scratch']['directory'], 'multi_alt_vgraph_out.txt')
    expected_variants_txt = 'tests/output/example_multi_alt_vgraph_out.txt'
    gene_ref_vcf = cyp2d6_gene_config.gene_ref_vcf
    het_range_normal = [tryfloat(config['parameters']['het_normal_lower']), tryfloat(config['parameters']['het_normal_upper'])]
    het_range_exon2 = [tryfloat(config['parameters']['het_exon2_lower']), tryfloat(config['parameters']['het_exon2_upper'])]
    parse_vgraph_out_cyp2d6(annotated_vcf, variants_txt, gene_ref_vcf, het_range_normal, het_range_exon2)

    df = pd.read_csv(variants_txt, sep='\t')
    expected_df = pd.read_csv(expected_variants_txt, sep='\t')
    assert df.equals(expected_df)


def test_merge_pgx_variant_records_empty_vcf(config):
    contig = 'chr10'
    contig_ref = '/src/reference/GCA_000001405.15_GRCh38_full-hs38d1-rsrs-homd100_analysis_set_chr10.fna'
    variants_gvcf = 'tests/data/variants_all_ref.g.vcf.gz'
    variants_vcf = 'tests/data/variants_all_ref.vcf.gz'
    merged_vcf =  os.path.join(config['scratch']['directory'], 'variants_all_ref_merged_test.vcf.gz')
    expected_merged_vcf = 'tests/data/merged_all_ref.vcf.gz'

    merge_pgx_variant_records(contig, contig_ref, variants_gvcf, variants_vcf, merged_vcf)

    with VariantFile(merged_vcf) as test, VariantFile(expected_merged_vcf) as expected:
        for test_record, expected_record in zip(test, expected):
            assert test_record == expected_record


def test_annotate_with_vgraph_cyp2d6(config, cyp2d6_gene_config, tmpdir):
    test_dir = str(tmpdir)
    input_vcf = 'tests/data/cyp2d6_adj_phased_variants_snpeff.vcf.gz'
    annotated_vcf = f'{test_dir}/cyp2d6_annotated_vgraph.vcf'

    annotate_with_vgraph(input_vcf, annotated_vcf, config, cyp2d6_gene_config)
    
    variants_found = []
    with open(annotated_vcf) as anno_vcf:
        for record in anno_vcf:
            if record.startswith('#'):
                continue
            fields = record.strip().split('\t')
            if fields[0] == 'chr22' and fields[1] == '42129132':
                if 'HGVSC_FOUND' in record and 'g.42129132C>T' in record:
                    variants_found.append(record)

    assert len(variants_found) == 1


def test_annotate_with_vgraph_cyp2c19(config, cyp2c19_gene_config, tmpdir):
    test_dir = str(tmpdir)
    input_vcf = 'tests/data/cyp2c19_adj_phased_variants_snpeff.vcf.gz'
    annotated_vcf = f'{test_dir}/cyp2c19_annotated_vgraph.vcf'

    annotate_with_vgraph(input_vcf, annotated_vcf, config, cyp2c19_gene_config)
    
    variants_found = []
    with open(annotated_vcf) as anno_vcf:
        for record in anno_vcf:
            if record.startswith('#'):
                continue
            fields = record.strip().split('\t')
            if fields[0] == 'chr10' and fields[1] == '94842866':
                if 'HGVSC_FOUND' in record and 'g.94842866A>G' in record:
                    variants_found.append(record)

    assert len(variants_found) == 1
