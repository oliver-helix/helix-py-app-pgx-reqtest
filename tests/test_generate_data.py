import os
import tempfile
import filecmp
from pathlib import Path
from os.path import join as path_join

import pysam
import pandas as pd

from fixtures import config, gene_configs, variant_gene_config, cyp2d6_gene_config
from helix_app_pgx.settings import ALIGNMENT_EXCLUSION_MASK
from helix_app_pgx.generate_data import (reassign_counts_cyp2d6, get_bcf_slice, parse_vcf, pseudo_post_alt, remap,
                                         select_readnames_by_pileup, rule_based_filtering_exon7, 
                                         rule_based_filtering_exon6, rule_based_filtering_exon4, 
                                         rule_based_filtering_exon2, is_left_mate, get_chromosome_slice)


def test_get_chromosome_slice():
    test_file = tempfile.NamedTemporaryFile()
    params = {
        "bam_url": "tests/data/NA12878_slice.bam",
        "bai_url": "tests/data/NA12878_slice.bam.bai",
        "chromosome": "chr22",
        "bam_flag": 1, 
        "bam_slice": test_file.name, 
        "ref_fasta": None, 
    }
    get_chromosome_slice(**params)  # copy slice from bam_url to bam_slice
    with pysam.AlignmentFile(test_file.name) as test_alignment, \
        pysam.AlignmentFile(params["bam_url"]) as correct_alignment:
        test_names = []
        for aligned_segment in test_alignment.fetch(params["chromosome"], 42126400, 42126500):
            if not (aligned_segment.flag & ALIGNMENT_EXCLUSION_MASK):
                test_names.append(aligned_segment.to_dict()["name"])
        correct_names = []
        for aligned_segment in correct_alignment.fetch(params["chromosome"], 42126400, 42126500):
            if not (aligned_segment.flag & ALIGNMENT_EXCLUSION_MASK):
                correct_names.append(aligned_segment.to_dict()["name"])

    assert test_names == correct_names

    os.remove(test_file.name + ".bai")
    test_file.close()

def select_readnames_by_region(bamfile, contig, start, end):
    readnames = set()
    with pysam.AlignmentFile(bamfile) as bam_in:
        for read in bam_in.fetch(contig, start, end):
            readnames.add(read.query_name)
    return readnames

def test_pseudo_post_alt_exon1(cyp2d6_gene_config):
    # This sample is *2/*35, where both haplotypes significantly differ from the ref genome in CYP2D6 intron1. This results in
    # many reads being mismapped to CYP2D7. The pseudo postalt step should reassign them to CYP2D6.
    gene_slice_bam = 'tests/data/CL-68-Plate_1-F01_133159728_gene_slice.bam'
    remap_ref = cyp2d6_gene_config.contig_mod2
    alt_ref = cyp2d6_gene_config.contig_alt
    scratch = cyp2d6_gene_config.scratch_dir
    prefix= 'CL-68-Plate_1-F01_133159728'
    readgroup = '{}{}'.format('@RG\tID:0\tPG:r2v\tPM:HiSeq\tSM:', prefix)
    primary_bam = path_join(scratch, 'CL-68-Plate_1-F01_133159728_primary.bam')
    alt_bam = path_join(scratch, 'CL-68-Plate_1-F01_133159728_alt.bam')

    remap(gene_slice_bam, remap_ref, readgroup, primary_bam, scratch, prefix)
    remap(gene_slice_bam, alt_ref, readgroup, alt_bam, scratch, prefix)
    pysam.index(alt_bam)
    reassign_to_cyp2d6_reads, reassign_to_cyp2d7_reads = pseudo_post_alt(primary_bam, alt_bam)

    # Check that a large number of reads in intron1, near exon junction (exon1 starts at 52953), are reassigned.
    alt_exon1_reads = select_readnames_by_region(alt_bam, 'chr22_KI270928v1_alt', 52893, 52953)
    assert len(reassign_to_cyp2d6_reads.intersection(alt_exon1_reads)) == 60
    assert 'E00520:268:hcj7kcfxx:1:2209:16967:23802' in reassign_to_cyp2d6_reads
    assert 'E00520:268:hcj7kcfxx:3:1216:14495:40013' in reassign_to_cyp2d6_reads

def test_pseudo_post_alt_exon7(cyp2d6_gene_config):
    # this sample is *5/*5, but with significant mismapping from CYP2D7 to exno7
    gene_slice_bam = 'tests/data/SA_homozygous_star5_153023872_gene_slice.bam'
    remap_ref = cyp2d6_gene_config.contig_mod2
    alt_ref = cyp2d6_gene_config.contig_alt
    scratch = cyp2d6_gene_config.scratch_dir
    prefix= 'SA_homozygous_star5_153023872'
    readgroup = '{}{}'.format('@RG\tID:0\tPG:r2v\tPM:HiSeq\tSM:', prefix)
    primary_bam = path_join(scratch, 'SA_homozygous_star5_153023872_primary.bam')
    alt_bam = path_join(scratch,'SA_homozygous_star5_153023872_alt.bam')

    remap(gene_slice_bam, remap_ref, readgroup, primary_bam, scratch, prefix)
    remap(gene_slice_bam, alt_ref, readgroup, alt_bam, scratch, prefix)
    pysam.index(alt_bam)
    reassign_to_cyp2d6_reads, reassign_to_cyp2d7_reads = pseudo_post_alt(primary_bam, alt_bam)

    # check a large fraction(46/74) of the original exon7 reads get reassigned
    primary_exon7_reads = select_readnames_by_region(primary_bam, 'chr22', 42127446, 42127634)
    assert len(reassign_to_cyp2d7_reads.intersection(primary_exon7_reads)) == 46

def test_select_readnames_by_pileup():
    # note that IGV somehow shows 9 reads supporting the ref allele
    primary_bam = 'tests/data/NA12877_exon4_subsample.bam'
    star4_ref_dict = {42128944: 'C'}

    star4_ref_reads = select_readnames_by_pileup(primary_bam, star4_ref_dict, 'chr22', 42128943, 42128945)
    assert len(star4_ref_reads) == 7

def test_rule_based_filtering_exon7():
    # this sample is *4/*4, the CYP2D7 sequence is closer to the CYP2D6 refgenome sequence than to the CYP2D7 refgenome
    # sequence, causing significant mismapping to CYP2D6
    primary_bam = 'tests/data/CL-10-Plate_1-B06_133615928_primary.bam'

    mismapped_reads = rule_based_filtering_exon7(primary_bam)

    # verify a large number of exon7 reads are identified as mismapped
    exon7_reads = select_readnames_by_region(primary_bam, 'chr22', 42127446, 42127634)
    assert len(mismapped_reads.intersection(exon7_reads)) == 84

def test_rule_based_filtering_exon6():
    #This SA sample(*15/*29) has low level of exon6 mismapping
    primary_bam = 'tests/data/SA_exon6_mismap_153025880_gene_slice.bam'

    mismapped_reads = rule_based_filtering_exon6(primary_bam)

    exon6_reads = select_readnames_by_region(primary_bam, 'chr22', 4212784, 42127983)
    assert len(mismapped_reads.intersection(exon6_reads)) == 21

def test_rule_based_filtering_exon4():
    #NA12877 *4/*4+*68, should be homozygous alt for *4 defining variant. All reads supporting ref are mismapped
    primary_bam = 'tests/data/NA12877-PC0000093-G05_153029878_primary_exon34.bam'

    mismapped_reads, ambiguous_reads = rule_based_filtering_exon4(primary_bam)

    # verify most of the star4 ref reads(78/94) are classified as either mismapped or ambiguous
    star4_ref_dict = {42128944: 'C'}
    star4_ref_reads = select_readnames_by_pileup(primary_bam, star4_ref_dict, 'chr22', 42128943, 42128945)
    assert len(mismapped_reads | ambiguous_reads) == 78

def test_rule_based_filtering_exon4_unmapped():
    #A read is unmapped and has no CIGAR string
    primary_bam = 'tests/data/ST_NA11892_exon4_unmapped.bam'

    mismapped_reads, ambiguous_reads = rule_based_filtering_exon4(primary_bam)

    # make sure the program doesn't die
    assert len(mismapped_reads) == 0
    assert len(ambiguous_reads) == 0

def test_rule_based_filtering_exon2():
    #A particular haplotype of CYP2D7 has intron 2 sequence(close to exon2) more similar to CYP2D6, causing mismapping
    primary_bam = 'tests/data/SA_exon2_mismap_145811723_primary.bam'

    mismapped_reads = rule_based_filtering_exon2(primary_bam)

    exon2_reads = select_readnames_by_region(primary_bam, 'chr22', 42129738, 42129909)
    assert len(mismapped_reads.intersection(exon2_reads)) == 13

def test_reassign_counts_cyp2d6():
    exon_counts_bed = 'tests/data/exon_counts.bed'
    exon_mismap_bed = 'tests/data/exon_mismap.bed'
    exon_alt_bed = 'tests/data/exon_alt.bed'
    test_output = 'tests/output/exon_corrected.bed'

    reassign_counts_cyp2d6(exon_counts_bed, exon_mismap_bed, exon_alt_bed, test_output)

    correct_output = 'tests/data/exon_corrected.bed'

    names = ['chr', 'start', 'end', 'label', 'count']
    test_df = pd.read_csv(test_output, names=names, header=None, sep='\t', index_col=3)
    correct_df = pd.read_csv(correct_output, names=names, header=None, sep='\t', index_col=3)

    assert test_df.equals(correct_df)

def test_get_bcf_slice(variant_gene_config):
    bcf_slice = 'tests/data/vkorc1_test_slice.bcf'
    bcf_slice_index = bcf_slice + '.csi'
    targets = variant_gene_config.gene_targets
    vcf_out = 'tests/output/vkorc1_noblocks_test.vcf'
    get_bcf_slice(targets, bcf_slice, bcf_slice_index, vcf_out)

    expected_vcf_out = 'tests/data/vkorc1_noblocks.vcf'
    assert filecmp.cmp(vcf_out, expected_vcf_out)

    Path.unlink(Path(vcf_out))

def test_parse_vcf(variant_gene_config):
    vcf = 'tests/data/vkorc1_noblocks.vcf'
    variants_txt = '/tests/scratch/vkorc1_variants_test.txt'
    variants_subset = variant_gene_config.variants_subset
    rsid2genomic = variant_gene_config.rsid2genomic
    parse_vcf(vcf, variants_txt, variants_subset, rsid2genomic)

    expected_variants_txt = 'tests/data/vkorc1_variants.txt'

    assert filecmp.cmp(variants_txt, expected_variants_txt)
    Path.unlink(Path(variants_txt))


def test_parse_vcf_imp(variant_gene_config):
    vcf = 'tests/data/vkorc1_imp.vcf'
    variants_txt = '/tests/scratch/vkorc1_variants_test_imp.txt'
    variants_subset = variant_gene_config.variants_subset
    rsid2genomic = variant_gene_config.rsid2genomic
    parse_vcf(vcf, variants_txt, variants_subset, rsid2genomic)

    expected_variants_txt = 'tests/data/vkorc1_variants_imp.txt'

    assert filecmp.cmp(variants_txt, expected_variants_txt)


def test_parse_vcf_nocall(variant_gene_config):
    vcf = 'tests/data/vkorc1_nocall.vcf'
    variants_txt = '/tests/scratch/vkorc1_variants_test_nocall.txt'
    variants_subset = variant_gene_config.variants_subset
    rsid2genomic = variant_gene_config.rsid2genomic
    parse_vcf(vcf, variants_txt, variants_subset, rsid2genomic)

    expected_variants_txt = 'tests/data/vkorc1_variants_nocall.txt'

    assert filecmp.cmp(variants_txt, expected_variants_txt)


def test_parse_vcf_norsid(variant_gene_config):
    vcf = 'tests/data/vkorc1_no_rsid.vcf'
    variants_txt = 'tests/output/vkorc1_variants_no_rsid_test.txt'
    variants_subset = variant_gene_config.variants_subset
    rsid2genomic = variant_gene_config.rsid2genomic
    parse_vcf(vcf, variants_txt, variants_subset, rsid2genomic)

    expected_variants_txt = 'tests/data/vkorc1_variants_no_rsid.txt'

    assert filecmp.cmp(variants_txt, expected_variants_txt)
    Path.unlink(Path(variants_txt))
