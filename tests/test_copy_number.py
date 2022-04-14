from helix_app_pgx.copy_number import call_whole_gene_copy, call_exon_copy, stepwise_poisson,\
                                      adjust_for_covariate
from fixtures import config, batch1, gene_configs


def test_adjust_for_covariate(batch1):
    gene_batch = batch1.gene_batches['CYP2D6']
    counts_df = gene_batch.counts_df
    baseline_sample_ids = counts_df.loc[counts_df['nearest_center']=='2'].index.tolist()
    adjusted = adjust_for_covariate(counts_df, 'gene_ratio', 'MEDIAN_INSERT_SIZE', baseline_sample_ids, 0.01)
    assert sum(adjusted - counts_df['gene_ratio']) == 0
    adjusted = adjust_for_covariate(counts_df, 'gene_ratio', 'AT_DROPOUT', baseline_sample_ids, 0.01)
    assert sum(adjusted - counts_df['gene_ratio']) != 0
    assert round(adjusted.loc[baseline_sample_ids].mean(), 8) == round(counts_df.loc[baseline_sample_ids, 'gene_ratio'].mean(), 8)


def test_call_whole_gene_copy(batch1, config):
    gene_batch = batch1.gene_batches['CYP2D6']
    call_whole_gene_copy(gene_batch, batch1.batch_id, config)
    assert gene_batch.qc['number_baseline'] == 51
    # This sample should be in baseline, after we tweak the algo
    assert 'CAP-2-Plate_1-E07_133615891' in gene_batch.baseline_sample_ids


def test_stepwise_poisson(batch1):
    gene_batch = batch1.gene_batches['CYP2D6']
    two_copy_counts = gene_batch.counts_df.loc[gene_batch.counts_df['nearest_center']=='2']
    two_copy_counts.loc[two_copy_counts.index, 'Intercept'] = 1
    covariates = ['Intercept', 'chr22', 'AT_DROPOUT', 'MEDIAN_INSERT_SIZE']

    fit, covariates = stepwise_poisson(two_copy_counts, 'CYP2D6.e1', covariates, 0.01, 'chr22')
    assert covariates == ['Intercept', 'chr22']
    assert all([ x<=0.01 for x in fit.pvalues])


def test_call_exon_copy(batch1, config):
    gene_batch = batch1.gene_batches['CYP2D6']
    counts_df = gene_batch.counts_df
    gene_batch.baseline_sample_ids = counts_df.loc[counts_df['nearest_center']=='2'].index.tolist()
    call_exon_copy(gene_batch, batch1.batch_id, config)
    na12878_df = gene_batch.samples['NA12878-CL-01-Plate_1-A03_132964731'].data['exon_emit_mat']
    # this drops the sample_id index
    na12878_df.set_index('exon', inplace=True, drop=False)
    assert na12878_df.loc['CYP2D6.e1']['mle_n'] == '3'
    assert na12878_df.loc['CYP2D6.e2']['mle_n'] == '2'
    assert na12878_df.loc['CYP2D7.e1']['mle_n'] == '2'
    assert na12878_df.loc['CYP2D7.e2']['mle_n'] == '3'
