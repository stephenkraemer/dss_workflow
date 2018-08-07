"""Pairwise DMR calling workflow

Args:

    [workflow parameters]
      metadata_table: tsv with columns: path_with_chrom_field, sample_id, group
      comparisons: comma-separated comparison tasks, e.g. group_name1-VS-group_name2,group_name1-VS-group_name3
      output_dir: directory where the computed results are stored
      chromosomes: comma_separated chromosomes, the order defines the chromosome sorting order of the DML and DMR dataframes. E.g. chromosomes=1,10,12,...,2,3,...,9 for alphabetically sorted mouse autosomes
    [dmr calling parameters]
      dmr_pvalue: argument for DSS::callDMR
      dmr_delta: argument for DSS::callDMR
      dmr_minlen=50 \
      dmr_minCG=3 \
      dmr_distance_merge_bp=50 \
      dmr_pct_significant=0.5 \

Notes:

  - take care to specify the chromosomes in the desired sorting order (alphabetical or numerical). The specified order will be used throughout the workflow
  - the metadata table is allowed to have additional columns.
  - the DMLtest function is called with hardcoded parameters (smoothing=True, window size=500bp, no equal dispersion for both groups). In the future it may be necessary to adjust this hardcoding to changed function signatures, or through expose these options through the cli.

Example call:

    snakemake \
      --snakefile /home/kraemers/projects/dmr_calling/dmr-calling.snakefile \
      --config \
        metadata_table=/home/kraemers/temp/metadata_table.csv \
        comparisons='hsc-VS-mpp1,hsc-VS-mpp2' \
        output_dir=/home/kraemers/temp/dmr_calling_test \
        chromosomes=1,10,11,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9 \
        dmr_pvalue=0.05 \
        dmr_delta=0.1 \
        dmr_minlen=50 \
        dmr_minCG=3 \
        dmr_distance_merge_bp=50 \
        dmr_pct_significant=0.5 \
      --jobs 100 \
      --latency-wait 180 \
      --jobscript /home/kraemers/projects/dmr_calling/local/jobscript.sh \
      --cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
      --dryrun \

      # other useful options
      --forcerun call_dmls \

"""
import pandas as pd

# Command line args handling
# --------------------------

metadata_table = pd.read_csv(config['metadata_table'], sep='\t', header=0)
required_cols = ['path_with_chrom_field', 'group', 'sample_id']
assert all([s in metadata_table.columns for s in required_cols])
metadata_table = metadata_table.set_index(['group', 'sample_id'], drop=False)
print(f'Reading metadata table, with head:\n{metadata_table.head()}')
assert metadata_table['path_with_chrom_field'].str.contains('{chrom}').all()

comparisons = [x.split('-VS-') for x in config['comparisons'].split(',')]
assert len(comparisons) >= 1
assert len(comparisons[0]) == 2
print(f'Running the following comparisons:\n{comparisons}')
print(f'Writing to {config["output_dir"]}')

# Path definitions
# ----------------
dml_by_group1_group2_chrom = config['output_dir'] + '/dmls/per_chrom/dmls_{group1}_vs_{group2}_chrom-{chrom}.rds'
dml_by_group1_group2 = config['output_dir'] + '/dmls/dmls_{group1}_vs_{group2}.rds'
dmr_by_group1_group2 = config['output_dir'] + '/dmrs/dmrs_{group1}_vs_{group2}.rds'


wildcard_constraints:
    chrom = "\d{1,2}"


rule all:
    input:
        # expand([dml_by_group1_group2_chrom.format(chrom='{chrom}', group1=group1, group2=group2) for group1, group2 in comparisons],
               # chrom=config['chromosomes']),
        [dmr_by_group1_group2.format(group1=group1, group2=group2) for group1, group2 in comparisons],


# "worker_scripts/create_bsseq.R"

rule call_dmls:
    input:
        group1_calls = lambda wildcards: metadata_table.loc[wildcards.group1, 'path_with_chrom_field'].str.replace('{chrom}', wildcards.chrom),
        group2_calls = lambda wildcards: metadata_table.loc[wildcards.group2, 'path_with_chrom_field'].str.replace('{chrom}', wildcards.chrom),
    params:
        group1_sample_ids = lambda wildcards: metadata_table.loc[wildcards.group1, 'sample_id'].tolist(),
        group2_sample_ids = lambda wildcards: metadata_table.loc[wildcards.group2, 'sample_id'].tolist(),
        walltime = "01:00",
        avg_mem = 2500,
        max_mem = 4000,
        name = "{group1}_vs_{group2}_chr{chrom}",
    threads: 1
    output:
        dml_by_group1_group2_chrom,
    script:
        "call_dmls.R"

rule collect_dml_calls:
    input:
        lambda wildcards: expand(dml_by_group1_group2_chrom,
                                 group1=wildcards.group1,
                                 group2=wildcards.group2,
                                 chrom=config['chromosomes'])
    params:
        walltime = "00:15",
        avg_mem = 6000,
        max_mem = 10000,
        name = "merge-dml_{group1}_vs_{group2}",
    threads: 1
    output:
        dml_by_group1_group2,
    script:
        "concat_dml_chrom_dfs.R"

rule call_dmrs:
    input:
        dml_by_group1_group2,
    params:
        walltime = "00:15",
        avg_mem = 4000,
        max_mem = 8000,
        name = "dmr-calling_{group1}_vs_{group2}",
        dmr_pvalue = config['dmr_pvalue'],
        dmr_delta = config['dmr_delta'],
        dmr_minlen = config['dmr_minlen'],
        dmr_minCG = config['dmr_minCG'],
        dmr_distance_merge_bp = config['dmr_distance_merge_bp'],
        dmr_pct_significant = config['dmr_pct_significant'],
    threads: 1
    output:
        dmr_by_group1_group2,
    script:
        "call_dmrs.R"

# Test code
# ====================
# import glob
# import pandas as pd

# glob_pattern = '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/*/meth/meth_calls/mcalls_*_CG_8_strands-merged.bed.gz'
# regex_pattern =  '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/(?P<group>.*)_(?P<rep>\d+)/meth/meth_calls/mcalls_.+_CG_{chrom}_strands-merged.bed.gz'
# metadata_table = pd.DataFrame(dict(path_with_chrom_field = glob.glob(glob_pattern)))
# metadata_table['path_with_chrom_field'] = metadata_table['path_with_chrom_field'].str.replace('_8_', '_{chrom}_')
# anno = metadata_table['path_with_chrom_field'].str.extract(regex_pattern)
# anno = anno.assign(sample_id = lambda df: df['group'] + '_' + df['rep']).drop('rep', axis=1)
# anno.head()
# metadata_table.head()
# metadata_table_cat = pd.concat([metadata_table, anno], axis=1)
# metadata_table_cat.head()
# metadata_table_cat['path'] = 'should not be used'
# metadata_table_cat['other_col'] = 'should not be a problem'
# print(metadata_table_cat.head())
# metadata_table_cat.to_csv('/home/kraemers/temp/metadata_table.csv', sep='\t', index=False)

# config = {}
# config['metadata_table'] = '/home/kraemers/temp/metadata_table.csv'
# config['comparisons'] = 'hsc-VS-mpp1,hsc-VS-mpp2'
# config['output_dir'] = '/home/kraemers/temp/dmr_calling_test'
# config['chromosomes'] = sorted([str(i) for i in range(1, 20)])
# config['dmr_pvalue'] = 0.05
# config['dmr_delta'] = 0.1
