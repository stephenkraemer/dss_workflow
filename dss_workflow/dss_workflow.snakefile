"""Pairwise DMR calling workflow

Example call:

    snakemake \
    --snakefile /home/kraemers/projects/dmr_calling/dmr-calling.snakefile \
    --configfile /path/to/config/file \
    --jobs 1000 \
    --latency-wait 180 \
    --jobscript /home/kraemers/projects/dmr_calling/local/jobscript.sh \
    --cluster "bsub -R rusage[mem={{params.avg_mem}}] -M {{params.max_mem}} -n {{threads}} -J {{params.name}} -W {{params.walltime}} -o /home/kraemers/temp/logs/" \
    --conda-prefix /icgc/dkfzlsdf/analysis/hs_ontogeny/envs/snakemake-managed \
    --use-conda \

"""

import pandas as pd
from dss_workflow.utils import subset_dict, dict_to_compact_str
from itertools import product

# Command line args handling
# --------------------------

metadata_table = pd.read_csv(config['metadata_table'], sep='\t', header=0)
required_cols = (subset_dict(config, ['group_column', 'sample_id_column', 'bed_by_chrom_column']).values())
assert all([s in metadata_table.columns for s in required_cols])
metadata_table = metadata_table.set_index([config['group_column'], config['sample_id_column']], drop=False)
print(f'Reading metadata table, with head:\n{metadata_table.head()}')
assert metadata_table[config['bed_by_chrom_column']].str.contains('{chrom}').all()

assert len(config['comparisons'][0]) == 2
print(f'Writing to {config["output_dir"]}')

# Path definitions
# ----------------
dml_by_group1_group2_params_chrom = config['output_dir'] + '/{dml_params}/dmls/per_chrom/dmls_{group1}_vs_{group2}_chrom-{chrom}.rds'
dml_by_group1_group2_params = config['output_dir'] + '/{dml_params}/dmls/dmls_{group1}_vs_{group2}.rds'
dmr_by_group1_group2_params = config['output_dir'] + '/{dml_params}/dmrs/{dmr_params}/dmrs_{group1}_vs_{group2}.bed'
dmr_coverage_bedgraph_by_group1_group2_params = config['output_dir'] + '/{dml_params}/dmrs/{dmr_params}/dmrs_{group1}_vs_{group2}_coverage.bedgraph'

# Construct targets
# ------------------

dmr_by_group1_group2_params = config['output_dir'] + '/{dml_params}/dmrs/{dmr_params}/dmrs_{group1}_vs_{group2}.bed'
dml_params_name_to_params_dict = {}
dmr_params_name_to_params_dict = {}
dmr_calls = []
dmr_coverage_bedgraphs = []
for group_tuple, parameters in product(config['comparisons'], config['parameters']):
    dml_params = subset_dict(parameters, ['smooth', 'smooth_span'])
    dml_params_str = dict_to_compact_str(dml_params)
    if dml_params_str not in dml_params_name_to_params_dict:
        dml_params_name_to_params_dict[dml_params_str] = dml_params
    dmr_params = subset_dict(parameters, 'pvalue delta minlen minCG merge_dist pct_sign'.split())
    dmr_params_str = dict_to_compact_str(dmr_params)
    if dmr_params_str not in dmr_params_name_to_params_dict:
        dmr_params_name_to_params_dict[dmr_params_str] = dmr_params
    dmr_calls.append(dmr_by_group1_group2_params.format(
            group1=group_tuple[0], group2=group_tuple[1],
            dml_params=dml_params_str, dmr_params=dmr_params_str))
    dmr_coverage_bedgraphs.append(dmr_coverage_bedgraph_by_group1_group2_params.format(
        group1=group_tuple[0], group2=group_tuple[1],
        dml_params=dml_params_str, dmr_params=dmr_params_str))
print(dmr_calls)

wildcard_constraints:
    chrom = "\d{1,2}"

targets = [dmr_calls]
if config['create_dmr_coverage_bedgraph']:
    targets.append(dmr_coverage_bedgraphs)

rule all:
    input:
        targets


rule call_dmls:
    input:
        group1_calls = lambda wildcards: metadata_table.loc[wildcards.group1, config['bed_by_chrom_column']].str.replace('{chrom}', wildcards.chrom),
        group2_calls = lambda wildcards: metadata_table.loc[wildcards.group2, config['bed_by_chrom_column']].str.replace('{chrom}', wildcards.chrom),
    params:
        group1_sample_ids = lambda wildcards: metadata_table.loc[wildcards.group1, 'sample_id'].tolist(),
        group2_sample_ids = lambda wildcards: metadata_table.loc[wildcards.group2, 'sample_id'].tolist(),
        dml_test_parameters = lambda wildcards: dml_params_name_to_params_dict[wildcards.dml_params],
        walltime = "03:00",
        avg_mem = 4000,
        max_mem = 6000,
        name = "{group1}_vs_{group2}_chr{chrom}",
    threads: 1
    conda: 'dss_workflow.yml'
    output:
        dml_by_group1_group2_params_chrom,
    script:
        "call_dmls.R"


rule collect_dml_calls:
    input:
        lambda wildcards: expand(dml_by_group1_group2_params_chrom,
                                 group1=wildcards.group1,
                                 group2=wildcards.group2,
                                 dml_params=wildcards.dml_params,
                                 chrom=config['chromosomes'])
    params:
        walltime = "00:15",
        avg_mem = 6000,
        max_mem = 10000,
        name = "merge-dml_{group1}_vs_{group2}",
    threads: 1
    conda: 'dss_workflow.yml'
    output:
        dml_by_group1_group2_params,
    script:
        "concat_dml_chrom_dfs.R"


rule call_dmrs:
    input:
        dml_by_group1_group2_params,
    params:
        walltime = "00:15",
        avg_mem = 4000,
        max_mem = 8000,
        name = "dmr-calling_{group1}_vs_{group2}",
        dmr_test_parameters = lambda wildcards: dmr_params_name_to_params_dict[wildcards.dmr_params],
        chromosomes = config['chromosomes'],
    threads: 1
    conda: 'dss_workflow.yml'
    output:
        dmr_by_group1_group2_params,
    script:
        "call_dmrs.R"

rule coverage_bedgraph:
    input:
        dmr_bed=dmr_by_group1_group2_params,
        index=config['cpg_index_file'],
    output:
        dmr_coverage_bedgraph_by_group1_group2_params,
    params:
        walltime = "00:10",
        avg_mem = 3000,
        max_mem = 7000,
        name = "dmr-coverage_{group1}_vs_{group2}",
    threads: 1
    conda: 'dss_workflow.yml'
    shell:
        """
        echo '#chr\tstart\tend\tis_in_dmr' > {output}
        bedtools coverage -sorted -a {input.index} -b {input.dmr_bed} | cut -f 1,2,3,4 >> {output}
        """
