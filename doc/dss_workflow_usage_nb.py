# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
# ---

# # Implementation notes

# DML parquet output: currently not possible from within the workflow
# (see comments on rule). Files are still present in metadata table.
# Use utils.dml_test_results_to_parquet to create parquet files,
# from env with r-arrow (not the workflow env, which is anyway automatically
# created and handled by the workflow)

# # Usage

# ## Requirements

# - conda env with
#   - snakemake
#   - pandas

# ## Clone repo

# git clone https://github.com/stephenkraemer/dss_workflow

# ## Install package with

# pip -e /path/to/dss_workflow

# # in a python script

import dss_workflow
import snakemake

# specify the conda prefix dir

snakemake_conda_prefix = "/path/to/envs"

# Create a config dict, for example

dmr_calling_config = {

    # Inputs
    # ======
    # you need a metadata table specifying where your input files are
    # look at the example here
    "metadata_table":      "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/metadata/hierarchy/meth-calls_v1_bistro-0.2.0_odcf-alignment_datafreeze-1.tsv",
    # Required are columns specifying the group (eg 'Tumor'), the sample_id (eg 'Patient1') and the path to methylation calls in BED format, split by chromosome
    # BED files are BED6 + beta_value n_meth n_total, with a header line
    # See here for an example:
    # /icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_pmcalls_b-cells_1_CG_chrom-merged_strands-merged.bed.gz
    "group_column":        "subject",
    "sample_id_column":    "sample_id",
    "bed_by_chrom_column": "bed_by_chrom_path",
    # Simple BED3 CpG index file, no header
    # See example - these are just the coordinates used in all the methylation
    # calling BEDs
    "cpg_index_file": "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/methylation_calling/indices/GRCm38mm10-phix-lambda/GRCm38mm10-phix-lambda_CG_strands-merged.bed.gz",

    # Outputs
    # =======
    # all results are placed within here
    "output_dir": "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/pairwise-dmr-calls/DSS/v1_bistro-0.2.0_odcf-alignment/ds1",
    # Optionally create coverage bedgraph - hasn't been run in a while
    # not sure if still functional
    "create_dmr_coverage_bedgraph": False,

    # Tasks
    # =====
    # will only run on specified chromosomes
    "chromosomes": [
        "1",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
    ],
    # The pairwise comparisons you want, (group1, group2)
    "comparisons": [
        ("hsc", "mpp1"),
        ("hsc", "mpp5"),
        # ...
    ],
    # The parameter sets you want to screen
    "parameters": [
        {
            "pvalue": 0.01,
            "delta": 0.1,
            "minlen": 50,
            "minCG": 2,
            "merge_dist": 50,
            "pct_sign": 0.5,
            "smooth": False,
            "smooth_span": 500,
        },
        {
            "pvalue": 0.001,
            "delta": 0.1,
            "minlen": 50,
            "minCG": 2,
            "merge_dist": 50,
            "pct_sign": 0.5,
            "smooth": False,
            "smooth_span": 500,
        },
    ],
}


# Run snakemake; this should work directly on the DKFZ LSF cluster using the supplied clutster interaction scripts

# Note
# - that we use the python package to find all resource files
# - that we can allow multiple restarts, at each restart, the requested resources will be enlarged

snakemake.snakemake(
    snakefile=dss_workflow.get_snakefile_path(),
    latency_wait=60,
    config=dmr_calling_config,
    nodes=5000,
    max_jobs_per_second=10,
    jobscript=dss_workflow.get_jobscript_path(),
    cluster=dss_workflow.get_submitscript_path(),
    cluster_status=dss_workflow.get_statusscript_path(),
    dryrun=True,
    conda_prefix=snakemake_conda_prefix,
    use_conda=True,
    restart_times=4,
    keepgoing=True,
    # may be necessary if run was aborted
    # force_incomplete=True,
)

# Optionally, create a record of all files that were just produced as metadata table

results_metadata_table = dss_workflow.create_dmr_metadata_table(
    dmr_calling_config
)

# Note that parquet files are currently not automatically produced, see impl. notes above. Parquet files are handy for collaboration with R users.
