import re
from itertools import product
from typing import Dict, Union, Hashable, Sequence, Any, List

import numpy as np
import pandas as pd
import pyranges as pr
from joblib import Parallel, delayed
from pkg_resources import resource_filename


def subset_dict(d: Dict, key_or_keys: Union[Hashable, Sequence[Hashable]]) -> dict:
    """Subset dict by one or more keys

    Args:
        key_or_keys: All hashable objects are interpreted as
            'one key'. This means that a tuple is *a single key*.
            Several keys should be given as sequence of hashables.
        d: the input dictionary

    Returns:
        A new dictionary (currently *not* a deepcopy), containing
        only the keys specified by key_or_keys
    """
    try:
        if isinstance(key_or_keys, Hashable):
            return {key_or_keys: d[key_or_keys]}
        elif isinstance(key_or_keys, Sequence):
            return {k: d[k] for k in key_or_keys}
        else:
            raise TypeError(f"Subset_dict cant handle input type {type(key_or_keys)}")
    except KeyError:
        raise KeyError(
            "Error while subsetting dict."
            f" Can't subset dict with keys {key_or_keys},"
            " because at least one key is missing"
        )


def dict_to_compact_str(d: Dict[str, Any]) -> str:
    return ",".join(f"{k}={v}" for k, v in d.items())


def get_snakefile_path() -> str:
    return resource_filename("dss_workflow", "dss_workflow.snakefile")


def get_jobscript_path() -> str:
    return resource_filename("dss_workflow", "lsf-jobscript.sh")


def get_statusscript_path() -> str:
    return resource_filename("dss_workflow", "lsf-status.py")


def get_submitscript_path() -> str:
    return resource_filename("dss_workflow", "lsf-submit.py")


def create_dmr_metadata_table(config: Dict[str, Any]) -> pd.DataFrame:
    """Use information in DMR calling workflow config dict to create a metadata table

    The metadata table gives the paths to all files created by the workflow,
    with metadata for each file.

    Currently, this only records DMR calls, not DML calls.

    Parquet output is already included in the metadata table, although it is not yet
    produced by the workflow, see comments in the workflow. To create parquet output,
    create env with r-arrow >= 0.17 and r-stringr, then run something like:

    # + {"language": "R"}
    # library(arrow)
    # library(stringr)
    # dml_test_results <- Sys.glob('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/pairwise-dmr-calls/DSS# /v1_bistro-0.2.0_odcf-alignment/ds1/smooth=False,smooth_span=500/dmls/dmls_hsc_vs_*.rds')
    # for (dml_test_result in dml_test_results) {
    #     parquet_path = str_replace(dml_test_result, '\\.rds', '.parquet')
    #     print(parquet_path)
    #     df <- readRDS(dml_test_result)
    #     write_parquet(df, parquet_path)
    # }
    # -

    Args:
        config: the dictionary passed to the DMR calling snakemake workflow

    Returns:
        A table with the following columns:
        - dmr_bed :: absolute filepath of DMR calls in BED format
        - dmr_coverage_bedgraph :: absolute filepath to coverage bedgraph file,
            if the file exists, otherwise 'NA'
        - parameters: Dict[str, Any] of parameters passed to DSS
        - parameters_str: str representation of the dict
        - group1: name of the first group in the comparison
        - group2: ...
        - comparison name: str label for the comparison
    """
    dmr_by_group1_group2_params = (
        config["output_dir"]
        + "/{dml_params}/dmrs/{dmr_params}/dmrs_{group1}_vs_{group2}.bed"
    )
    dml_by_group1_group2_params_parquet = (
        config["output_dir"] + "/{dml_params}/dmls/dmls_{group1}_vs_{group2}.parquet"
    )

    dmr_coverage_bedgraph_by_group1_group2_params = (
        config["output_dir"]
        + "/{dml_params}/dmrs/{dmr_params}/dmrs_{group1}_vs_{group2}_coverage.bedgraph"
    )

    file_metadata_l: List[Dict[str, Any]] = []
    for group_tuple, parameters in product(config["comparisons"], config["parameters"]):
        # Missing: Assert that the parameters dict contains all and only the expected parameters
        dml_params = subset_dict(parameters, ["smooth", "smooth_span"])
        dml_params_str = dict_to_compact_str(dml_params)
        dmr_params = subset_dict(
            parameters, "pvalue delta minlen minCG merge_dist pct_sign".split()
        )
        dmr_params_str = dict_to_compact_str(dmr_params)
        if config["create_dmr_coverage_bedgraph"]:
            dmr_coverage_bedgraph_fp = dmr_coverage_bedgraph_by_group1_group2_params.format(
                group1=group_tuple[0],
                group2=group_tuple[1],
                dml_params=dml_params_str,
                dmr_params=dmr_params_str,
            )
        else:
            dmr_coverage_bedgraph_fp = "NA"
        file_metadata_l.append(
            {
                "dmr_bed": dmr_by_group1_group2_params.format(
                    group1=group_tuple[0],
                    group2=group_tuple[1],
                    dml_params=dml_params_str,
                    dmr_params=dmr_params_str,
                ),
                "dmr_coverage_bedgraph": dmr_coverage_bedgraph_fp,
                "dml_parquet": dml_by_group1_group2_params_parquet.format(
                    group1=group_tuple[0],
                    group2=group_tuple[1],
                    dml_params=dml_params_str,
                ),
                "parameters": parameters,
                "parameters_str": dict_to_compact_str(parameters),
                "group1": group_tuple[0],
                "group2": group_tuple[1],
                "comparison_name": f"dmrs_{group_tuple[0]}_vs_{group_tuple[1]}",
            }
        )
    metadata_table = pd.DataFrame(file_metadata_l)
    return metadata_table


def get_dml_test_dfs_d(dml_test_files, query_cpgs, n_jobs=1) -> Dict:
    """

    Parameters
    ----------
    dml_test_files
        parquet files of DSS DML result
    query_cpgs
        CpGs for which the pvalues should be retrieved
    n_jobs

    Returns
    -------
    dict
        keys: p_value, q_value, p_value_posterior (for delta > 0.1), mean, se

    """

    dmr_cpg_test_results_d = concat_dml_stats(
        dml_test_files=dml_test_files, cpg_query_df=query_cpgs, n_jobs=n_jobs
    )

    dmr_cpg_test_results_d["p_value_posterior"] = dss_posterior_probability(
        delta_df=dmr_cpg_test_results_d["delta"],
        se_df=dmr_cpg_test_results_d["se"],
        min_delta=0.1,
    )


from scipy.stats import norm


def dss_posterior_probability(delta_df, se_df, min_delta: float):
    # https://github.com/haowulab/DSS/blob/911b1bd091a6a0e39b66284c23e8a972ac716961/R/DML.R
    # - DML calling is either based on the p-values (if no delta is specified)
    # or on the posterior p-values.
    # - No mht correction is applied for DML calling, if I read this correctly
    p1 = norm.cdf(delta_df - min_delta, loc=0, scale=se_df)
    p2 = 1 - norm.cdf(delta_df + min_delta, loc=0, scale=se_df)
    return pd.DataFrame(1 - (p1 + p2), index=delta_df.index, columns=delta_df.columns)


def concat_dml_stats(
    dml_test_files: Dict[str, str], cpg_query_df: pd.DataFrame, n_jobs=1
) -> Dict[str, pd.DataFrame]:
    """Concatenate DML test statistics and p-values on query cpgs

    Collects p-values, q-values, delta, se into individual cpg x pop dataframes

    Parameters
    ----------
    dml_test_files
        Dict: pop_name -> path to dml test result (as parquet file)
        pop_name is used as column name in the resulting dataframes
    cpg_query_df
        Chromosome Start End [optional_col, ...] for all CpGs of interest
        Chromosome must be category, the chromosome dtype is used in the function
    n_jobs
        parallelize over dml_test files

    Returns
    -------
    Dict
        dict with dataframes (cpg, pop), keys: ["p_value", "q_value", "delta", "se"]
    """

    # Implementation notes
    # --------------------
    # - DMLtest results discard CpGs (probably those with no coverage).
    #   - First reindex test results to get same order of CpGs for all files, conforming to cpg_index_df.
    #   - Then extract the query CpGs.
    #   - Finally, create and return DFs for p_value, q_value and delta

    assert cpg_query_df.Chromosome.dtype.name == "category"

    # Extract test results for the query CpGs
    # get a list with one dataframe per pop, query cpgs x [p_value q_value delta]
    dml_stats_per_pop_l = Parallel(n_jobs, backend="multiprocessing", max_nbytes=None)(
        delayed(_get_single_test_stat_df)(dml_test_file, cpg_query_df)
        for name, dml_test_file in dml_test_files.items()
    )

    # get individual dfs (cpg, pop) for p_value, q_value, delta
    # use dml_test_files.keys() as pop column names
    return {
        k: pd.concat(
            [df[k] for df in dml_stats_per_pop_l], axis=1, keys=dml_test_files.keys()
        ).reset_index(drop=True)
        for k in ["p_value", "q_value", "delta", "se"]
    }


def _get_single_test_stat_df(
    dml_test_file: str, query_cpg_df: pd.DataFrame
) -> pd.DataFrame:
    """Retrieve DF of DML test stats for query cpgs from a single DMLtest result

    Used in concat_dml_stats
    """

    # Implementation notes
    # --------------------
    # - DMLtest result discards CpGs (probably those with no coverage).

    chrom_dtype = query_cpg_df.Chromosome.dtype
    query_cpg_gr = pr.PyRanges(query_cpg_df)

    # DF: chr pos mu1 mu2 diff diff.se stat phi1 phi2 pval fdr
    dml_test_df = (
        pd.read_parquet(
            dml_test_file, columns=["chr", "pos", "diff", "pval", "fdr", "diff.se"]
        ).rename(
            columns={
                "chr": "Chromosome",
                "pos": "Start",
                "diff": "delta",
                "pval": "p_value",
                "fdr": "q_value",
                "diff.se": "se",
            }
        )
        # 1-based start coords
        # no End column present
        .assign(
            Start=lambda df: df.Start - 1,
            End=lambda df: df.Start + 2,
            Chromosome=lambda df: df.Chromosome.astype(chrom_dtype),
        )[["Chromosome", "Start", "End", "p_value", "delta", "q_value", "se"]]
    )
    dml_test_gr = pr.PyRanges(dml_test_df)

    # join query cpg GR and dml test GR
    # sort result according to correct chrom dtype
    reindexed_dml_test_gr = query_cpg_gr.join(dml_test_gr, how="left")
    reindexed_dml_test_df = (
        reindexed_dml_test_gr.df.assign(
            Chromosome=lambda df: df.Chromosome.astype(chrom_dtype)
        )
        .sort_values(["Chromosome", "Start", "End"])
        .reset_index(drop=True)
    )

    # assert that joined result matches query cpgs
    pd.testing.assert_frame_equal(
        query_cpg_df[["Chromosome", "Start", "End"]],
        reindexed_dml_test_df[["Chromosome", "Start", "End"]],
        check_names=False,
        check_dtype=False,
    )

    # currently, join uses -1 instead of NA for unmatched rows
    # find them, replace with NA
    matched_rows = ~reindexed_dml_test_df.Start_b.eq(-1)
    # noinspection PyUnresolvedReferences
    assert (
        reindexed_dml_test_df.loc[matched_rows, ["Start", "End"]].to_numpy()
        == reindexed_dml_test_df.loc[matched_rows, ["Start_b", "End_b"]].to_numpy()
    ).all()
    reindexed_dml_test_df.loc[
        ~matched_rows, ["delta", "se", "p_value", "q_value"]
    ] = np.nan

    # only return stat columns of interest
    dmr_cpg_dml_test_results = reindexed_dml_test_df[
        ["delta", "se", "p_value", "q_value"]
    ]

    return dmr_cpg_dml_test_results

