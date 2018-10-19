from pkg_resources import resource_filename
from itertools import product

import pandas as pd
from typing import Dict, List, Any

from dss_workflow.utils import subset_dict, dict_to_compact_str


def get_snakefile_path() -> str:
    return resource_filename('dss_workflow', 'dmr-calling.snakefile')


def create_dmr_metadata_table(config: Dict[str, Any]) -> pd.DataFrame:
    """Use information in DMR calling workflow config dict to create a metadata table

    The metadata table gives the paths to all files created by the workflow,
    with metadata for each file.

    Currently, this only records DMR calls, not DML calls.

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
            config['output_dir']
            + '/{dml_params}/dmrs/{dmr_params}/dmrs_{group1}_vs_{group2}.bed')
    dmr_coverage_bedgraph_by_group1_group2_params = (
            config['output_dir']
            + '/{dml_params}/dmrs/{dmr_params}/dmrs_{group1}_vs_{group2}_coverage.bedgraph')

    file_metadata_l: List[Dict[str, Any]] = []
    for group_tuple, parameters in product(config['comparisons'],
                                           config['parameters']):
        # Missing: Assert that the parameters dict contains all and only the expected parameters
        dml_params = subset_dict(parameters, ['smooth', 'smooth_span'])
        dml_params_str = dict_to_compact_str(dml_params)
        dmr_params = subset_dict(parameters, 'pvalue delta minlen minCG merge_dist pct_sign'.split())
        dmr_params_str = dict_to_compact_str(dmr_params)
        if config['create_dmr_coverage_bedgraph']:
            dmr_coverage_bedgraph_fp = dmr_coverage_bedgraph_by_group1_group2_params.format(
                    group1=group_tuple[0], group2=group_tuple[1],
                    dml_params=dml_params_str, dmr_params=dmr_params_str)
        else:
            dmr_coverage_bedgraph_fp = 'NA'
        file_metadata_l.append({
            'dmr_bed': dmr_by_group1_group2_params.format(
                    group1=group_tuple[0], group2=group_tuple[1],
                    dml_params=dml_params_str, dmr_params=dmr_params_str),
            'dmr_coverage_bedgraph': dmr_coverage_bedgraph_fp,
            'parameters': parameters,
            'parameters_str': dict_to_compact_str(parameters),
            'group1': group_tuple[0],
            'group2': group_tuple[1],
            'comparison_name': f'dmrs_{group_tuple[0]}_vs_{group_tuple[1]}',
        })
    metadata_table = pd.DataFrame(file_metadata_l)
    return metadata_table
