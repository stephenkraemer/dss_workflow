from itertools import product
from typing import Dict, Union, Hashable, Sequence, Any, List

import pandas as pd
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
            raise TypeError(f'Subset_dict cant handle input type {type(key_or_keys)}')
    except KeyError:
        raise KeyError("Error while subsetting dict."
                       f" Can't subset dict with keys {key_or_keys},"
                       " because at least one key is missing")


def dict_to_compact_str(d: Dict[str, Any]) -> str:
    return ','.join(f'{k}={v}' for k, v in d.items())


def get_snakefile_path() -> str:
    return resource_filename('dss_workflow', 'dss_workflow.snakefile')


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