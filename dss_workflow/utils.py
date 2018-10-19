from typing import Dict, Union, Hashable, Sequence, Any

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
