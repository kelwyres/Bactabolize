# pylint: disable=redefined-outer-name


import json
import pathlib
import sys


# Set below at runtime
_DATA = {
    'fba_specs': None,
    'media_definitions': None,
}


def _create_map(datatype):
    # pylint: disable=no-else-return
    assert datatype in _DATA

    # Discover files
    file_fp = pathlib.Path(__file__).absolute()
    datatype_dir = file_fp.parent / 'data' / datatype
    datatype_fps = list(datatype_dir.glob('*json'))

    # Process files appropriately
    if datatype == 'fba_specs':
        return {fp.stem.replace('_spec', ''): fp for fp in datatype_fps}
    elif datatype == 'media_definitions':
        return {fp.stem.replace('_media', ''): fp for fp in datatype_fps}
    else:
        assert False


def available(datatype):
    assert datatype in _DATA
    return list(_DATA[datatype].keys())


def get_fp(datatype, name):
    # pylint: disable=no-else-return
    assert datatype in _DATA

    if datatype_fp := _DATA[datatype].get(name):
        return datatype_fp
    else:
        defs_str = '\n\t'.join(available(datatype))
        msg = f'Could not find {name} {datatype}, available definitions are:\n\t{defs_str}'
        print(msg, file=sys.stderr)
        sys.exit(1)


def get_data(datatype, name):
    with get_fp(datatype, name).open('r') as fh:
        return json.load(fh)


# Get available files for each defined datatype
for datatype in _DATA.copy():
    _DATA[datatype] = _create_map(datatype)
