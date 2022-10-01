import json
import pathlib
import sys


# Set below at runtime
_MEDIA_DEFS = None


def _discover_package_media_files():
    file_fp = pathlib.Path(__file__).absolute()
    media_dir = file_fp.parent / 'data' / 'media_definitions'
    return [fp for fp in media_dir.glob('*json')]


def _create_media_map():
    media_fps = _discover_package_media_files()
    return {fp.stem.replace('_media', ''): fp for fp in media_fps}


def available():
    return list(_MEDIA_DEFS.keys())


def get(name):
    if (media_fp := _MEDIA_DEFS.get(name)):
        with media_fp.open('r') as fh:
            return json.load(fh)
    else:
        defs_str = '\n\t'.join(available())
        msg = f'Could not find {name} definition, available definitions are:\n\t{defs_str}'
        print(msg, file=sys.stderr)
        sys.exit(1)


_MEDIA_DEFS = _create_media_map()
