import argparse
import pathlib
import shutil


from . import __version__


class WideHelpFormatter(argparse.HelpFormatter):

    def __init__(self, *args, **kwargs):
        terminal_width = shutil.get_terminal_size().columns
        help_width = min(terminal_width, 140)
        super().__init__(*args, **kwargs, max_help_position=help_width, width=help_width)


def get_args():
    parser = argparse.ArgumentParser(formatter_class=WideHelpFormatter, add_help=False)
    parser.add_argument('--isolate_fp', required=True, type=pathlib.Path,
            help='Isolate genbank filepath')
    parser.add_argument('--ref_gbk_fp', required=True, type=pathlib.Path,
            help='Reference genbank filepath')
    parser.add_argument('--ref_model_fp', required=True, type=pathlib.Path,
            help='Reference model filepath')

    parser.add_argument('--output_fp', required=True, type=pathlib.Path,
            help='Output filepath')

    parser.add_argument('-v', '--version', action='version',
            version='%(prog)s {}'.format(__version__), help='Show version number and exit')
    parser.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()
    if not args.isolate_fp.exists():
        parser.error(f'Input file {args.isolate_fp} does not exist')
    if not args.ref_gbk_fp.exists():
        parser.error(f'Input file {args.ref_gbk_fp} does not exist')
    if not args.ref_model_fp.exists():
        parser.error(f'Input file {args.ref_model_fp} does not exist')
    return args
