import argparse
import pathlib
import sys


from . import __program_name__
from . import __version__


class ArgumentParserCustomHelp(argparse.ArgumentParser):

    def parse_args(self, args=None, namespace=None):
        args, argv = self.parse_known_args(args, namespace)
        if argv:
            print(help_text(args.command), file=sys.stderr)
            plurality = 'argument' if len(argv) == 1 else 'arguments'
            msg = f'\n{__program_name__}: error: unrecognized {plurality}: {" ".join(argv)}'
            print(msg, file=sys.stderr)
            sys.exit(1)
        return args


    def _check_value(self, action, value):
        if action.choices is not None and value not in action.choices:
            if action.dest == 'command':
                print(help_text(None), file=sys.stderr)
                print(f'\n{__program_name__}: error: unrecognized stage: {value}', file=sys.stderr)
                sys.exit(1)
            else:
                args = {'value': value,
                        'choices': ', '.join(map(repr, action.choices))}
                msg = 'invalid choice: %(value)r (choose from %(choices)s)'
                raise argparse.ArgumentError(action, msg % args)


def parse():
    parser = ArgumentParserCustomHelp(add_help=False)
    subparsers = parser.add_subparsers(dest='command')

    parser.add_argument('--assembly_fps', type=pathlib.Path, nargs='+')
    parser.add_argument('--ref_gbk_fp', type=pathlib.Path)
    parser.add_argument('--ref_model_fp', type=pathlib.Path)
    parser.add_argument('--output_fp', type=pathlib.Path)

    parser.add_argument('-v', '--version', action='store_true')
    parser.add_argument('-h', '--help', action='store_true')

    parser_identorth = subparsers.add_parser('draft_model', add_help=False)
    parser_identorth.add_argument('--assembly_fp', type=pathlib.Path)
    parser_identorth.add_argument('--ref_gbk_fp', type=pathlib.Path)
    parser_identorth.add_argument('--ref_model_fp', type=pathlib.Path)
    parser_identorth.add_argument('--output_fp', type=pathlib.Path)
    parser_identorth.add_argument('-h', '--help', action='store_true')

    args = parser.parse_args()
    check_arguments(args)
    return args


def check_arguments(args):
    if args.help:
        print(help_text(args.command), file=sys.stdout)
        sys.exit(0)
    if args.version:
        print(f'{__program_name__} {__version__}', file=sys.stdout)
        sys.exit(0)
    # Check we have required arguments, this is purposely decouped from argparse
    required_args = {
        'base': ('assembly_fps', 'ref_gbk_fp', 'ref_model_fp', 'output_fp'),
        'draft_model': ('assembly_fp', 'ref_gbk_fp', 'ref_model_fp', 'output_fp')
    }
    command = 'base' if not args.command else args.command
    assert command in required_args
    missing_args = list()
    for required_arg in required_args[command]:
        if args.__dict__[required_arg]:
            continue
        missing_args.append(required_arg)
    if missing_args:
        missing_args_str = ', '.join(f'--{arg}' for arg in missing_args)
        plurality = 'argument is' if len(missing_args) == 1 else 'arguments are'
        msg = f'{__program_name__}: error: the following {plurality} required {missing_args_str}'
        print(help_text(args.command), file=sys.stderr)
        print(msg, file=sys.stderr)
        sys.exit(1)
    # Check all input file objects exist
    for arg, value in args.__dict__.items():
        if not value or arg == 'output_fp':
            continue
        elif isinstance(value, pathlib.Path):
            if not value.exists():
                print(f'{__program_name__}: error: input {value} does not exist', file=sys.stderr)
                sys.exit(1)
        elif isinstance(value, list):
            for v in value:
                if not isinstance(value, pathlib.Path):
                    continue
                if not v.exists():
                    print(f'{__program_name__}: error: input {v} does not exist', file=sys.stderr)
                    sys.exit(1)


def help_text(command):
    info_text = f'\nProgram: {__program_name__}\nVersion: {__version__}\n'
    if not command:
        help_text = (f'Pipeline usage: {__program_name__} [options]\n'
                     f'Single stage usage: {__program_name__} <stage> [options]\n\n'
                      'Pipeline options:\n'
                      '  --assembly_fps FILES        Isolate genbank filepaths\n'
                      '  --ref_gbk_fp FILE           Reference gebank filepath\n'
                      '  --ref_model_fp FILE         Reference gebank filepath (JSON)\n'
                      '  --output_fp DIR             Output directory\n\n'
                      'Other options:\n'
                      '  --version                   Print program name and version, and exit\n'
                      '  --help                      Print this message and exit)\n\n'
                      'Single stage subcommands:\n'
                      '  assembly_qc                 QC for input assembly\n'
                      '  annotation                  Annotate assembly ORFs\n'
                      '  draft_model                 Create a draft model\n'
                      '  model_fba                   Simulate growth on media with FBA\n\n'
                     f'For more information about single stage subcommands run: {__program_name__} <stage> --help\n')
    elif command == 'draft_model':
        help_text = (f'Usage: {__program_name__} {command} [options]\n'
                      'Options:\n'
                      '  --assembly_fp FILES         Isolate genbank filepath\n'
                      '  --ref_gbk_fp FILE           Reference gebank filepath\n'
                      '  --ref_model_fp FILE         Reference gebank filepath (JSON)\n'
                      '  --output_fp DIR             Output directory\n')
    else:
        assert False
    return f'{info_text}\n{help_text}'
