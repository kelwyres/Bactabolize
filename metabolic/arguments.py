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
    parser.add_argument('-v', '--version', action='store_true')
    parser.add_argument('-h', '--help', action='store_true')

    subparsers = parser.add_subparsers(dest='command')

    parser_draft = subparsers.add_parser('draft_model', add_help=False)
    parser_draft.add_argument('--assembly_fp', type=pathlib.Path)
    parser_draft.add_argument('--ref_genbank_fp', type=pathlib.Path)
    parser_draft.add_argument('--ref_model_fp', type=pathlib.Path)
    parser_draft.add_argument('--output_fp', type=pathlib.Path)
    parser_draft.add_argument('-h', '--help', action='store_true')

    parser_draft.add_argument('--prodigal_model_fp', type=pathlib.Path)
    parser_draft.add_argument('--no_reannotation', action='store_true')

    parser_patch = subparsers.add_parser('patch_model', add_help=False)
    parser_patch.add_argument('--draft_model_fp', type=pathlib.Path)
    parser_patch.add_argument('--ref_model_fp', type=pathlib.Path)
    parser_patch.add_argument('--patch_fp', type=pathlib.Path)
    parser_patch.add_argument('--output_fp', type=pathlib.Path)
    parser_patch.add_argument('-h', '--help', action='store_true')

    parser_fba = subparsers.add_parser('fba', add_help=False)
    parser_fba.add_argument('--model_fp', type=pathlib.Path)
    parser_fba.add_argument('--fba_spec_fp', type=pathlib.Path)
    parser_fba.add_argument('--output_fp', type=pathlib.Path)
    parser_fba.add_argument('-h', '--help', action='store_true')

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
        'draft_model': ('assembly_fp', 'ref_genbank_fp', 'ref_model_fp', 'output_fp'),
        'patch_model': ('draft_model_fp', 'ref_model_fp', 'patch_fp', 'output_fp'),
        'fba': ('model_fp', 'fba_spec_fp', 'output_fp')
    }
    if not args.command:
        msg = f'{__program_name__}: error: you must provide a command to execute'
        print(help_text(args.command), file=sys.stderr)
        print(msg, file=sys.stderr)
        sys.exit(1)
    else:
        assert args.command in required_args
    missing_args = list()
    for required_arg in required_args[args.command]:
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


def help_text(command):
    info_text = f'\nProgram: {__program_name__}\nVersion: {__version__}\n'
    if not command:
        help_text = (f'Usage: {__program_name__} <command> [options]\n\n'
                      'Commands:\n'
                      '  draft_model                 Create a draft model\n'
                      '  patch_model                 Patch a draft model\n'
                      '  fba                         Simulate growth on media with FBA\n\n'
                     f'For more information about commands, run: {__program_name__} <command> --help\n')
    elif command == 'draft_model':
        help_text = (f'Usage: {__program_name__} {command} [options]\n'
                      'Options:\n'
                      '  --assembly_fp FILE          Isolate assembly filepath (GenBank format)\n'
                      '  --ref_genbank_fp FILE       Reference genbank filepath\n'
                      '  --ref_model_fp FILE         Reference model filepath (JSON)\n'
                      '  --output_fp FILE            Output filepath\n'
                      '\nOther:\n'
                      '  --prodigal_model_fp FILE    Prodigal model to use\n'
                      '  --no_reannotation           Do not reannotate genbank file\n')
    elif command == 'patch_model':
        help_text = (f'Usage: {__program_name__} {command} [options]\n'
                      'Options:\n'
                      '  --draft_model_fp FILE       Isolate model filepath\n'
                      '  --ref_model_fp FILE         Reference model filepath\n'
                      '  --patch_fp FILE             Patch file (JSON)\n'
                      '  --output_fp FILE            Output filepath\n')
    elif command == 'fba':
        help_text = (f'Usage: {__program_name__} {command} [options]\n'
                      'Options:\n'
                      '  --model_fp FILE       Isolate model filepath\n'
                      '  --fba_spec_fp FILE          FBA spec filepath (JSON format)\n'
                      '  --output_fp FILE            Output filepath\n')
    else:
        assert False
    return f'{info_text}\n{help_text}'
