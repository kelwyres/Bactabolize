import argparse
import pathlib
import sys


from . import __program_name__
from . import __version__
from . import media_definitions


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
                args = {'value': value, 'choices': ', '.join(map(repr, action.choices))}
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
    parser_draft.add_argument('--ref_proteins_fp', type=pathlib.Path)
    parser_draft.add_argument('--ref_genes_fp', type=pathlib.Path)
    parser_draft.add_argument('--ref_model_fp', type=pathlib.Path)
    parser_draft.add_argument('--min_coverage', type=float, default=25)
    parser_draft.add_argument('--min_pident', type=float, default=80)
    parser_draft.add_argument('--min_ppos', type=float)
    parser_draft.add_argument('--media_type', type=str, default='m9', choices=media_definitions.available())
    parser_draft.add_argument('--atmosphere_type', type=str, choices=['aerobic', 'anaerobic'])
    parser_draft.add_argument('--biomass_reaction_id', type=str, default='BIOMASS_')
    parser_draft.add_argument('--output_fp', type=pathlib.Path)
    parser_draft.add_argument('--memote_report_fp', type=pathlib.Path)
    parser_draft.add_argument('-h', '--help', action='store_true')

    parser_draft.add_argument('--no_reannotation', action='store_true')

    parser_patch = subparsers.add_parser('patch_model', add_help=False)
    parser_patch.add_argument('--draft_model_fp', type=pathlib.Path)
    parser_patch.add_argument('--ref_model_fp', type=pathlib.Path)
    parser_patch.add_argument('--patch_fp', type=pathlib.Path)
    parser_patch.add_argument('--media_type', type=str, default='m9', choices=media_definitions.available())
    parser_patch.add_argument('--atmosphere_type', type=str, choices=['aerobic', 'anaerobic'])
    parser_patch.add_argument('--output_fp', type=pathlib.Path)
    parser_patch.add_argument('--biomass_reaction_id', type=str, default='BIOMASS_')
    parser_patch.add_argument('--memote_report_fp', type=pathlib.Path)
    parser_patch.add_argument('-h', '--help', action='store_true')

    parser_fba = subparsers.add_parser('fba', add_help=False)
    parser_fba.add_argument('--model_fp', type=pathlib.Path)
    parser_fba.add_argument('--fba_open_value', type=float, default=-1000)
    parser_fba.add_argument('--fba_spec_fp', type=pathlib.Path)
    parser_fba.add_argument('--output_fp', type=pathlib.Path)
    parser_fba.add_argument('-h', '--help', action='store_true')

    args = parser.parse_args()
    check_arguments(args)
    return args


def check_arguments(args):
    # pylint: disable=no-else-continue,too-many-branches
    if args.help:
        print(help_text(args.command), file=sys.stdout)
        sys.exit(0)
    if args.version:
        print(f'{__program_name__} {__version__}', file=sys.stdout)
        sys.exit(0)

    # Check we have required arguments, this is purposely decouped from argparse
    required_args = {
        'draft_model': {
            'single': ('assembly_fp', 'ref_model_fp', 'output_fp'),
            'exactly_one': (
                ('ref_genbank_fp', 'ref_proteins_fp'),
                ('ref_genbank_fp', 'ref_genes_fp'),
            ),
            'all': (
                ('ref_proteins_fp', 'ref_genes_fp'),
            ),
        },
        'patch_model': {
            'single': ('draft_model_fp', 'ref_model_fp', 'patch_fp', 'output_fp'),
        },
        'fba': {
            'single': ('model_fp', 'fba_spec_fp', 'output_fp'),
        },
    }

    if not args.command:
        msg = f'{__program_name__}: error: you must provide a command to execute'
        print(help_text(args.command), file=sys.stderr)
        print(msg, file=sys.stderr)
        sys.exit(1)
    else:
        assert args.command in required_args

    args_error_msgs = list()
    for req_type, required_set in required_args[args.command].items():
        for required_item in required_set:
            if req_type == 'single':
                if not args.__dict__[required_item]:
                    msg = f'--{required_item} is missing'
                    args_error_msgs.append(msg)
            elif req_type == 'exactly_one':
                args_present = sum(args.__dict__[arg] is not None for arg in required_item)
                if args_present != 1:
                    msg_part = ', '.join(f'--{arg}' for arg in required_item[:-1])
                    msg_part = f'{msg_part} or --{required_item[-1]}'
                    if args_present == 0:
                        msg = f'{msg_part} is required'
                    elif args_present > 1:
                        msg = f'can only specify one of {msg_part}'
                    args_error_msgs.append(msg)
            elif req_type == 'all':
                args_present = list()
                args_missing = list()
                for arg in required_item:
                    if args.__dict__[arg] is None:
                        args_missing.append(arg)
                    else:
                        args_present.append(arg)
                if args_present and args_missing:
                    args_missing_str = ', '.join(f'--{arg}' for arg in args_missing)
                    args_present_str = ', '.join(f'--{arg}' for arg in args_present)
                    args_error_msgs.append(f'the use of {args_present_str} also requires {args_missing_str}')
            else:
                print(f'error: argument requirement type: {req_type}')
                sys.exit(1)

    if args_error_msgs:
        msg_initial = f'{__program_name__}: error: the following argument errors were found:'
        print(help_text(args.command), file=sys.stderr)
        print(msg_initial, file=sys.stderr)
        for msg in args_error_msgs:
            print('\t', msg, sep='', file=sys.stderr)
        sys.exit(1)

    # Check all input file objects exist
    for arg, value in args.__dict__.items():
        if not value or arg in {'output_fp', 'memote_report_fp'}:
            continue
        elif isinstance(value, pathlib.Path):
            if not value.exists():
                print(f'{__program_name__}: error: input {value} does not exist', file=sys.stderr)
                sys.exit(1)
    # Check that output directory in output filepath exists
    if not args.output_fp.parent.exists():
        print(f'Output directory {args.output_fp.parent} for --output_fp does not exist', file=sys.stderr)
        sys.exit(1)


def help_text(command):
    info_text = f'\nProgram: {__program_name__}\nVersion: {__version__}\n'
    if not command:
        help_text_str = (
            f'Usage: {__program_name__} <command> [options]\n\n'
            'Commands:\n'
            '  draft_model                 Create a draft model\n'
            '  patch_model                 Patch a draft model\n'
            '  fba                         Simulate growth on media with FBA\n\n'
            f'For more information about commands, run: {__program_name__} <command> --help\n'
        )
    elif command == 'draft_model':
        help_text_str = (
            f'Usage: {__program_name__} {command} [options]\n'
            'Options:\n'
            '  --assembly_fp FILE          Isolate assembly filepath (GenBank)\n'
            '  --ref_genbank_fp FILE       Reference genbank filepath\n'
            '  --ref_proteins_fp FILE      Reference proteins filepath (FASTA)\n'
            '  --ref_genes_fp FILE         Reference genes filepath (FASTA)\n'
            '  --ref_model_fp FILE         Reference model filepath (JSON, XML [SMBL v3.1])\n'
            '  --min_coverage FLOAT        Alignment minimum coverage [default: 25]\n'
            '  --min_pident FLOAT          Alignment minimum percentage identity [default: 80]\n'
            '  --min_ppos FLOAT            Alignment minimum percentage positive matches\n'
            '  --media_type STR            Media type used to validate model [default: m9]\n'
            '  --atmosphere_type STR       Atmosphere type used to validate model\n'
            '  --biomass_reaction_id STR   Identifier of the biomass reaction [default: BIOMASS_]\n'
            '  --memote_report_fp FILE     MEMOTE report output filepath\n'
            '  --output_fp FILE            Output filepath\n'
            '\nOther:\n'
            '  --no_reannotation           Do not reannotate genbank file\n'
        )
    elif command == 'patch_model':
        help_text_str = (
            f'Usage: {__program_name__} {command} [options]\n'
            'Options:\n'
            '  --draft_model_fp FILE       Isolate model filepath\n'
            '  --ref_model_fp FILE         Reference model filepath\n'
            '  --patch_fp FILE             Patch file (JSON, XML [SMBL v3.1])\n'
            '  --media_type STR            Media type used to validate model [default: m9]\n'
            '  --atmosphere_type STR       Atmosphere type used to validate model\n'
            '  --biomass_reaction_id STR   Identifier of the biomass reaction [default: BIOMASS_]\n'
            '  --memote_report_fp FILE     MEMOTE report output filepath\n'
            '  --output_fp FILE            Output filepath\n'
        )
    elif command == 'fba':
        help_text_str = (
            f'Usage: {__program_name__} {command} [options]\n'
            'Options:\n'
            '  --model_fp FILE             Isolate model filepath\n'
            '  --fba_open_value FLOAT      Open reaction value to use in FBA [default: -1000]\n'
            '  --fba_spec_fp FILE          FBA spec filepath (JSON, XML [SMBL v3.1])\n'
            '  --output_fp FILE            Output filepath\n'
        )
    else:
        assert False
    return f'{info_text}\n{help_text_str}'
