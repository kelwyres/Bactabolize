import tempfile
import sys


from . import annotate
from . import arguments
from . import model_fba
from . import draft_model
from . import patch_model
from . import util


def entry():
    # Get command line arguments
    args = arguments.parse()

    # Execute workflows
    if args.command == 'draft_model':
        run_draft_model(args)
    elif args.command == 'patch_model':
        patch_model.run(args.draft_model_fp, args.ref_model_fp, args.patch_fp, args.output_fp)
    elif args.command == 'fba':
        model_fba.run(args.model_fp, args.fba_spec_fp, args.output_fp)
    else:
        assert False


def run_draft_model(args):
    # Get input assembly format and convert if needed
    dh = tempfile.TemporaryDirectory()
    assembly_filetype = util.determine_assembly_filetype(args.assembly_fp)
    # If we have a FASTA input, require that we annotate
    if assembly_filetype == 'fasta' and args.no_reannotation:
        print('error: cannot specify --no_reannotation with a FASTA input assembly', file=sys.stderr)
        exit(1)
    # Run annotation if requested
    if not args.no_reannotation:
        assembly_genbank_fp = args.output_fp.parent / f'{args.output_fp.stem}.gbk'
        annotate.run(args.assembly_fp, assembly_genbank_fp, model_fp=args.prodigal_model_fp)
    else:
        assembly_genbank_fp = args.assembly_fp
    # Explicitly remove temporary directory
    dh.cleanup()
    # Create draft model
    draft_model_fp = args.output_fp.parent / f'{args.output_fp.stem}_model.json'
    model = util.read_model_and_check(args.ref_model_fp, args.ref_genbank_fp)
    draft_model.run(assembly_genbank_fp, args.ref_genbank_fp, model, draft_model_fp)


if __name__ == '__main__':
    entry()
