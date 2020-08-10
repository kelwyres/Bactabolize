import tempfile
import sys


from . import annotate
from . import arguments
from . import assembly_stats
from . import model_fba
from . import draft_model
from . import patch_model
from . import util


def entry():
    # Get command line arguments
    args = arguments.parse()

    # Execute workflows
    if args.command is None:
        run_complete_workflow(args)
    elif args.command == 'assembly_qc':
        assembly_stats.run(args.assembly_fp, args.output_fp)
    elif args.command == 'annotate':
        annotate.run(args.assembly_fp, args.prodigal_model_fp, args.output_fp)
    elif args.command == 'draft_model':
        model = util.read_model_and_check(args.ref_model_fp, args.ref_genbank_fp)
        draft_model.run(args.assembly_fp, args.ref_genbank_fp, model, args.output_fp)
    elif args.command == 'patch_model':
        patch_model.run(args.draft_model_fp, args.ref_model_fp, args.patch_fp, args.output_fp)
    elif args.command == 'model_fba':
        model_fba.run(args.model_fp, args.fba_spec_fp, args.output_fp)
    else:
        assert False


def run_complete_workflow(args):
    # Get input assembly format and convert if needed
    dh = tempfile.TemporaryDirectory()
    assembly_filetype = util.determine_assembly_filetype(args.assembly_fp)
    # If we have a FASTA input, require that we annotate
    if assembly_filetype == 'fasta' and args.no_reannotation:
        print('error: cannot specify --no_reannotation with a FASTA input assembly', file=sys.stderr)
        exit(1)
    # Run optional stages if requested
    if not args.no_qc:
        stats_fp = args.output_dir / f'{args.assembly_fp.stem}_stats.tsv'
        assembly_stats.run(args.assembly_fp, stats_fp)
    if not args.no_reannotation:
        assembly_genbank_fp = args.output_dir / f'{args.assembly_fp.stem}.gbk'
        annotate.run(args.assembly_fp, args.prodigal_model_fp, assembly_genbank_fp)
    else:
        assembly_genbank_fp = args.assembly_fp
    # Explicitly remove temporary directory
    dh.cleanup()
    # Create draft model
    draft_model_fp = args.output_dir / f'{args.assembly_fp.stem}_model.json'
    model = util.read_model_and_check(args.ref_model_fp, args.ref_genbank_fp)
    draft_model.run(assembly_genbank_fp, args.ref_genbank_fp, model, draft_model_fp)
    # Run FBA
    if not args.no_fba:
        fba_results_fp = args.output_dir / f'{args.assembly_fp.stem}_fba.tsv'
        model_fba.run(draft_model_fp, args.fba_spec_fp, fba_results_fp)


if __name__ == '__main__':
    entry()
