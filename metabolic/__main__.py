import pathlib
import tempfile
import sys


from . import annotate
from . import arguments
from . import assembly_stats
from . import model_fba
from . import draft_model
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
        annotate.run(args.assembly_fp, args.model_fp, args.output_fp)
    elif args.command == 'draft_model':
        model = util.read_model_and_check(args.ref_model_fp, args.ref_gbk_fp)
        draft_model.run(args.assembly_fp, args.ref_gbk_fp, model, args.output_fp)
    elif args.command == 'model_fba':
        model_fba.run(args.model_fp)
    else:
        assert False


def run_complete_workflow(args):
    # TODO determine input assembly filetype and handle, currently assuming genbank
    dh = tempfile.TemporaryDirectory()  # fs directory deleted this is out of scope
    assembly_fasta_fp = util.write_gbk_to_fasta(args.assembly_fp, dh.name)
    # TODO: apply some QC to these stats
    if not args.no_qc:
        # TODO: using named arg to show that we need to set output_fp to something later
        assembly_stats.run(assembly_fasta_fp, output_fp=None)
    # TODO: clean this up - don't use assembly_fp as a name for easier op below
    if not args.no_annotation:
        assembly_fp = pathlib.Path(dh.name, f'{args.assembly_fp.stem}_reanno.gbk')
        annotate.run(assembly_fasta_fp, args.prodigal_model_fp, assembly_fp)
    else:
        assembly_fp = args.assembly_fp
    # Create draft model
    draft_model_fp = pathlib.Path(args.output_dir / f'{assembly_fp.stem}_model.json')
    model = util.read_model_and_check(args.ref_model_fp, args.ref_gbk_fp)
    draft_model.run(assembly_fp, args.ref_gbk_fp, model, draft_model_fp)
    # Simulate growth on media
    model_fba.run(draft_model_fp)



if __name__ == '__main__':
    entry()
