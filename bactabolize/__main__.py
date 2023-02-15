import pathlib
import tempfile
import sys


from . import annotate
from . import arguments
from . import configuration
from . import draft_model
from . import model_fba
from . import model_sgk
from . import patch_model
from . import util


def entry():
    # Get command line arguments
    args = arguments.parse()

    # Execute workflows
    if args.command == 'draft_model':
        config = configuration.ConfigDraftModel(args)
        run_draft_model(config)
    elif args.command == 'patch_model':
        config = configuration.ConfigPatchModel(args)
        patch_model.run(config)
    elif args.command == 'fba':
        config = configuration.ConfigFba(args)
        model_fba.run(config)
    elif args.command == 'sgk':
        config = configuration.ConfigSgk(args)
        model_sgk.run(config)
    else:
        assert False


def run_draft_model(config):
    # pylint: disable=consider-using-with
    # Get input assembly format and convert if needed
    dh = tempfile.TemporaryDirectory()
    assembly_filetype = util.determine_assembly_filetype(config.assembly_fp)
    # If we have a FASTA input, require that we annotate
    if assembly_filetype == 'fasta' and config.no_reannotation:
        print('error: cannot specify --no_reannotation with a FASTA input assembly', file=sys.stderr)
        sys.exit(1)
    # Run annotation if requested
    if not config.no_reannotation:
        config.assembly_genbank_fp = config.output_fp.parent / f'{config.output_fp.stem}.gbk'
        annotate.run(config.assembly_fp, config.assembly_genbank_fp)
    else:
        config.assembly_genbank_fp = config.assembly_fp
    # If model is provided as a genbank, convert to FASTA
    if config.ref_genbank_fp:
        config.model_ref_genes_fp = pathlib.Path(dh.name, 'ref_genes.fasta')
        config.model_ref_proteins_fp = pathlib.Path(dh.name, 'ref_proteins.fasta')
        util.write_genbank_coding(config.ref_genbank_fp, config.model_ref_genes_fp, seq_type='nucl')
        util.write_genbank_coding(config.ref_genbank_fp, config.model_ref_proteins_fp, seq_type='prot')
    else:
        config.model_ref_genes_fp = config.ref_genes_fp
        config.model_ref_proteins_fp = config.ref_proteins_fp
    # Create draft model
    config.model_output_fp = config.output_fp.parent / f'{config.output_fp.stem}_model.json'
    config.model = util.read_model_and_check(
        config.ref_model_fp,
        config.model_ref_genes_fp,
        config.model_ref_proteins_fp,
    )
    config.alignment_thresholds = {
        'min_coverage': config.min_coverage,
        'min_pident': config.min_pident,
        'min_ppos': config.min_ppos,
    }
    draft_model.run(config)
    # Explicitly remove temporary directory
    dh.cleanup()


if __name__ == '__main__':
    entry()
