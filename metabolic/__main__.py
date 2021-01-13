import pathlib
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
    # If model is provided as a genbank, convert to FASTA
    if args.ref_genbank_fp:
        ref_genes, ref_proteins = create_fasta_from_genbank(args.ref_genbank_fp)
        ref_genes_fp = pathlib.Path(dh.name, 'ref_genes.fasta')
        ref_proteins_fp = pathlib.Path(dh.name, 'ref_proteins.fasta')
        util.write_dict_to_fasta(ref_genes, ref_genes_fp)
        util.write_dict_to_fasta(ref_proteins, ref_proteins_fp)
    else:
        ref_genes_fp = args.ref_genes_fp
        ref_proteins_fp = args.ref_proteins_fp
    # Create draft model
    output_fp = args.output_fp.parent / f'{args.output_fp.stem}_model.json'
    model = util.read_model_and_check(args.ref_model_fp, ref_genes_fp, ref_proteins_fp)
    draft_model.run(assembly_genbank_fp, ref_genes_fp, ref_proteins_fp, model, output_fp)
    # Explicitly remove temporary directory
    dh.cleanup()


import Bio.SeqIO

def create_fasta_from_genbank(genbank_fp):
    genes = dict()
    proteins = dict()
    locus_tag_dups = dict()
    with genbank_fp.open('r') as fh:
        for record in Bio.SeqIO.parse(fh, 'genbank'):
            for feature in record.features:
                if feature.type != 'CDS':
                    continue
                [locus_tag] = feature.qualifiers['locus_tag']
                if locus_tag in genes:
                    n = locus_tag_dups.get(locus_tag, 1)
                    locus_tag_dups[locus_tag] = n + 1
                    msg = f'warning: duplicate locus tag {locus_tag}, renaming to {locus_tag}_{n}'
                    print(msg, file=sys.stderr)
                    locus_tag = f'{locus_tag}_{n}'
                genes[locus_tag] = feature.extract(record.seq)
                if 'translation' not in feature.qualifiers:
                    protein_seq = genes[locus_tag].translate(stop_symbol='')
                else:
                    [protein_seq] = feature.qualifiers['translation']
                proteins[locus_tag] = protein_seq
    return genes, proteins


if __name__ == '__main__':
    entry()
