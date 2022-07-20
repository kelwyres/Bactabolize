import contextlib
import pathlib
import subprocess
import sys


import Bio.SeqIO
import Bio.SeqIO.FastaIO
import cobra.io
from cobra.io import read_sbml_model


def read_model_and_check(model_fp, genes_fp, proteins_fp):
    print('\n========================================')
    print('reading reference model data')
    print('========================================')
    with model_fp.open('r') as fh:
        if model_fp.suffix == '.json':
            model = cobra.io.load_json_model(fh)
        elif model_fp.suffix == '.xml':
            model = read_sbml_model(fh)
    # Collect genes/proteins
    model_genes = {gene.id for gene in model.genes}
    with genes_fp.open('r') as fh:
        genes = {desc for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
    with proteins_fp.open('r') as fh:
        proteins = {desc for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
    # Check gene counts
    check_genes_proteins(model_genes, genes, 'genes')
    check_genes_proteins(model_genes, proteins, 'proteins')
    return model


def check_genes_proteins(model_genes, other_genes, other_type):
    missing = model_genes.difference(other_genes)
    if missing:
        plurality = 'entry' if len(missing) == 1 else 'entries'
        if len(missing) > 10:
            msg = f'warning: could not find model {len(missing)} {plurality} in reference {other_type}'
            print(msg, file=sys.stderr)
        else:
            msg = f'warning: could not find model {len(missing)} {plurality} in reference {other_type}:'
            print(msg, ', '.join(missing), file=sys.stderr)


def determine_assembly_filetype(filepath):
    accepted_types = ('genbank', 'fasta')
    expected_map = {
        'fasta': 'fasta',
        'fna': 'fasta',
        'fa': 'fasta',
        'gbk': 'genbank',
        'gb': 'genbank',
    }
    file_extension = filepath.suffix[1:]
    expected_filetype = expected_map.get(file_extension, 'unknown')
    for filetype in accepted_types:
        with filepath.open('r') as fh:
            if list(Bio.SeqIO.parse(fh, filetype)):
                break
    else:
        print(f'error: could not parse {filepath} as either GenBank or FASTA format', file=sys.stderr)
        sys.exit(1)
    if expected_filetype != filetype:
        msg = f'warning: parsed {filepath} as {filetype}'
        if expected_filetype == 'unknown':
            print(f'{msg} but had an unknown file extension (.{file_extension})', file=sys.stderr)
        else:
            print(f'{msg} but expected {expected_filetype}', file=sys.stderr)
    return filetype


def write_genbank_coding(filepath, output_fp, *, seq_type):
    # pylint: disable=too-many-branches
    with contextlib.ExitStack() as stack:
        fh_in = stack.enter_context(filepath.open('r'))
        fh_out = stack.enter_context(output_fp.open('w'))
        # Iterate coding features
        for record in Bio.SeqIO.parse(fh_in, 'genbank'):
            for feature in iterate_coding_features(record):
                nucleotide_seq = feature.extract(record.seq)
                if seq_type == 'prot':
                    if 'translation' not in feature.qualifiers:
                        # Append trailing N if we have a partial codon
                        seq_n = (3 - len(nucleotide_seq)) % 3
                        nucleotide_seq = nucleotide_seq + 'N' * seq_n
                        seq = nucleotide_seq.translate(stop_symbol='')
                    else:
                        [seq] = feature.qualifiers['translation']
                elif seq_type == 'nucl':
                    seq = nucleotide_seq
                else:
                    assert False
                # Write to disk
                [locus_tag] = feature.qualifiers['locus_tag']
                desc = f'>{locus_tag}'
                seq_lines = [seq[i : i + 80] for i in range(0, len(seq), 80)]
                print(desc, file=fh_out)
                print(*seq_lines, sep='\n', file=fh_out)
    return output_fp


def write_genbank_seq(filepath, dirpath):
    fasta_fp = str()
    output_fp = pathlib.Path(dirpath, f'{filepath.stem}.fasta')
    with contextlib.ExitStack() as stack:
        fh_in = stack.enter_context(filepath.open('r'))
        fh_out = stack.enter_context(output_fp.open('w'))
        fasta_fp = fh_out.name
        for record in Bio.SeqIO.parse(fh_in, 'genbank'):
            desc = f'>{record.name}'
            seq_lines = [record.seq[i : i + 80] for i in range(0, len(record.seq), 80)]
            print(desc, file=fh_out)
            print(*seq_lines, sep='\n', file=fh_out)
    return pathlib.Path(fasta_fp)


def iterate_coding_features(record):
    for feature in record.features:
        if feature.type != 'CDS':
            continue
        assert 'locus_tag' in feature.qualifiers
        yield feature


def extract_nucleotides_from_ref(hit, fasta):
    if hit.sstart > hit.send:
        seq = fasta[hit.sseqid][hit.send - 1 : hit.sstart].reverse_complement()
    else:
        seq = fasta[hit.sseqid][hit.sstart - 1 : hit.send]
    return seq


def execute_command(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')
    if result.returncode != 0:
        print('Failed to run command:', result.args, file=sys.stderr)
        print('stdout:', result.stdout, file=sys.stderr)
        print('stderr:', result.stderr, file=sys.stderr)
        sys.exit(1)
    return result
