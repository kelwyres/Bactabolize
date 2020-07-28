import contextlib
import pathlib
import subprocess
import sys
import tempfile


import Bio.SeqIO
import cobra.io


def read_model_and_check(model_fp, genbank_fp):
    print('\n========================================')
    print('reading reference model and genbank')
    print('========================================')
    with model_fp.open('r') as fh:
        model = cobra.io.load_json_model(fh)
    # Get genbank genes
    genbank_genes = set()
    with genbank_fp.open('r') as fh:
        for record in Bio.SeqIO.parse(fh, 'genbank'):
            for feature in record.features:
                if feature.type != 'CDS':
                    continue
                [locus_tag] = feature.qualifiers['locus_tag']
                genbank_genes.add(locus_tag)
    # Check all model genes are in genbank
    model_genes = {gene.id for gene in model.genes}
    missing_genes = model_genes.difference(genbank_genes)
    if missing_genes:
        if len(missing_genes) > 10:
            msg = f'warning: could not find {len(missing_genes)} model {plurality} in reference genbank'
            print(msg, file=sys.stderr)
        else:
            plurality = 'gene' if len(missing_genes) == 1 else 'genes'
            msg = f'warning: could not find {len(missing_genes)} model {plurality} in reference genbank:'
            print(msg, ', '.join(missing_genes), file=sys.stderr)
    return model


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


def write_genbank_coding_sequence(filepath, dirpath, retain=None, *, seq_type):
    assert isinstance(retain, set) or retain is None
    assert seq_type in {'prot', 'nucl'}
    sequence_fp = str()
    output_fp = pathlib.Path(dirpath, f'{filepath.stem}.fasta')
    with contextlib.ExitStack() as stack:
        fh_in = stack.enter_context(filepath.open('r'))
        fh_out = stack.enter_context(output_fp.open('w'))
        sequence_fp = fh_out.name
        # Iterate coding features
        for record in Bio.SeqIO.parse(fh_in, 'genbank'):
            for feature, nucleotide_seq in iterate_coding_features(record):
                # Get appropriate output sequence type
                if seq_type == 'prot':
                    # Append trailing N if we have a partial codon
                    seq_n = (3 - len(nucleotide_seq)) % 3
                    nucleotide_seq = nucleotide_seq + 'N' * seq_n
                    seq = nucleotide_seq.translate()
                elif seq_type == 'nucl':
                    seq = nucleotide_seq
                else:
                    assert False
                [locus_tag] = feature.qualifiers['locus_tag']
                # Skip if we do not intend to retain
                if retain and locus_tag not in retain:
                    continue
                # Write to disk
                desc = f'>{locus_tag}'
                seq_lines = [seq[i:i+80] for i in range(0, len(seq), 80)]
                print(desc, file=fh_out)
                print(*seq_lines, sep='\n', file=fh_out)
    return sequence_fp


def iterate_coding_features(record, sequence=True):
    for feature in record.features:
        if feature.type != 'CDS':
            continue
        assert 'locus_tag' in feature.qualifiers
        if sequence:
            yield feature, feature.extract(record.seq)
        else:
            yield feature


def write_genbank_to_fasta(filepath, dirpath):
    fasta_fp = str()
    output_fp = pathlib.Path(dirpath,  f'{filepath.stem}.fasta')
    with contextlib.ExitStack() as stack:
        fh_in = stack.enter_context(filepath.open('r'))
        fh_out = stack.enter_context(output_fp.open('w'))
        fasta_fp = fh_out.name
        for record in Bio.SeqIO.parse(fh_in, 'genbank'):
            desc = f'>{record.name}'
            seq_lines = [record.seq[i:i+80] for i in range(0, len(record.seq), 80)]
            print(desc, file=fh_out)
            print(*seq_lines, sep='\n', file=fh_out)
    return pathlib.Path(fasta_fp)


def extract_nucleotides_from_ref(hit, fasta):
    if hit.sstart > hit.send:
        seq = fasta[hit.sseqid][hit.send-1:hit.sstart].reverse_complement()
    else:
        seq = fasta[hit.sseqid][hit.sstart-1:hit.send]
    return seq


def execute_command(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')
    if result.returncode != 0:
        print('Failed to run command:', result.args, file=sys.stderr)
        print('stdout:', result.stdout, file=sys.stderr)
        print('stderr:', result.stderr, file=sys.stderr)
        sys.exit(1)
    return result
