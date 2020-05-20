import contextlib
import subprocess
import sys
import tempfile


import Bio.SeqIO


def write_gbk_sequence(filepath, dirpath, retain=None, *, seq_type):
    assert isinstance(retain, set) or retain is None
    assert seq_type in {'prot', 'nucl'}
    sequence_fp = str()
    with contextlib.ExitStack() as stack:
        fh_in = stack.enter_context(filepath.open('r'))
        fh_out = stack.enter_context(tempfile.NamedTemporaryFile('wt', delete=False, dir=dirpath))
        sequence_fp = fh_out.name
        # Iterate coding features, record seen genes if retaining
        genes_seen = set()
        for record in Bio.SeqIO.parse(fh_in, 'genbank'):
            for feature, nucleotide_seq in iterate_coding_features(record):
                # Get appropriate output sequence type
                if seq_type == 'prot':
                    seq = nucleotide_seq.translate()
                elif seq_type == 'nucl':
                    seq = nucleotide_seq
                else:
                    assert False
                [locus_tag] = feature.qualifiers['locus_tag']
                # Skip if we do not intend to retain
                if retain and locus_tag not in retain:
                    continue
                elif retain:
                    genes_seen.add(locus_tag)
                # Write to disk
                desc = f'>{locus_tag}'
                seq_lines = [seq[i:i+80] for i in range(0, len(seq), 80)]
                print(desc, file=fh_out)
                print(*seq_lines, sep='\n', file=fh_out)
    # TODO: move this to pre-flight reference and model validation
    # If we were retaining genes ensure all were found, otherwise warn
    if retain and len(retain ^ genes_seen) > 0:
        missing_genes = retain ^ genes_seen
        plurality = 'gene' if len(missing_genes) == 1 else 'genes'
        gene_list_str = [f'\t{gene}' for gene in missing_genes]
        print(f'warning: could not find {len(missing_genes)} {plurality}:', file=sys.stderr)
        print(*gene_list_str, sep='\n', file=sys.stderr)
    return sequence_fp


def iterate_coding_features(record):
    for feature in record.features:
        if feature.type != 'CDS':
            continue
        assert 'locus_tag' in feature.qualifiers
        yield feature, feature.extract(record.seq)


def write_gbk_to_fasta(filepath, dirpath):
    fasta_fp = str()
    with contextlib.ExitStack() as stack:
        fh_in = stack.enter_context(filepath.open('r'))
        fh_out = stack.enter_context(tempfile.NamedTemporaryFile('wt', delete=False, dir=dirpath))
        fasta_fp = fh_out.name
        for record in Bio.SeqIO.parse(fh_in, 'genbank'):
            desc = f'>{record.name}'
            seq_lines = [record.seq[i:i+80] for i in range(0, len(record.seq), 80)]
            print(desc, file=fh_out)
            print(*seq_lines, sep='\n', file=fh_out)
    return fasta_fp


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
