import re


import Bio.Seq
import Bio.SeqFeature
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqRecord
import Bio.Alphabet.IUPAC


from . import util


def run(assembly_fp, model_fp, output_fp):
    print('\n==============================')
    print('running annotation')
    print('==============================')
    #pickle_mode = 'read'
    pickle_mode = 'write'
    #pickle_mode = 'noop'
    prodigal_data = run_prodigal(assembly_fp, model_fp)

    # TEMP: store/load prodigal results
    import pickle
    if pickle_mode == 'write':
        with open(f'pickled/prodigal_{assembly_fp.stem}.bin', 'wb') as fh:
            pickle.dump(prodigal_data, fh)
    elif pickle_mode == 'read':
        with open(f'pickled/prodigal_{assembly_fp.stem}.bin', 'rb') as fh:
            prodigal_data = pickle.load(fh)
    elif pickle_mode == 'noop':
        pass
    else:
        assert False

    prodigal_orfs = parse_prodigal_output(prodigal_data)
    genbank_records = create_genbank(prodigal_orfs, assembly_fp)
    with output_fp.open('w') as fh:
        Bio.SeqIO.write(genbank_records, fh, 'genbank')


def run_prodigal(assembly_fp, model_fp):
    command = f'prodigal -f sco -i {assembly_fp} -m -t {model_fp}'
    result = util.execute_command(command)
    # Prodigal includes \r from FASTAs, causing problems with the output. Remove \r here
    return result.stdout.replace('\r', '')


def parse_prodigal_output(prodigal_data):
    prodigal_result_re = re.compile(r'^>[0-9]+_([0-9]+)_([0-9]+)_([-+])$')
    prodigal_contig_re = re.compile(r'^# Sequence.+?seqhdr="(.+?)"(?:;|$)')

    orfs = list()
    for line in prodigal_data.rstrip().split('\n'):
        if line.startswith('# Sequence Data'):
            contig = prodigal_contig_re.match(line).group(1)
        elif line.startswith('# Model Data'):
            continue
        else:
            orf_info = prodigal_result_re.match(line).groups()
            orfs.append((contig, orf_info))
    return orfs


def create_genbank(orfs, assembly_fp):
    # Create unannotated gebnank records
    genbank_records = dict()
    with assembly_fp.open('r') as fh:
        for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh):
            contig_id = desc.split(' ', maxsplit=1)[0]
            assert contig_id not in genbank_records
            sequence_record = Bio.Seq.Seq(''.join(seq), Bio.Alphabet.IUPAC.ambiguous_dna)
            genbank_records[contig_id] = Bio.SeqRecord.SeqRecord(seq=sequence_record, id=contig_id, name=contig_id)
    # Annotate records with prodigal ORFs
    gene_n = 0
    for contig, (start_str, end_str, strand_str) in orfs:
        gene_n += 1
        contig_id = contig.split(' ', maxsplit=1)[0]
        strand = +1 if strand_str == '+' else -1
        quals = {'gene': gene_n, 'locus_tag': gene_n}
        feature_loc = Bio.SeqFeature.FeatureLocation(start=int(start_str), end=int(end_str), strand=strand)
        feature = Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=quals)
        genbank_records[contig_id].features.append(feature)
    return [record for record in genbank_records.values()]
