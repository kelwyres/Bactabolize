import re
import tempfile


import Bio.Seq
import Bio.SeqFeature
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqRecord


# Backwards compatible Bio.Alphabet import
try:
    import Bio.Alphabet.IUPAC as bio_alphabet
except ImportError:
    bio_alphabet = None


from . import util


def run(assembly_fp, output_fp, *, model_fp=None):
    print('\n========================================')
    print('running annotation')
    # Check assembly filetype and convert if needed
    dh = tempfile.TemporaryDirectory()
    assembly_filetype = util.determine_assembly_filetype(assembly_fp)
    if assembly_filetype == 'genbank':
        assembly_genbank_fp = assembly_fp
        assembly_fasta_fp = util.write_genbank_seq(assembly_fp, dh.name)
    elif assembly_filetype == 'fasta':
        assembly_genbank_fp = None
        assembly_fasta_fp = assembly_fp
    print('========================================')
    prodigal_data = run_prodigal(assembly_fasta_fp, model_fp)
    prodigal_orfs = parse_prodigal_output(prodigal_data)

    print(f'Found {len(prodigal_orfs)} open-reading frames')
    genbank_records = create_genbank(prodigal_orfs, assembly_fasta_fp)
    with output_fp.open('w') as fh:
        Bio.SeqIO.write(genbank_records, fh, 'genbank')

    # Explicitly remove temporary directory
    dh.cleanup()

    # Match ORFs if we have an input genbank
    if assembly_filetype == 'genbank':
        match_existing_orfs_updated_annotations(output_fp, assembly_genbank_fp)


def run_prodigal(assembly_fp, model_fp):
    command = f'prodigal -f sco -i {assembly_fp} -m -o /dev/null -d /dev/stdout'
    if model_fp:
        command = f'{command} -t {model_fp}'
    result = util.execute_command(command)
    # Prodigal includes \r from FASTAs, causing problems with the output. Remove \r here
    return result.stdout.replace('\r', '')


def parse_prodigal_output(prodigal_data):
    prodigal_result_re = re.compile(r'^>(.+)_([0-9]+) # ([0-9]+) # ([0-9]+) # (-?1) # .+partial=([0-1]+).+$')
    orfs = list()
    for line in prodigal_data.rstrip().split('\n'):
        if not line.startswith('>'):
            continue
        orf_data = prodigal_result_re.match(line).groups()
        orfs.append(orf_data)
    return orfs


def create_genbank(orfs, assembly_fp):
    # Create unannotated gebnank records
    genbank_records = dict()
    with assembly_fp.open('r') as fh:
        for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh):
            contig_id = desc.split(' ', maxsplit=1)[0]
            assert contig_id not in genbank_records
            # Backwards compatiblity Seq init
            if bio_alphabet:
                sequence_record = Bio.Seq.Seq(''.join(seq), bio_alphabet)
            else:
                sequence_record = Bio.Seq.Seq(''.join(seq))
            genbank_records[contig_id] = Bio.SeqRecord.SeqRecord(
                seq=sequence_record,
                id=contig_id,
                name=contig_id,
                annotations={'molecule_type': 'DNA'}
            )
    # Annotate records with prodigal ORFs
    for contig, orf_n, posl_str, posr_str, strand_str, partial_str in orfs:
        if strand_str == '1':
            strand = 1
        elif strand_str == '-1':
            strand = -1
        else:
            assert False
        quals = {'gene': orf_n, 'locus_tag': orf_n, 'note': get_qual_note(partial_str)}
        feature_loc = Bio.SeqFeature.FeatureLocation(start=int(posl_str)-1, end=int(posr_str), strand=strand)
        feature = Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=quals)
        genbank_records[contig].features.append(feature)
    return [record for record in genbank_records.values()]


def get_qual_note(partial_str):
    if partial_str == '00':
        note = 'left:complete;right:complete'
    elif partial_str == '01':
        note = 'left:complete;right:unfinished'
    elif partial_str == '10':
        note = 'left:unfinished;right:complete'
    elif partial_str == '11':
        note = 'left:unfinished;right:unfinished'
    else:
        assert False
    return note


def match_existing_orfs_updated_annotations(new_fp, existing_fp, overlap_min=0.80):
    # Get features and create list of start and end objects for each
    features_new = collect_all_features(new_fp)
    features_existing = collect_all_features(existing_fp)
    contig_positions_new = create_positions(features_new, 'new')
    contig_positions_existing = create_positions(features_existing, 'existing')

    # Match features
    contig_features_updated = dict()
    for contig in set(contig_positions_new) | set(contig_positions_existing):
        # Find overlaps
        positions = contig_positions_new[contig] + contig_positions_existing[contig]
        features_matched = discover_overlaps(positions, overlap_min)
        # Discover those not matched
        features_matched_flat = set()
        for features in features_matched:
            features_matched_flat.update(features)
        new_unmatched = set(features_new[contig]).difference(features_matched_flat)
        existing_unmatched = set(features_existing[contig]).difference(features_matched_flat)
        # For each matched update bounds update locus tag, product, gene (if present) to match existing
        quals = ('locus_tag', 'product', 'gene')
        features_updated = list()
        for feature_new, feature_existing in features_matched:
            for qual in quals:
                if qual not in feature_existing.qualifiers:
                    continue
                feature_new.qualifiers[qual] = feature_existing.qualifiers[qual]
            features_updated.append(feature_new)
        # Add existing ORFs that had no match
        features_updated.extend(new_unmatched)
        features_updated.extend(existing_unmatched)
        contig_features_updated[contig] = sorted(features_updated, key=lambda k: k.location.start)
        print(f'\n{contig}:')
        print(f'\t{len(features_existing[contig])} features in existing genbank')
        print(f'\t{len(features_new[contig])} features in reannotation')
        print(f'\t{len(features_matched)} reannotated features matched to existing')
        print(f'\t{len(existing_unmatched)} existing features unmatched')
        print(f'\t{len(new_unmatched)} re-annotated features unmatched')
        print(f'\t{len(features_updated)} total features')
    # Update new genbank with new feature set
    update_genbank_annotations(new_fp, contig_features_updated)


def collect_all_features(fp):
    features = dict()
    with fp.open('r') as fh:
        for record in Bio.SeqIO.parse(fh, 'genbank'):
            features[record.name] = list(util.iterate_coding_features(record))
    return features


def create_positions(contig_features, source):
    positions = dict()
    for contig, features in contig_features.items():
        positions[contig] = list()
        for feature in features:
            start = {'position': feature.location.start, 'type': 'start', 'feature': feature, 'source': source}
            end = {'position': feature.location.end, 'type': 'end', 'feature': feature, 'source': source}
            positions[contig].extend((start, end))
    return positions


def discover_overlaps(positions, overlap_min):
    in_new = list()
    in_existing = list()
    features_matched = set()
    for position in sorted(positions, key=lambda k: k['position']):
        # Add features we're entering and remove those we're exiting
        if position['type'] == 'start':
            if position['source'] == 'new':
                in_new.append(position['feature'])
            elif position['source'] == 'existing':
                in_existing.append(position['feature'])
        elif position['type'] == 'end':
            if position['source'] == 'new':
                in_new.remove(position['feature'])
            elif position['source'] == 'existing':
                in_existing.remove(position['feature'])
        if not (in_new and in_existing):
            continue
        # Iterate new/existing features that overlap
        for feature_new in in_new:
            for feature_existing in in_existing:
                if feature_new.strand != feature_existing.strand:
                    continue
                if (feature_new, feature_existing) in features_matched:
                    continue
                # Get overlap
                start = max(feature_new.location.start, feature_existing.location.start)
                end = min(feature_new.location.end, feature_existing.location.end)
                overlap_count = end - start
                overlap_new = overlap_count / len(feature_new)
                overlap_existing = overlap_count / len(feature_existing)
                # Record if overlap over threshold
                if overlap_new > overlap_min and overlap_existing > overlap_min:
                    # Update note to include overlap information
                    [note_new] = feature_new.qualifiers['note']
                    feature_new.qualifiers['note'][0] = f'{note_new};overlap:{overlap_new:.2f}'
                    features_matched.add((feature_new, feature_existing))
    return features_matched


def update_genbank_annotations(filepath, features_updated):
    # Read in data then update
    with filepath.open('r') as fh:
        records = {record.name: record for record in Bio.SeqIO.parse(fh, 'genbank')}
    for contig, record in records.items():
        record.features = features_updated[contig]
    # Write
    with filepath.open('w') as fh:
        Bio.SeqIO.write(records.values(), fh, 'genbank')
