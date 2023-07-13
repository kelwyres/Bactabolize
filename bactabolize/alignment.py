import pathlib


import filelock


from . import util


BlastFormat = {
    'qseqid': str,
    'sseqid': str,
    'qlen': int,
    'slen': int,
    'qstart': int,
    'qend': int,
    'sstart': int,
    'send': int,
    'length': int,
    'evalue': float,
    'bitscore': float,
    'pident': float,
    'nident': int,
    'ppos': float,
    'mismatch': int,
    'gaps': int,
}


class BlastResult:
    def __init__(self, *values):
        for (attr, attr_type), value in zip(BlastFormat.items(), values):
            setattr(self, attr, attr_type(value))

    def __str__(self):
        data_gen = (getattr(self, attr) for attr in BlastFormat)
        return '\t'.join(str(data) for data in data_gen)


def run_blastp(query_fp, subject_fp):
    # Create a database
    create_blast_database(subject_fp, 'prot')
    # Run alignment
    command_opts = f'-evalue 0.001 -outfmt \'6 {" ".join(BlastFormat)}\''
    command_run = f'blastp -db {subject_fp} -query {query_fp} {command_opts}'
    result = util.execute_command(command_run)
    if result.stdout != "":
        return parse_results(result.stdout)
    else:
        return {}


def run_blastn(query_fp, subject_fp):
    # Create database
    create_blast_database(subject_fp, 'nucl')
    # Run alignment
    command_opts = f'-evalue 0.001 -outfmt \'6 {" ".join(BlastFormat)}\''
    command_run = f'blastn -db {subject_fp} -query {query_fp} {command_opts}'
    result = util.execute_command(command_run)
    if result.stdout != "":
        return parse_results(result.stdout)
    else:
        return {}


def create_blast_database(subject_fp, db_type):
    # Check if we already have a database
    database_exts = ['pdb', 'psq', 'pto', 'ptf', 'pot', 'pin', 'phr']
    database_fps = [pathlib.Path(f'{subject_fp}.{ext}') for ext in database_exts]
    if any(not fp.exists() for fp in database_fps):
        # Create a lock to prevent race condition; multiple concurrent processes could eval to here
        with filelock.FileLock(f'{subject_fp}.lock', timeout=60):
            # A second check for multiple concurrent processes to ensure database is created only once
            if any(not fp.exists() for fp in database_fps):
                command_db = f'makeblastdb -in {subject_fp} -out {subject_fp} -dbtype {db_type}'
                util.execute_command(command_db)


def filter_results(results, *, min_coverage=None, min_pident=None, min_ppos=None):
    results_filtered = dict()
    for hits in results.values():
        for hit in hits:
            if min_coverage and hit.length / hit.qlen * 100 < min_coverage:
                continue
            if min_pident and hit.pident < min_pident:
                continue
            if min_ppos and hit.ppos < min_ppos:
                continue
            if hit.qseqid not in results_filtered:
                results_filtered[hit.qseqid] = list()
            results_filtered[hit.qseqid].append(hit)
    return results_filtered


def parse_results(results):
    # pylint: disable=no-member
    query_hits = dict()
    line_token_gen = (line.split() for line in results.rstrip().split('\n'))
    for line_tokens in line_token_gen:
        hit = BlastResult(*line_tokens)
        if hit.qseqid not in query_hits:
            query_hits[hit.qseqid] = list()
        query_hits[hit.qseqid].append(hit)
    return query_hits
