from . import util



BlastFormat = {'qseqid': str,
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
               'mismatch': int,
               'gaps': int}


class BlastResult:

    def __init__(self, *values):
        for (attr, attr_type), value in zip(BlastFormat.items(), values):
            setattr(self, attr, attr_type(value))

    def __str__(self):
        data_gen = (getattr(self, attr) for attr in BlastFormat)
        return '\t'.join(str(data) for data in data_gen)


def run_blastp(query_fp, subject_fp):
    # Create database
    command_db = f'makeblastdb -in {subject_fp} -out {subject_fp} -dbtype prot'
    util.execute_command(command_db)
    # Run alignment
    command_opts = f'-evalue 0.001 -outfmt \'6 {" ".join(BlastFormat)}\''
    command_run = f'blastp -db {subject_fp} -query {query_fp} {command_opts}'
    result = util.execute_command(command_run)
    return parse_results(result.stdout)


def run_blastn(query_fp, subject_fp):
    # Create database
    command_db = f'makeblastdb -in {subject_fp} -out {subject_fp} -dbtype nucl'
    util.execute_command(command_db)
    # Run alignment
    command_opts = f'-evalue 0.001 -outfmt \'6 {" ".join(BlastFormat)}\''
    command_run = f'blastn -db {subject_fp} -query {query_fp} {command_opts}'
    result = util.execute_command(command_run)
    return parse_results(result.stdout)


def filter_results(results, *, min_coverage, min_pident):
    results_filtered = dict()
    for qseqid, hits in results.items():
        for hit in hits:
            # Filter on coverage and pident
            if hit.length / hit.qlen * 100 < min_coverage:
                continue
            if hit.pident < min_pident:
                continue
            if hit.qseqid not in results_filtered:
                results_filtered[hit.qseqid] = list()
            results_filtered[hit.qseqid].append(hit)
    return results_filtered


def parse_results(results):
    query_hits = dict()
    line_token_gen = (line.split() for line in results.rstrip().split('\n'))
    for line_tokens in line_token_gen:
        hit = BlastResult(*line_tokens)
        if hit.qseqid not in query_hits:
            query_hits[hit.qseqid] = list()
        query_hits[hit.qseqid].append(hit)
    return query_hits
