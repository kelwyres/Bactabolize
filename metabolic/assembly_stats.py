import math
import statistics


from Bio.SeqIO.FastaIO import SimpleFastaParser


def run(assembly_fp):
    # Get contig lengths
    stats = dict()
    with assembly_fp.open('r') as f:
        contig_lengths = [len(s) for d, s in SimpleFastaParser(f)]
    # Calculate stats
    q1, q2, q3 = calculate_quartiles(contig_lengths)
    stats['contig_number'] = len(contig_lengths)
    # If we only have one contig, use that as n50 (otherwise n50 calc fails)
    if len(contig_lengths) == 1:
        stats['n50'] = largest
    else:
        stats['n50'] = calculate_n50(contig_lengths, sum(contig_lengths)/2)
    stats['q1'] = q1
    stats['q2'] = q2
    stats['q3'] = q3
    stats['mean'] = int(round(statistics.mean(contig_lengths), 0))
    stats['smallest'] = min(contig_lengths)
    stats['largest'] = max(contig_lengths)
    stats['length'] = sum(contig_lengths)

    # TODO: apply some thresholds for QC

    return stats


def calculate_quartiles(lengths):
    # Set up
    np = [0.25, 0.50, 0.75]
    x = sorted(lengths)
    n = len(x)
    # Get bounds
    indices = [(n - 1) * p for p in np]
    lo = [math.floor(i) for i in indices]
    hi = [math.ceil(i) for i in indices]
    qs = [x[i] for i in lo]
    # Update if required and then return
    for i in range(len(indices)):
        if not indices[i] > lo[i]:
            continue
        h = indices[i] - lo[i]
        qs[i] = (1 - h) * qs[i] + h * x[hi[i]]
    return qs


def calculate_n50(lengths, median):
    csum = 0
    for length in sorted(lengths, reverse=True):
        # If cumulative sum exceeds median, return last length
        if csum > median:
            return prev_length
        # Cumulative sum
        csum += length
        # Set current to prev_length and then loop
        prev_length = length

