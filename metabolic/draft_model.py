import pathlib
import sys
import tempfile


import Bio.SeqIO
import cobra.flux_analysis
import cobra.io
import cobra.manipulation
import cobra.manipulation.modify


from . import alignment
from . import media_definitions
from . import util


def run(assembly_fp, ref_genbank_fp, model, output_fp):
    print('\n========================================')
    print('running draft model creation')
    print('========================================')
    # Get orthologs of genes in model
    model_genes = {gene.id for gene in model.genes}
    isolate_orthologs, blast_results = identify(assembly_fp, ref_genbank_fp, model_genes)

    # Remove genes from model that have no ortholog in the isolate
    missing_genes = list()
    for gene in model_genes - set(isolate_orthologs):
        # TODO: handle artificial genes better
        if gene == 'KPN_SPONT':
            continue
        missing_genes.append(model.genes.get_by_id(gene))
    # Mutate a copy of the model and rename genes
    model_draft = model.copy()
    model_draft.id = assembly_fp.stem
    cobra.manipulation.remove_genes(model_draft, missing_genes, remove_reactions=True)
    cobra.manipulation.modify.rename_genes(model_draft, isolate_orthologs)

    # Write model to disk and assess model
    with output_fp.open('w') as fh:
        cobra.io.save_json_model(model_draft, fh)
    assess_model(model, model_draft, blast_results, output_fp)


def assess_model(model, model_draft, blast_results, output_fp):
    # Assess model by observing whether the objective function for biomass optimises
    # We perform an FBA on minimal media (m9)
    for reaction in model_draft.exchanges:
        reaction.lower_bound = 0
    for reaction_id, lower_bound in media_definitions.m9.items():
        try:
            reaction = model_draft.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: draft model does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound
    solution = model_draft.optimize()

    # Somewhat arbitrary threshold for whether a model produces biomass
    if solution.objective_value < 1e-4:
        msg = (f'error: model failed to produce biomass on minimal media, '
                'manual intervention is required to fix the draft model')
        print(msg, file=sys.stderr)
        create_troubleshooter(model, model_draft, blast_results, f'{output_fp}.troubleshoot')
        sys.exit(1)
    else:
        print('model produces biomass on minimal media')


def create_troubleshooter(model, model_draft, blast_results, prefix):
    # Determine which genes are missing
    gapfilled = cobra.flux_analysis.gapfill(model_draft, model, demand_reactions=False, iterations=5)
    reaction_counts = dict()
    for result in gapfilled:
        for reaction in result:
            if reaction not in reaction_counts:
                reaction_counts[reaction] = 0
            reaction_counts[reaction] += 1
    # Collect BLAST results
    blastp_hits = dict()
    blastn_hits = dict()
    for reaction in reaction_counts:
        for gene in reaction.genes:
            assert (reaction, gene) not in blastp_hits
            hitsp = blast_results['blastp_ref'].get(gene.id, list())
            blastp_hits[(reaction, gene)] = hitsp
            if len(hitsp) > 0:
                continue
            blastn_hits[(reaction, gene)] = blast_results['blastn'].get(gene.id, list())
    # Set base URLs for reaction annotations
    reaction_external_urls = {
        'bigg.reaction': 'http://bigg.ucsd.edu/universal/reactions/',
        'biocyc': 'http://identifiers.org/biocyc/',
        'ec-code': 'http://identifiers.org/ec-code/',
        'kegg.reaction': 'http://identifiers.org/kegg.reaction/',
        'metanetx.reaction': 'http://identifiers.org/metanetx.reaction/',
        'rhea': 'http://identifiers.org/rhea/',
        'seed.reaction': 'http://identifiers.org/seed.reaction/',
    }
    # Write summary info
    output_fp = pathlib.Path(f'{prefix}_summary.txt')
    with output_fp.open('w') as fh:
        print('Gapfilling results (5 iterations)', sep='', file=fh)
        for reaction, count in reaction_counts.items():
            print(reaction.id, f'{count}/5', file=fh)

        print('\n', 'BLASTp hits', sep='', file=fh)
        for (reaction, gene), hits in blastp_hits.items():
            print(reaction.id, gene.id, len(hits), file=fh)

        print('\n', 'BLASTn hits (only done for ORFs with no BLASTp result)', sep='', file=fh)
        for (reaction, gene), hits in blastn_hits.items():
            print(reaction.id, gene.id, len(hits), file=fh)

        print('\n', 'BiGG info', sep='', end='', file=fh)
        for reaction in reaction_counts:
            print('\nname: ', reaction.id, sep='', file=fh)
            print('reaction:', reaction.reaction, file=fh)
            print('urls:', file=fh)
            for url_type, url_ids in reaction.annotation.items():
                if url_type not in reaction_external_urls:
                    continue
                for url_id in url_ids:
                    print('\t', reaction_external_urls[url_type] + url_id, sep='', file=fh)
    # Write BLAST results
    write_blast_results(blastp_hits, pathlib.Path(f'{prefix}_blastp.tsv'))
    write_blast_results(blastn_hits, pathlib.Path(f'{prefix}_blastn.tsv'))


def write_blast_results(data, output_fp):
    with output_fp.open('w') as fh:
        print(*alignment.BlastFormat, sep='\t', file=fh)
        for hits in data.values():
            print(*hits, sep='\n', file=fh)


def identify(iso_fp, ref_fp, model_genes):
    pickle_mode = 'read'
    #pickle_mode = 'write'
    #pickle_mode = 'noop'
    # First we perform a standard best bi-directional hit analysis to identify orthologs
    # Extract protein sequences from both genomes but only keep model genes from the reference
    dh = tempfile.TemporaryDirectory()
    #iso_protein_fp = util.write_genbank_coding_sequence(iso_fp, dh.name, seq_type='prot')
    #ref_protein_fp = util.write_genbank_coding_sequence(ref_fp, dh.name, model_genes, seq_type='prot')
    ## Run BLASTp bidirectionally (filtering with evalue <=1e-3, coverage >=25%, and pident >=80%)
    #blastp_iso_all = alignment.run_blastp(iso_protein_fp, ref_protein_fp, dh.name)
    #blastp_ref_all = alignment.run_blastp(ref_protein_fp, iso_protein_fp, dh.name)
    #blastp_iso = alignment.filter_results(blastp_iso_all, min_coverage=25, min_pident=80)
    #blastp_ref = alignment.filter_results(blastp_ref_all, min_coverage=25, min_pident=80)

    # TEMP: store/load blastp results to/from disk so we dont have to recompute
    import pickle
    if pickle_mode == 'write':
        with open(f'pickled/{iso_fp.stem}_blastp_iso_all.bin', 'wb') as fh:
            pickle.dump(blastp_iso_all, fh)
        with open(f'pickled/{iso_fp.stem}_blastp_ref_all.bin', 'wb') as fh:
            pickle.dump(blastp_ref_all, fh)
        with open(f'pickled/{iso_fp.stem}_blastp_iso.bin', 'wb') as fh:
            pickle.dump(blastp_iso, fh)
        with open(f'pickled/{iso_fp.stem}_blastp_ref.bin', 'wb') as fh:
            pickle.dump(blastp_ref, fh)
    elif pickle_mode == 'read':
        with open(f'pickled/{iso_fp.stem}_blastp_iso_all.bin', 'rb') as fh:
            blastp_iso_all = pickle.load(fh)
        with open(f'pickled/{iso_fp.stem}_blastp_ref_all.bin', 'rb') as fh:
            blastp_ref_all = pickle.load(fh)
        with open(f'pickled/{iso_fp.stem}_blastp_iso.bin', 'rb') as fh:
            blastp_iso = pickle.load(fh)
        with open(f'pickled/{iso_fp.stem}_blastp_ref.bin', 'rb') as fh:
            blastp_ref = pickle.load(fh)
    elif pickle_mode == 'noop':
        pass
    else:
        assert False

    # Find orthologs from BLASTp results
    model_orthologs = dict()
    for ref_gene_name, hits in blastp_ref.items():
        # Get the best hit in isolate for this reference protein
        best_iso_hit = max(hits, key=lambda k: k.pident)
        # For the isolate protein collect the best hit in the reference, if any
        if best_iso_hit.sseqid not in blastp_iso:
            continue
        best_ref_hit = max(blastp_iso[best_iso_hit.sseqid], key=lambda k: k.pident)
        # If they are reciprocal consider them orthologous
        if best_ref_hit.sseqid == ref_gene_name:
            model_orthologs[ref_gene_name] = best_ref_hit.qseqid

    # Next we check for unannotated genes by comparing reference genes to isolate genome
    # Get reference genes that do not have an ortholog
    model_genes_no_orth = model_genes.difference(set(model_orthologs))
    # Extract gene sequences from reference and write isolate sequence as fasta
    ref_gene_fp = util.write_genbank_coding_sequence(ref_fp, dh.name, model_genes_no_orth, seq_type='nucl')
    iso_fasta_fp = util.write_genbank_to_fasta(iso_fp, dh.name)
    # Run BLASTn (filtering with evalue <=1e-3, coverage >=80%, and pident >=80%)
    blastn_res_all = alignment.run_blastn(ref_gene_fp, iso_fasta_fp, dh.name)
    blastn_res = alignment.filter_results(blastn_res_all, min_coverage=80, min_pident=80)

    # Explicitly remove temporary directory
    dh.cleanup()

    # Discover unannotated model genes in isolate
    with iso_fp.open('r') as fh:
        iso_fasta = {record.name: record.seq for record in Bio.SeqIO.parse(iso_fp, 'genbank')}
    for ref_gene_name, hits in blastn_res.items():
        assert ref_gene_name not in model_orthologs
        for hit in hits:
            # Get nucleotide sequence and check for premature stop codons
            nucleotide_seq = util.extract_nucleotides_from_ref(hit, iso_fasta)
            # Append trailing N if we have a partial codon
            seq_n = (3 - len(nucleotide_seq)) % 3
            nucleotide_seq = nucleotide_seq + 'N' * seq_n
            if '*' in nucleotide_seq.translate()[:-1]:
                continue
            else:
                model_orthologs[ref_gene_name] = f'{hit.qseqid}_unannotated'
                break

    # Return orthologs and BLAST results
    blast_results = {'blastp_iso': blastp_iso_all, 'blastp_ref': blastp_ref_all, 'blastn': blastn_res_all}
    return model_orthologs, blast_results
