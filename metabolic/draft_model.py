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

    # Write model to disk
    with output_fp.open('w') as fh:
        cobra.io.save_json_model(model_draft, fh)

    # Assess model by observing whether the objective function for biomass optimises
    # We perform an FBA on minimal media (m9)
    for reaction in model_draft.reactions:
        if reaction.id.startswith('EX_'):
            reaction.lower_bound = 0
    for reaction_id, lower_bound in media_definitions.m9.items():
        try:
            reaction = model_draft.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: draft model {assembly_fp.stem} does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound
    solution = model_draft.optimize()

    # Somewhat arbitrary threshold for whether a model produces biomass
    if solution.objective_value < 1e-4:
        msg = (f'error: model for {assembly_fp.stem} failed to produce biomass on minimal media, '
                'manual intervention is required to fix the draft model')
        print(msg, file=sys.stderr)

        # TODO: place this troubleshooting information in an output
        gapfilled = cobra.flux_analysis.gapfill(model_draft, model, demand_reactions=False, iterations=5)
        reaction_counts = dict()
        for result in gapfilled:
            for reaction in result:
                if reaction not in reaction_counts:
                    reaction_counts[reaction] = 0
                reaction_counts[reaction] += 1
        print('\n', 'Gapfilling results', sep='')
        for reaction, count in reaction_counts.items():
            print(reaction.id, count)

        print('\n', 'BLASTp hits', sep='')
        for reaction in reaction_counts:
            for gene in reaction.genes:
                hits = blast_results['blastp_ref'].get(gene.id, list())
                print(reaction.id, gene.id, len(hits))

        print('\n', 'BLASTn hits', sep='')
        for reaction in reaction_counts:
            for gene in reaction.genes:
                hits = blast_results['blastn'].get(gene.id, list())
                print(reaction.id, gene.id, len(hits))

        with open('data/bigg_models_reactions.txt', 'r') as fh:
            line_token_gen = (line.rstrip().split('\t') for line in fh)
            header_tokens = next(line_token_gen)
            records = dict()
            for line_tokens in line_token_gen:
                record = {k: v for k, v in zip(header_tokens, line_tokens)}
                assert record['bigg_id'] not in records
                records[record['bigg_id']] = record

        print('\n', 'BiGG info', sep='')
        for reaction in reaction_counts:
            record = records[reaction.id]
            print('name: ', reaction.id, sep='')
            print('reaction:', record['reaction_string'])
            print('urls:')
            for url_entry in record['database_links'].split('; '):
                print('\t', url_entry, sep='')
            print()

        # Other information about these genes
        #       - is the gene ORF complete
        #           - see genbank notes if reannotated
        #           - check if at bounds

        sys.exit(1)
    else:
        # TODO: report in more meaningful way
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            print()
            print(model.summary())


def identify(iso_fp, ref_fp, model_genes):
    #pickle_mode = 'read'
    #pickle_mode = 'write'
    pickle_mode = 'noop'
    # First we perform a standard best bi-directional hit analysis to identify orthologs
    # Extract protein sequences from both genomes but only keep model genes from the reference
    dh = tempfile.TemporaryDirectory()
    iso_protein_fp = util.write_genbank_coding_sequence(iso_fp, dh.name, seq_type='prot')
    ref_protein_fp = util.write_genbank_coding_sequence(ref_fp, dh.name, model_genes, seq_type='prot')
    # Run BLASTp bidirectionally (filtering with evalue <=1e-3, coverage >=25%, and pident >=80%)
    blastp_iso = alignment.run_blastp(iso_protein_fp, ref_protein_fp, dh.name)
    blastp_ref = alignment.run_blastp(ref_protein_fp, iso_protein_fp, dh.name)

    # TEMP: store/load blastp results to/from disk so we dont have to recompute
    import pickle
    if pickle_mode == 'write':
        with open(f'pickled/{iso_fp.stem}_blastp_iso.bin', 'wb') as fh:
            pickle.dump(blastp_iso, fh)
        with open(f'pickled/{iso_fp.stem}_blastp_ref.bin', 'wb') as fh:
            pickle.dump(blastp_ref, fh)
    elif pickle_mode == 'read':
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
    blastn_res = alignment.run_blastn(ref_gene_fp, iso_fasta_fp, dh.name)

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
    blast_results = {'blastp_iso': blastp_iso, 'blastp_ref': blastp_ref, 'blastn': blastn_res}
    return model_orthologs, blast_results
