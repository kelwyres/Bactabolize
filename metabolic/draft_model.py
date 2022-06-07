import contextlib
import pathlib
import sys
import tempfile


import Bio.Seq
import Bio.SeqIO
import cobra.core.reaction
import cobra.flux_analysis
import cobra.io
import cobra.manipulation
import cobra.manipulation.modify


from . import alignment
from . import media_definitions
from . import util


def run(assembly_fp, ref_genes_fp, ref_proteins_fp, model, alignment_thresholds, output_fp):
    print('\n========================================')
    print('running draft model creation')
    print('========================================')
    # Get orthologs of genes in model
    model_genes = {gene.id for gene in model.genes}
    isolate_orthologs, blast_results = identify(
        assembly_fp,
        ref_genes_fp,
        ref_proteins_fp,
        model_genes,
        alignment_thresholds
    )
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

    # Threshold for whether a model produces biomass
    if solution.objective_value < 1e-4:
        msg = ('error: model failed to produce biomass on minimal media, '
                'manual intervention is required to fix the draft model')
        print(msg, file=sys.stderr)
        create_troubleshooter(model, model_draft, blast_results, f'{output_fp}.troubleshoot')
        sys.exit(101)
    else:
        print('model produces biomass on minimal media')


def create_troubleshooter(model, model_draft, blast_results, prefix):
    # Determine what required products model cannot product and missing reactions/genes
    reactions_missing = gapfill_model(model, model_draft)
    metabolites_missing = check_biomass_metabolites(model_draft.copy())
    # Collect BLAST results
    blastp_hits = dict()
    blastn_hits = dict()
    for reaction in reactions_missing:
        for gene in reaction.genes:
            assert (reaction, gene) not in blastp_hits
            hitsp = blast_results['blastp_ref'].get(gene.id, list())
            blastp_hits[(reaction, gene)] = hitsp
            if len(hitsp) > 0:
                continue
            blastn_hits[(reaction, gene)] = blast_results['blastn'].get(gene.id, list())
    # Write summary info
    output_fp = pathlib.Path(f'{prefix}_summary.txt')
    write_troubleshoot_summary(model, metabolites_missing, reactions_missing, blastp_hits, blastn_hits, output_fp)
    # Write BLAST results
    write_blast_results(blastp_hits, pathlib.Path(f'{prefix}_blastp.tsv'))
    write_blast_results(blastn_hits, pathlib.Path(f'{prefix}_blastn.tsv'))


def check_biomass_metabolites(model):
    metabolites_missing = list()
    reaction_biomass = model.reactions.get_by_id('BIOMASS_')
    for metabolite in reaction_biomass.metabolites:
        # Create a drain reaction for the metabolite, add to model, and set as objective
        reaction = cobra.core.reaction.Reaction(
            id=f'{metabolite.id}_drain',
            subsystem='Transport/Exchange',
            lower_bound=0,
            upper_bound=1000,
        )
        reaction.add_metabolites({metabolite: -1})
        model.add_reaction(reaction)
        model.objective = reaction
        # Optimise model to new objective and check for biomass production
        solution = model.optimize()
        if solution.objective_value < 0.0001:
            metabolites_missing.append(metabolite)
        # Remove reaction to reset model
        model.remove_reactions([reaction])
    return metabolites_missing


def gapfill_model(model, model_draft):
    gapfilled = cobra.flux_analysis.gapfill(model_draft, model, demand_reactions=False, iterations=5)
    reactions_missing = dict()
    for result in gapfilled:
        for reaction in result:
            if reaction not in reactions_missing:
                reactions_missing[reaction] = 0
            reactions_missing[reaction] += 1
    return reactions_missing


def write_troubleshoot_summary(model, metabolites_missing, reactions_missing, blastp_hits, blastn_hits, output_fp):
    # Set base URLs for reaction annotations
    reaction_external_urls = {
        'bigg.metabolite': 'http://bigg.ucsd.edu/universal/metabolites/',
        'bigg.reaction': 'http://bigg.ucsd.edu/universal/reactions/',
        'biocyc': 'http://identifiers.org/biocyc/',
        'ec-code': 'http://identifiers.org/ec-code/',
        'kegg.reaction': 'http://identifiers.org/kegg.reaction/',
        'metanetx.reaction': 'http://identifiers.org/metanetx.reaction/',
        'rhea': 'http://identifiers.org/rhea/',
        'seed.reaction': 'http://identifiers.org/seed.reaction/',
    }
    with output_fp.open('w') as fh:
        print('Missing metabolites required for biomass production:', end='', file=fh)
        for metabolite in metabolites_missing:
            reactions = model.reactions.query(lambda r: metabolite.id in {m.id for m in r.metabolites})
            if 'BIOMASS_' in reactions:
                reactions.remove('BIOMASS_')
            print('\nid:', metabolite.id, file=fh)
            print('name:', metabolite.name, file=fh)
            print('reactions:', ', '.join(r.id for r in reactions), file=fh)
            for url_type, url_ids in metabolite.annotation.items():
                if url_type not in reaction_external_urls:
                    continue
                for url_id in url_ids:
                    print('\t', reaction_external_urls[url_type] + url_id, sep='', file=fh)
        print('\nMissing reactions required to fix model (5 gapfill iterations)', sep='', file=fh)
        for reaction, count in reactions_missing.items():
            print(reaction.id, f'{count}/5', file=fh)
        print('\n', 'Reaction info', sep='', end='', file=fh)
        for reaction in reactions_missing:
            print('\nid: ', reaction.id, sep='', file=fh)
            print('name: ', reaction.name, sep='', file=fh)
            print('subsystem: ', reaction.subsystem, sep='', file=fh)
            print('reaction:', reaction.reaction, file=fh)
            print('genes:', ', '.join(gene.id for gene in reaction.genes), file=fh)
            print('urls:', file=fh)
            for url_type, url_ids in reaction.annotation.items():
                if url_type not in reaction_external_urls:
                    continue
                for url_id in url_ids:
                    print('\t', reaction_external_urls[url_type] + url_id, sep='', file=fh)
        print('\nBLASTp hits', file=fh)
        for (reaction, gene), hits in blastp_hits.items():
            print(reaction.id, gene.id, len(hits), file=fh)
        print('\nBLASTn hits (only done for ORFs with no BLASTp result)', file=fh)
        for (reaction, gene), hits in blastn_hits.items():
            print(reaction.id, gene.id, len(hits), file=fh)


def write_blast_results(data, output_fp):
    with output_fp.open('w') as fh:
        print(*alignment.BlastFormat, sep='\t', file=fh)
        for hits in data.values():
            print(*hits, sep='\n', file=fh)


def identify(iso_fp, ref_genes_fp, ref_proteins_fp, model_genes, alignment_thresholds):
    # First we perform a standard best bi-directional hit analysis to identify orthologs
    # Extract protein sequences from both genomes but only keep model genes from the reference
    dh = tempfile.TemporaryDirectory()
    iso_proteins_fp = pathlib.Path(dh.name, 'isolate_proteins.fasta')
    util.write_genbank_coding(iso_fp, iso_proteins_fp, seq_type='prot')
    # Run BLASTp bidirectionally (filtering with evalue <=1e-3, and user defined coverage, pident, ppos)
    blastp_iso_all = alignment.run_blastp(iso_proteins_fp, ref_proteins_fp)
    blastp_ref_all = alignment.run_blastp(ref_proteins_fp, iso_proteins_fp)
    blastp_iso = alignment.filter_results(blastp_iso_all, **alignment_thresholds)
    blastp_ref = alignment.filter_results(blastp_ref_all, **alignment_thresholds)

    # Find orthologs from BLASTp results
    model_orthologs = discover_orthologs(blastp_ref, blastp_iso)

    # For reference genes without orthologs, we check for unannotated hits
    # Write isolate sequence as fasta
    model_genes_no_orth = model_genes.difference(set(model_orthologs))
    iso_fasta_fp = pathlib.Path(dh.name, 'isolate_genes.fasta')
    util.write_genbank_coding(iso_fp, iso_fasta_fp, seq_type='nucl')
    # Write reference gene sequences with no ortholog as fasta
    ref_genes_noorth_fp = pathlib.Path(dh.name, 'ref_genes_noorth.fasta')
    with contextlib.ExitStack() as stack:
        fin = stack.enter_context(ref_genes_fp.open('r'))
        fout = stack.enter_context(ref_genes_noorth_fp.open('w'))
        for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fin):
            if desc not in model_genes_no_orth:
                continue
            print(f'>{desc}', file=fout)
            print(*[seq[i:i+80] for i in range(0, len(seq), 80)], sep='\n', file=fout)
    # Run BLASTn (filtering with evalue <=1e-3, coverage >=80%, and pident >=80%)
    blastn_res_all = alignment.run_blastn(ref_genes_noorth_fp, iso_fasta_fp)
    blastn_res = alignment.filter_results(blastn_res_all, min_coverage=80, min_pident=80)
    # Discover unannotated model genes in isolate
    model_orthologs = discover_unannotated_orthologs(blastn_res, iso_fasta_fp, model_orthologs)

    # Return orthologs and BLAST results, explicitly remove temp directory
    dh.cleanup()
    blast_results = {'blastp_iso': blastp_iso_all, 'blastp_ref': blastp_ref_all, 'blastn': blastn_res_all}
    return model_orthologs, blast_results


def discover_orthologs(blastp_ref, blastp_iso):
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
    return model_orthologs


def discover_unannotated_orthologs(blastn_res, iso_fasta_fp, model_orthologs):
    with iso_fasta_fp.open('r') as fh:
        iso_fasta = {des: Bio.Seq.Seq(seq) for des, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
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
    # Returning mutable to be explicit
    return model_orthologs
