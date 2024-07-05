import contextlib
import math
import pathlib
import sys
import tempfile
import csv

import Bio.Seq
import Bio.SeqIO
import cobra.core.reaction
import cobra.flux_analysis
import cobra.io
import cobra.manipulation
import cobra.manipulation.modify


from . import alignment
from . import package_data
from . import util


def run(config):
    print('\n========================================')
    print('running draft model creation of ' + config.assembly_genbank_fp.stem)
    print('========================================')
    # Get orthologs of genes in model
    model_genes = {gene.id for gene in config.model.genes}
    isolate_orthologs, blast_results, unannotated_sequences = identify(
        config.assembly_genbank_fp,
        config.model_ref_genes_fp,
        config.model_ref_proteins_fp,
        model_genes,
        config.alignment_thresholds,
    )

    # Writing all identified unannotated sequences to a fasta file (can't be matched to genbank)
    fasta_output = config.output_fp.parent / f'{config.output_fp.stem}_unannotated_sequences.fasta'
    with open(fasta_output, "w") as output_handle:
        Bio.SeqIO.write(unannotated_sequences, output_handle, "fasta")

    # Remove genes from model that have no ortholog in the isolate
    missing_genes = list()
    for gene in model_genes - set(isolate_orthologs):
        # NOTE: must handle artificial genes better
        if gene == 'KPN_SPONT':
            continue
        missing_genes.append(config.model.genes.get_by_id(gene))
    # Mutate a copy of the model and rename genes
    model_draft = config.model.copy()
    model_draft.id = config.assembly_genbank_fp.stem
    cobra.manipulation.remove_genes(model_draft, missing_genes, remove_reactions=True)

    # Save original gene IDs prior to replacing with genome annotations
    original_genes = []
    for gene in model_draft.genes:
        original_genes.append(gene.id)
    model_draft.notes['Original_Genes'] = original_genes

    # Same gene dictionary of reference model and genome annotations to csv 
    gene_dict_fp = config.output_fp.parent / f'{config.output_fp.stem}_gene_dictionary.csv'
    with open(gene_dict_fp, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in isolate_orthologs.items():
            writer.writerow([key, value])
            
    # Mutate a copy of the model and rename genes
    cobra.manipulation.modify.rename_genes(model_draft, isolate_orthologs)

    # Write model to disk and assess model
    with config.model_output_fp.open('w') as fh:
        cobra.io.save_json_model(model_draft, fh)
        cobra.io.write_sbml_model(model_draft, str(config.model_output_fp).rsplit('.', 1)[0] + '.xml')  # .xml output
    assess_model(
        config.model,
        model_draft,
        blast_results,
        config.media_type,
        config.atmosphere_type,
        config.biomass_reaction_id,
        config.model_output_fp,
    )

    # Generate MEMOTE report file if requested
    if config.memote_report_fp:
        util.generate_memote_report(model_draft, config.memote_report_fp)


def assess_model(
    model,
    model_draft,
    blast_results,
    media_type,
    atmosphere_type,
    biomass_reaction_id,
    output_fp,
):
    # pylint: disable=too-many-branches
    # Assess model by observing whether the objective function for biomass optimises
    # We perform an set media
    for reaction in model_draft.exchanges:
        reaction.lower_bound = 0
    media = package_data.get_data('media_definitions', media_type)
    for reaction_id, lower_bound in media['exchanges'].items():
        try:
            reaction = model_draft.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: draft model does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound

    # Set atmospheric conditions
    reaction_oxygen = model_draft.reactions.get_by_id('EX_o2_e')
    if atmosphere_type == 'aerobic':
        reaction_oxygen.lower_bound = -20
    elif atmosphere_type == 'anaerobic':
        reaction_oxygen.lower_bound = 0
    elif atmosphere_type is not None:
        raise ValueError

    solution = model_draft.optimize()

    # Threshold for whether a model produces biomass
    if solution.objective_value < 1e-4:
        msg = (
            f'error: {model_draft} model failed to produce biomass on minimal media, '
            f'construct patch.json manually using {model_draft}'
            f'_model.json.troubleshoot_summary.txt and fix draft {model_draft}'
            f' model via patch_model command'
        )

        print(msg, file=sys.stderr)
        create_troubleshooter(model, model_draft, blast_results, biomass_reaction_id, f'{output_fp}.troubleshoot')
        sys.exit(101)
    else:
        print(f'{model_draft} model produces biomass on minimal media')
        msg = (
            'Please cite:\n'
            '  - Vezina, B., Watts, S.C. et al. '
            'Bactabolize: A tool for high-throughput generation of bacterial strain-specific metabolic models. '
            'https://doi.org/10.1101/2023.02.26.530115 \n'
            '  - Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. '
            'COBRApy: COnstraints-Based Reconstruction and Analysis for Python. '
            'BMC Syst Biol 7, 74 (2013). https://doi.org/10.1186/1752-0509-7-74'
        )
        print(msg)


def create_troubleshooter(model, model_draft, blast_results, biomass_reaction_id, prefix):
    # Determine what required products model cannot product and missing reactions/genes
    reactions_missing, gapfilled, nonzero_threshold = gapfill_model(model, model_draft)
    metabolites_missing = check_biomass_metabolites(model_draft.copy(), biomass_reaction_id)
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
            if len(blast_results['blastn']) != 0:
                blastn_hits[(reaction, gene)] = blast_results['blastn'].get(gene.id, list())
            else:
                blastn_hits[(reaction, gene)] = []
    # Write summary info
    output_fp = pathlib.Path(f'{prefix}_summary.txt')
    write_troubleshoot_summary(
        model,
        metabolites_missing,
        reactions_missing,
        gapfilled,
        nonzero_threshold,
        blastp_hits,
        blastn_hits,
        biomass_reaction_id,
        output_fp,
    )
    # Write BLAST results
    write_blast_results(blastp_hits, pathlib.Path(f'{prefix}_blastp.tsv'))
    write_blast_results(blastn_hits, pathlib.Path(f'{prefix}_blastn.tsv'))


def check_biomass_metabolites(model, biomass_reaction_id):
    metabolites_missing = list()
    reaction_biomass = model.reactions.get_by_id(biomass_reaction_id)
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
    # pylint: disable=undefined-loop-variable
    # Attempt gapfilling with thresholds ranging from default 1e-6 to 0
    integer_thresholds = [math.pow(10, y) for y in (-6, -7, -10, -20, -50, -math.inf)]
    for threshold in integer_thresholds:
        gapfiller = cobra.flux_analysis.gapfilling.GapFiller(
            model=model_draft,
            universal=model,
            demand_reactions=False,
            integer_threshold=threshold,
        )
        try:
            gapfilled = gapfiller.fill(
                iterations=5,
            )
            break
        except RuntimeError:
            continue

    # Gather missing reactions
    reactions_missing = dict()
    for result in gapfilled:
        for reaction in result:
            if reaction not in reactions_missing:
                reactions_missing[reaction] = 0
            reactions_missing[reaction] += 1
    return reactions_missing, gapfilled, threshold


def write_troubleshoot_summary(
    model,
    metabolites_missing,
    reactions_missing,
    gapfilled,
    nonzero_threshold,
    blastp_hits,
    blastn_hits,
    biomass_reaction_id,
    output_fp,
):
    # pylint: disable=cell-var-from-loop,consider-using-with,too-many-branches
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
            if biomass_reaction_id in reactions:
                reactions.remove(biomass_reaction_id)
            print('\nid:', metabolite.id, file=fh)
            print('name:', metabolite.name, file=fh)
            print('reactions:', ', '.join(r.id for r in reactions), file=fh)
            for url_type, url_ids in metabolite.annotation.items():
                if url_type not in reaction_external_urls:
                    continue
                for url_id in url_ids:
                    print('\t', reaction_external_urls[url_type] + url_id, sep='', file=fh)

        msg = f'Missing reactions required to fix model (iterations: 5; threshold: {nonzero_threshold})'
        print(f'\n{msg}', sep='', file=fh)

        for reaction, count in reactions_missing.items():
            print(reaction.id, f'{count}/5', file=fh)
        # Reformat gapfilled information in troubleshoot summary to include all possible combinations
        count = 1
        for iteration in gapfilled:
            print('\nIteration ', count, sep='', file=fh)
            for reaction in iteration:
                print(reaction.id, file=fh)
            count += 1
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
    # pylint: disable=consider-using-with
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
            print(*[seq[i : i + 80] for i in range(0, len(seq), 80)], sep='\n', file=fout)
    # Run BLASTn (filtering with evalue <=1e-3, coverage >=80%, and pident >=80%)
    blastn_res_all = alignment.run_blastn(ref_genes_noorth_fp, iso_fasta_fp)
    if len(blastn_res_all) != 0:
        blastn_res = alignment.filter_results(blastn_res_all, min_coverage=80, min_pident=80)
    else:
        blastn_res = {}
    # Discover unannotated model genes in isolate
    model_orthologs, unannotated_sequences = discover_unannotated_orthologs(blastn_res, iso_fasta_fp, model_orthologs)

    # Return orthologs and BLAST results, explicitly remove temp directory
    dh.cleanup()
    blast_results = {'blastp_iso': blastp_iso_all, 'blastp_ref': blastp_ref_all, 'blastn': blastn_res_all}
    return model_orthologs, blast_results, unannotated_sequences


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
    # pylint: disable=no-else-continue
    with iso_fasta_fp.open('r') as fh:
        iso_fasta = {des: Bio.Seq.Seq(seq) for des, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
    unannotated_sequences = []
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
                # Record information for unannotated sequences
                record = Bio.SeqIO.SeqRecord(nucleotide_seq, id=f'{hit.qseqid}_unannotated', description='')
                unannotated_sequences.append(record)
                break
    # Returning mutable to be explicit
    return model_orthologs, unannotated_sequences
