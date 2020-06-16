import copy
import pathlib
import tempfile


import Bio.SeqIO
import cobra.io
import cobra.manipulation
import cobra.manipulation.modify


from . import alignment
from . import util


def run(assembly_fp, ref_gbk_fp, model, output_fp):
    # Get orthologs of genes in model
    model_genes = {gene.id for gene in model.genes}
    isolate_orthologs = identify(assembly_fp, ref_gbk_fp, model_genes)

    # Remove genes from model that have no ortholog in the isolate
    missing_genes = list()
    for gene in model_genes - set(isolate_orthologs):
        # TODO: handle artificial genes better
        if gene == 'KPN_SPONT':
            continue
        missing_genes.append(model.genes.get_by_id(gene))
    # Mutate model inplace and rename genes
    model.id = copy.deepcopy(assembly_fp.stem)
    cobra.manipulation.remove_genes(model, missing_genes, remove_reactions=True)
    cobra.manipulation.modify.rename_genes(model, isolate_orthologs)

    # Write model to disk
    with output_fp.open('w') as fh:
        cobra.io.save_json_model(model, fh)

    # TODO: draft model FBA in basic media for QC
    # TODO: troubleshooting information if FBA fails to optimise to target


def identify(iso_fp, ref_fp, model_genes):
    #pickle_mode = 'read'
    #pickle_mode = 'write'
    pickle_mode = 'read'
    # First we perform a standard best bi-directional hit analysis to identify orthologs
    # Extract protein sequences from both genomes but only keep model genes from the reference
    dh = tempfile.TemporaryDirectory()  # fs directory deleted this is out of scope
    iso_protein_fp = util.write_gbk_sequence(iso_fp, dh.name, seq_type='prot')
    ref_protein_fp = util.write_gbk_sequence(ref_fp, dh.name, model_genes, seq_type='prot')
    # Run BLASTp bidirectionally (filtering with evalue <=1e-3, coverage >=25%, and pident >=80%)
    #blastp_iso = alignment.run_blastp(iso_protein_fp, ref_protein_fp, dh.name)
    #blastp_ref = alignment.run_blastp(ref_protein_fp, iso_protein_fp, dh.name)

    # TEMP: store/load blastp results to/from disk so we dont have to recompute
    import pickle
    if pickle_mode == 'write':
        with open('pickled/blastp_iso.bin', 'wb') as fh:
            pickle.dump(blastp_iso, fh)
        with open('pickled/blastp_ref.bin', 'wb') as fh:
            pickle.dump(blastp_ref, fh)
    elif pickle_mode == 'read':
        with open('pickled/blastp_iso.bin', 'rb') as fh:
            blastp_iso = pickle.load(fh)
        with open('pickled/blastp_ref.bin', 'rb') as fh:
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
    ref_gene_fp = util.write_gbk_sequence(ref_fp, dh.name, model_genes_no_orth, seq_type='nucl')
    iso_fasta_fp = util.write_gbk_to_fasta(iso_fp, dh.name)
    # Run BLASTn (filtering with evalue <=1e-3, coverage >=80%, and pident >=80%)
    blastn_res = alignment.run_blastn(ref_gene_fp, iso_fasta_fp, dh.name)

    # Discover unannotated model genes in isolate
    with iso_fp.open('r') as fh:
        iso_fasta = {record.name: record.seq for record in Bio.SeqIO.parse(iso_fp, 'genbank')}
    for ref_gene_name, hits in blastn_res.items():
        for hit in hits:
            # Get nucleotide sequence and check for premature stop codons
            nucleotide_seq = util.extract_nucleotides_from_ref(hit, iso_fasta)
            if '*' in nucleotide_seq.translate()[:-1]:
                continue
            assert ref_gene_name not in model_orthologs
            model_orthologs[ref_gene_name] = f'{hit.qseqid}_unannotated'
    return model_orthologs
