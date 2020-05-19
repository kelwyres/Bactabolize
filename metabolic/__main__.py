import cobra.io
import cobra.manipulation
import cobra.manipulation.modify


from . import arguments
from . import orthologs


def entry():
    # Get command line arguments
    args = arguments.get_args()

    # Parse model and get list of genes
    with args.ref_model_fp.open('r') as fh:
        model = cobra.io.load_json_model(fh)
    model_genes = {gene.id for gene in model.genes}

    # Get orthologs of genes in model
    isolate_orthologs = orthologs.identify(args.ref_gbk_fp, args.isolate_fp, model_genes)

    # Remove genes from model that have no ortholog in the isolate
    missing_genes = list()
    for gene in model_genes - set(isolate_orthologs):
        # TODO: handle artifical genes better
        if gene == 'KPN_SPONT':
            continue
        missing_genes.append(model.genes.get_by_id(gene))
    # Mutate model inplace and rename genes
    # NOTE: will need deep copy if we're to process more than one isolate
    model.id = args.isolate_fp.stem
    cobra.manipulation.remove_genes(model, missing_genes, remove_reactions=True)
    cobra.manipulation.modify.rename_genes(model, isolate_orthologs)

    # Write model to disk
    with args.output_fp.open('w') as fh:
        cobra.io.save_json_model(model, fh)


if __name__ == '__main__':
    entry()
