class ConfigDraftModel:

    def __init__(self, args):
        self.ref_model_fp = args.ref_model_fp
        self.assembly_fp = args.assembly_fp
        self.ref_genbank_fp = args.ref_genbank_fp
        self.ref_genes_fp = args.ref_genes_fp
        self.ref_proteins_fp = args.ref_proteins_fp
        self.media_type = args.media_type
        self.min_coverage = args.min_coverage
        self.min_pident = args.min_pident
        self.min_ppos = args.min_ppos
        self.no_reannotation = args.no_reannotation
        self.memote_report_fp = args.memote_report_fp
        self.output_fp = args.output_fp

        self.alignment_thresholds = None
        self.assembly_genbank_fp = None
        self.model = None
        self.model_genes_fp  = None
        self.model_proteins_fp = None
        self.model_output_fp = None


class ConfigPatchModel:

    def __init__(self, args):

        self.draft_model_fp = args.draft_model_fp
        self.ref_model_fp = args.ref_model_fp
        self.patch_fp = args.patch_fp
        self.media_type = args.media_type
        self.output_fp = args.output_fp
        self.memote_report_fp = args.memote_report_fp


class ConfigFba:

    def __init__(self, args):

        self.model_fp = args.model_fp
        self.fba_open_value = args.fba_open_value
        self.fba_spec_fp = args.fba_spec_fp
        self.output_fp = args.output_fp
