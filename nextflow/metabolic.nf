#!/usr/bin/env nextflow
// Enable DSL2
nextflow.preview.dsl=2

// Import processes
include assembly_qc from './processes.nf'
include annotate from './processes.nf'
include draft_model from './processes.nf'
include model_fba from './processes.nf'



// Print splash
log.info('--------------------------------------------------------------------')
log.info("""
‌‌

                      dP            dP                dP oo
                      88            88                88
88d8b.d8b. .d8888b. d8888P .d8888b. 88d888b. .d8888b. 88 dP .d8888b.
88'`88'`88 88ooood8   88   88'  `88 88'  `88 88'  `88 88 88 88'  `""
88  88  88 88.  ...   88   88.  .88 88.  .88 88.  .88 88 88 88.  ...
dP  dP  dP `88888P'   dP   `88888P8 88Y8888' `88888P' dP dP `88888P'

‌‌
""".stripIndent())
log.info('--------------------------------------------------------------------')


// Require some variables to be boolean
// We must check and change values if needed. The global param variables are immutable so instead we declare new ones
def check_boolean_option(option, name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "ERROR: ${name} option must be true or false"
}
run_assembly_qc = check_boolean_option(params.assembly_qc, 'assembly_qc')
run_annotation = check_boolean_option(params.reannotation, 'reannotation')
run_fba = check_boolean_option(params.model_fba, 'model_fba')

// Create file input objects and check they exist
def check_file_exists(filepath, name) {
  if (! filepath.exists()) {
    exit 1, "ERROR: reference input '${filepath}' for '${name}' does not exist"
  }
}
ref_genbank_fp = file(params.ref_genbank_fp)
ref_model_fp = file(params.ref_model_fp)
prodigal_model_fp = file(params.prodigal_model_fp)
check_file_exists(ref_genbank_fp, 'ref_genbank_fp')
check_file_exists(ref_model_fp, 'ref_model_fp')
check_file_exists(prodigal_model_fp, 'prodigal_model_fp')

// Check output directory is empty
output_dir_files = []
output_dir = file(params.output_dir)
output_dir.eachFile { output_dir_files.add(it.name) }
run_info_dirname = file(params.run_info_dir).simpleName
output_dir_files.remove(run_info_dirname)
if (output_dir_files.size() > 0 && ! params.force) {
  exit 1, "ERROR: output directory '${output_dir}' already exists and contains other files, remove or use --force to overwrite"
}

// Create input channel for assemblies to emit tuples of (isolate_id, assembly_fp)
assembly_ch = Channel.fromPath(params.assembly_fps).ifEmpty {
    exit 1, "ERROR: did not find any assembly files with '${params.assembly_fps}'"
  }.map { filepath ->
    [filepath.simpleName, filepath]
}

log.info('--------------------------------------------------------------------')

workflow {
  main:
    if (run_assembly_qc) {
      assembly_qc(assembly_ch)
    }
    if (run_annotation) {
      assembly_ch = annotate(assembly_ch, prodigal_model_fp)
    }
    model_draft_ch = draft_model(assembly_ch, ref_genbank_fp, ref_model_fp)
    if (run_fba) {
      model_fba(model_draft_ch)
    }
}
