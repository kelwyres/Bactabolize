process assembly_qc {
  publishDir "${params.output_dir}/${isolate_id}/"

  input:
  tuple isolate_id, path(assembly_fp)

  output:
  path('*_stats.tsv')

  script:
  """
  bactabolize assembly_qc --assembly_fp ${assembly_fp} --output_fp ${isolate_id}_stats.tsv
  """
}

process annotate {
  publishDir "${params.output_dir}/${isolate_id}/", saveAs: { "${isolate_id}.gbk" }

  input:
  tuple isolate_id, path(assembly_fp)

  output:
  tuple val(isolate_id), path('*_annotated.gbk')

  script:
  """
  bactabolize annotate --assembly_fp ${assembly_fp} --output_fp ${isolate_id}_annotated.gbk
  """
}

process draft_model {
  publishDir "${params.output_dir}/${isolate_id}/"
  validExitStatus 0, 101

  input:
  tuple isolate_id, path(assembly_fp)
  path(ref_genbank_fp)
  path(ref_model_fp)

  output:
  tuple val(isolate_id), path("${isolate_id}_model*")

  script:
  """
  bactabolize draft_model --assembly_fp ${assembly_fp} --ref_genbank_fp ${ref_genbank_fp} --ref_model_fp ${ref_model_fp} --output_fp ${isolate_id}_model.json
  """
}

process model_fba {
  publishDir "${params.output_dir}/${isolate_id}/"

  input:
  tuple isolate_id, path(model_fp)

  output:
  path('*_fba.tsv')

  script:
  fba_spec_opt = (params.fba_spec_fp) ? "--fba_spec_fp ${params.fba_spec_fp}" : ''
  fba_types_opt = (params.fba_types) ? "--fba_types ${params.fba_types}" : ''
  """
  bactabolize model_fba --model_fp ${model_fp} ${fba_spec_opt} ${fba_types_opt} --output_fp ${isolate_id}_fba.tsv
  """
}
