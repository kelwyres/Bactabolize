process assembly_qc {
  publishDir "${params.output_dir}/${isolate_id}/"

  input:
  tuple isolate_id, path(assembly_fp)

  output:
  path('*_stats.tsv')

  script:
  """
  metabolic assembly_qc --assembly_fp ${assembly_fp} --output_fp ${isolate_id}_stats.tsv
  """
}

process annotate {
  publishDir "${params.output_dir}/${isolate_id}/"

  input:
  tuple isolate_id, path(assembly_fp)
  path(prodigal_model_fp)

  output:
  tuple val(isolate_id), path('*_reannotated.gbk')

  script:
  """
  metabolic annotate --assembly_fp ${assembly_fp} --prodigal_model_fp ${prodigal_model_fp} --output_fp ${isolate_id}_reannotated.gbk
  """
}

process draft_model {
  publishDir "${params.output_dir}/${isolate_id}/"

  input:
  tuple isolate_id, path(assembly_fp)
  path(ref_genbank_fp)
  path(ref_model_fp)

  output:
  tuple val(isolate_id), path('*_model.json') optional true

  script:
  """
  metabolic draft_model --assembly_fp ${assembly_fp} --ref_genbank_fp ${ref_genbank_fp} --ref_model_fp ${ref_model_fp} --output_fp ${isolate_id}_model.json
  """
}

process model_fba {
  publishDir "${params.output_dir}/${isolate_id}/"

  input:
  tuple isolate_id, path(model_fp)

  output:
  path('*_fba.txt')

  script:
  """
  metabolic model_fba --model_fp ${model_fp} --fba_spec_fp /dev/null --output_fp /dev/null > ${isolate_id}_fba.txt
  """
}
