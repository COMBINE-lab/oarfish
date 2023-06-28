process map_samples {
  debug true
  cpus 9

  input:
    tuple val(sample_id), path(sample), val(ref_type), path(ref)
  output:
    tuple val(sample_id), val(ref_type), path("${sample_id}_${ref_type}.bam")

    """
    minimap2 -t 8 -ax map-ont -N 100 ${ref} ${sample} | samtools view -@4 -h -F 2052 -bS > ${sample_id}_${ref_type}.bam
    """
}

process build_sirv_index {
  input:
    tuple path(ref_fasta), path(ref_gtf), val(ref_type)
  output:
    tuple val(ref_type), path("sirv_isoforms_${ref_type}.fa")

  """
  gffread -E -v -w sirv_isoforms_${ref_type}.fa -g ${ref_fasta} ${ref_gtf}
  """
}

process quant_samples_oarfish {
  publishDir "sirv/quants/oarfish"

  input:
    tuple val(sample_id), val(ref_type), path(ref_bam)
  output:
    path "${sample_id}_${ref_type}_quant.tsv"

  """
  ${params.oarfish} -a ${ref_bam} -t 50 -o ${sample_id}_${ref_type}_quant.tsv
  """
}

process quant_samples_nanocount {
  publishDir "sirv/quants/nanocount"

  input:
    tuple val(sample_id), val(ref_type), path(ref_bam)
  output:
    path "${sample_id}_${ref_type}_quant.tsv"

  """
  NanoCount -i ${ref_bam} -x -o ${sample_id}_${ref_type}_quant.tsv
  """
}

workflow {
  def sirv_refs = Channel.fromList(params.sirv_ref_inputs)
  build_sirv_index(sirv_refs)

  def sirv_sample_paths = Channel.fromList(params.sirv_samples)
    .map{ it -> [it, "${params.sirv_input_path}/${it}_1.fastq.gz"] }

  def sirv_map_input = sirv_sample_paths
    .combine(build_sirv_index.out)

  // map the samples with minimap2
  map_samples(sirv_map_input)

  // quantify with both oarfish and NanoCount
  quant_samples_nanocount(map_samples.out)
  quant_samples_oarfish(map_samples.out)
}
