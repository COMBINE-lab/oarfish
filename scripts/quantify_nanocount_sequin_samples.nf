process map_samples {
  debug true
  cpus 9

  input:
    tuple val(sample_name), path(sample)
  output:
    tuple val(sample_name), path("${sample_name}.bam")

    """
    minimap2 -t 8 -ax map-ont -N 100 ${params.gencode_31_with_sequin} ${sample} | samtools view -@4 -h -F 2052 -bS > ${sample_name}.bam
    """
}

process quant_samples_oarfish {
  publishDir "sequin/quants/oarfish"

  input:
    tuple val(sample_name), path(ref_bam)
  output:
    path "${sample_name}_quant.tsv"

  """
  ${params.oarfish} -a ${ref_bam} --max-em-iter 100 -t 50 -o ${sample_name}_quant.tsv
  """
}

process quant_samples_nanocount {
  publishDir "sequin/quants/nanocount"

  input:
    tuple val(sample_name), path(ref_bam)
  output:
    path "${sample_name}_quant.tsv"

  """
  NanoCount -i ${ref_bam} -x -o ${sample_name}_quant.tsv
  """
}

workflow {
  def sequin_sample_paths = Channel.fromList(params.sequin_samples)
    .map{ it -> [it, "${params.sequin_input_path}/${it}_MinION_sequencing.fastq.gz"] }

  // map the samples with minimap2
  map_samples(sequin_sample_paths)

  // quantify with both oarfish and NanoCount
  quant_samples_nanocount(map_samples.out)
  quant_samples_oarfish(map_samples.out)
}
