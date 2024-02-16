process map_samples {
  debug true
  cpus 9

  input:
    tuple val(seq_type), val(sample_name), path(sample)
  output:
    tuple val(seq_type), val(sample_name), path("${sample_name}.bam")

    """
    echo ${seq_type} ${params.bambu.ref} ${sample_name} ${sample}
    minimap2 -t 8 -ax map-ont -N 100 ${params.bambu.ref} ${sample} | samtools view -@2 -bS > ${sample_name}.bam
    """
}

process quant_samples_oarfish {
  publishDir "bambu_data/quants/oarfish"

  input:
    tuple val(seq_type), val(sample_name), path(ref_bam)
  output:
    path "${sample_name}_quant.tsv"

  script:
    if (seq_type == "drna")
      """
      ${params.oarfish} -a ${ref_bam} -t 50 -o ${sample_name}_quant.tsv
      """
    else if (seq_type == "cdna")
       """
      ${params.oarfish} -a ${ref_bam} -n -o ${sample_name}_quant.tsv
      """     
}

process quant_samples_nanocount {
  publishDir "bambu_data/quants/nanocount"

  input:
    tuple val(seq_type), val(sample_name), path(ref_bam)
  output:
    path "${sample_name}_quant.tsv"

  script:
    if (seq_type == "drna")
      """
      NanoCount -x -i ${ref_bam} -o ${sample_name}_quant.tsv
      """
    else if (seq_type == "cdna")
      """
      NanoCount -t -1 -n -x -i ${ref_bam} -o ${sample_name}_quant.tsv
      """     
}

workflow {
  def cdna_sample_paths = Channel.fromList(params.bambu.samples.cdna)
    .map{ it -> ["cdna", it, "${params.bambu.data_path}/${it}.fastq.gz"] }
  def drna_sample_paths = Channel.fromList(params.bambu.samples.drna)
    .map{ it -> ["drna", it, "${params.bambu.data_path}/${it}.fastq.gz"] }
  
  def sample_paths = cdna_sample_paths.concat(drna_sample_paths)

  // map the samples with minimap2
  map_samples(sample_paths)

  // quantify with both oarfish and NanoCount
  quant_samples_nanocount(map_samples.out)
  quant_samples_oarfish(map_samples.out)
}
