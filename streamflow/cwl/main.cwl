#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: Workflow
$namespaces:
  sf: "https://streamflow.org/cwl#"
inputs:
  cc_rds: File
  genome: Directory
  inp_dir: Directory
  max_pc_mito: int
  max_feature: int
  min_feature: int
  mito_prefix: string
  num_pc: int
  organism: string
  samplesheet: File
  singler_threads: int
requirements:
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
outputs:
  analyses:
    type:
      type: array
      items: File
    outputSource: subworkflow/analysis
steps:
  mkfastq:
    run: mkfastq.cwl
    in:
      inp_dir: inp_dir
      samplesheet: samplesheet
    out: [fastqs]
  subworkflow:
    run:
      class: Workflow
      inputs:
        cc_rds: File
        genome: Directory
        max_pc_mito: int
        max_feature: int
        min_feature: int
        mito_prefix: string
        num_pc: int
        organism: string
        singler_threads: int
        fastq: Directory
      outputs:
        analysis:
          type: File
          outputSource: singler/analysis
      steps:
        count:
          run: count.cwl
          in:
            genome: genome
            fastq: fastq
          out: [matrix, sample_id]
        seurat:
          run: seurat.cwl
          in:
            cc_rds: cc_rds
            matrix_dir: count/matrix
            sample_id: count/sample_id
            max_pc_mito: max_pc_mito
            max_feature: max_feature
            min_feature: min_feature
            mito_prefix: mito_prefix
            num_pc: num_pc
          out: [rds]
        singler:
          run: singler.cwl
          in:
            rds: seurat/rds
            sample_id: count/sample_id
            species: organism
            threads: singler_threads
          out: [analysis]
    scatter: fastq
    in:
      cc_rds: cc_rds
      fastq: mkfastq/fastqs
      genome: genome
      max_pc_mito: max_pc_mito
      max_feature: max_feature
      min_feature: min_feature
      mito_prefix: mito_prefix
      num_pc: num_pc
      organism: organism
      singler_threads: singler_threads
    out: [analysis]