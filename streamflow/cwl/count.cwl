#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
baseCommand: [cellranger, count]
arguments:
  - position: 1
    prefix: --id
    valueFrom: $(inputs.fastq.basename)
  - position: 3
    prefix: --fastqs
    valueFrom: $(inputs.fastq.path)
  - position: 4
    prefix: --sample
    valueFrom: $(inputs.fastq.basename)
inputs:
  genome:
    type: Directory
    inputBinding:
      position: 2
      prefix: --transcriptome
  fastq:
    type: Directory
outputs:
  matrix:
    type: Directory
    outputBinding:
      glob: '$(inputs.fastq.basename)/outs/filtered_feature_bc_matrix'
  sample_id:
    type: string
    outputBinding:
      outputEval: '$(inputs.fastq.basename)'