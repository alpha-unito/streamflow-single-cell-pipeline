#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
baseCommand: [cellranger, mkfastq]
arguments:
  - position: 3
    prefix: --output-dir
    valueFrom: $(runtime.outdir)
  - position: 4
    prefix: --project
    valueFrom: "single-cell"
inputs:
  samplesheet:
    type: File
    inputBinding:
      position: 1
      prefix: --csv
  inp_dir:
    type: Directory
    inputBinding:
      position: 2
      prefix: --run
outputs:
  fastqs:
    type: Directory[]
    outputBinding:
      glob: 'single-cell/C[0-9]'
