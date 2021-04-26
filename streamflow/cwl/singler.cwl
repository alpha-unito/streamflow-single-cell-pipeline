#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
baseCommand: ["Rscript", "/r-scripts/doSingleR.R"]
arguments:
  - position: 4
    valueFrom: '$(runtime.outdir)/$(inputs.sample_id)'
inputs:
  rds:
    type: File
    inputBinding:
      position: 1
  sample_id:
    type: string
    inputBinding:
      position: 2
  species:
    type: string
    inputBinding:
      position: 3
  threads:
    type: int
    inputBinding:
      position: 5
outputs:
  analysis:
    type: File
    outputBinding:
      glob: '$(inputs.sample_id)/singleR.rds'
