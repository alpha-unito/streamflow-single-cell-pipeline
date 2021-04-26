$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: seurat
label: seurat
class: CommandLineTool
cwlVersion: v1.0
baseCommand: ["Rscript", "/r-scripts/Template_basic_example.R"]
arguments:
  - position: 2
    valueFrom: "$(runtime.outdir)"
  - position: 4
    valueFrom: "/r-scripts/seurattools"
  - position: 10
    valueFrom: "/usr/bin/python3"
inputs:
  matrix_dir:
    type: Directory
    inputBinding:
      position: 1
  sample_id:
    type: string
    inputBinding:
      position: 3
  max_pc_mito:
    type: int
    inputBinding:
      position: 5
  min_feature:
    type: int
    inputBinding:
      position: 6
  max_feature:
    type: int
    inputBinding:
      position: 7
  num_pc:
    type: int
    inputBinding:
      position: 8
  cc_rds:
    type: File
    inputBinding:
      position: 9
  mito_prefix:
    type: string
    inputBinding:
      position: 11
outputs:
  rds:
    type: File
    outputBinding:
      glob: "$(inputs.sample_id)/$(inputs.sample_id)_logNormalized.rds"
