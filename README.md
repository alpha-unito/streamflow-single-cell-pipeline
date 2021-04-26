# Single-cell RNA sequencing pipeline

This example showcases how [StreamFlow](https://streamflow.di.unito.it) can be used to manage a single-cell RNA sequencing workflow in a hybrid HPC-Cloud environment. The pipeline was originally described in [this article](https://doi.org/10.1109/TETC.2020.3019202):

```text
I. Colonnelli, B. Cantalupo, I. Merelli and M. Aldinucci, “StreamFlow: cross-breeding cloud with HPC,” in IEEE Transactions on Emerging Topics in Computing, doi: 10.1109/TETC.2020.3019202.
```

## Preliminary steps

### Prepare data

All the pipeline inputs are saved in the `streamflow/cwl/data` folder, with the exception of the `genome` folder and the `input` folder. In order to download them, use the related script

```bash
streamflow/cwl/data/download_data.sh
```

### Docker containers

In the original setting, the pipeline has been executed on two Docker containers:

- The `cellranger` container, which has been executed on the Occam HPC facility at University of Torino;
- The `r-environment` container, that has been deployed on a Kubernetes cluster.

Both containers can be recreated using the corresponding `build` script in the `docker/cellranger` and `docker/r-environment` folders, respectively.

Please note that, in order to build the CellRanger container, you must download the correct version of the [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.1) software (distributed under [this license](https://support.10xgenomics.com/docs/license)), and the correct version of the [bcl2fastq2](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/downloads.html) software (distributed under [this license](https://jp.support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-eula.pdf)).

### Running the HPC part

Unfortunately, the Occam facility is not publicly available, and this prevents the reproducibility of the pipeline in the same exact configuration reported in the original article. Nevertheless, it is still possible to run the first part of the pipeline on a different HPC facility and then move the data to a Kubernetes cluster.

Since normally HPC facilities do not support Docker containers, we also provide a SLURM configuration under `streamflow/environment/slurm`, which relies on a Singularity image. The Singularity image can be built using the recipe in the `singularity/cellranger` folder, which only converts the Docker container to the Singularity format.

Then, the resulting `*.sif` image can be manually copied to the HPC facility and the `streamflow/environment/slurm/slurm_template.jinja2` file can be updated to point to the correct image path. Finally, it is necessary to modify the `slurm-single-cell` model in the `streamflow/streamflow.yml` file with the correct connection details to access the data center.

It is also worth noting that the StreamFlow `workdir` should point to a shared portion of the HPC file-system. Indeed, StreamFlow only transfers files from the local machine to the frontend node, but SLURM-managed worker nodes must be able to see them in order to work properly. By default, StreamFlow stages all inputs and outputs under the `/tmp/streamflow` path. If this path is mounted as a shared file system, please modify it accordingly in the StreamFlow file.

## Run the example

In order to run the example, you have to [install StreamFlow](https://github.com/alpha-unito/streamflow#use-streamflow) on your machine. Then, you can launch the execution with the following command

```bash
streamflow streamflow/streamflow.yml
```

In the original paper, we pre-loaded all the data to the HPC facility, in order to simulate a scenario in which a sequencing machine writes its output files directly on the data center file-system.

Nevertheless, it is also possible to let StreamFlow transfer the files from the local machine to the HPC login node, and the `streamflow/cwl/config.yml` file points to local files by default. If data have been manually pre-loaded, simply put remote paths in the `streamflow/cwl/config.yml` before starting the execution.

For debugging and experimentation purposes, it is also possible to run the entire pipeline on the local machine using the `docker-compose` environment. To do that, simply modify the `streamflow/streamflow.yml` file to bind each workflow step with the `docker-compose-single-cell` model. Please be aware that CellRanger requires a large amount of RAM (at least 64GB, 128GB recommended) to perform the `count` step, making it impossible to complete the pipeline on a typical desktop machine.
