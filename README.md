# DIOspyros scans for REpeated SELection (dioresel)

## Introduction

This repository is part of my Master's project on finding repeated adapted genes during an adaptive radiation. The involved workflow includes different tests for adaptation (pcadapt, Fst, McDonald-Kreitman test).

## Usage

1. Install [miniforge3 (version: 25.3.0)](https://github.com/conda-forge/miniforge/releases/tag/25.3.0-3).

2. Clone this repository:

```
git clone -b dev git@github.com:flotji-potji/dioresel.git
cd dioresel
```

3. Create the environment:

```
mamba env create -f workflow/envs/conda_env.yaml
```

4. Activate the environment:

```
conda activate dioresel
```

5. Run the analysis locally:

```
snakemake -c 1 --use-conda
```

6. Run the analysis on a HPC:

```
snakemake -c 1 --use-conda --profile config/slurm
```

It is advised to adjust resource usage of jobs to match their cluster performance. This can be done in the `config.yaml` file in `config/slurm`, just adjust the rules individual resources. This workflow was tested on RedHat Linux 11 using the Life Science Compute Cluster at the University of Vienna.
