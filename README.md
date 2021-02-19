# ITS_PIPITS_BROCC
This is a Snakemake pipeline for analyzing fungal internal transcribed spacer (ITS) sequences using PIPITS (https://github.com/hsgweon/pipits) and BROCC (https://github.com/kylebittinger/brocc)

## Installation
To install, we assume you already have installed `Miniconda3` (https://docs.conda.io/en/latest/miniconda.html)
- Clone the repository:
```bash
git clone https://github.com/PennChopMicrobiomeProgram/ITS_PIPITS_BROCC.git
```
- Create a conda environment and install the required packages:
```bash
cd ITS_PIPITS_BROCC
conda create -n ITS_PIPITS_BROCC --channel bioconda --channel conda-forge --channel defaults python=3.6
conda install --name ITS_PIPITS_BROCC --file requirements.txt
/anaconda/envs/venv_name/bin/pip install brocc #brocc needs to be installed through your environment's pip
```

- The following software also need to be installed:
  - `dnabc` (https://github.com/PennChopMicrobiomeProgram/dnabc)
  - To install:
  ```bash
  git clone https://github.com/PennChopMicrobiomeProgram/dnabc
  cd dnabc
  conda activate ITS_PIPITS_BROCC
  pip install -e ./
  ```

## Required input files for the pipeline
To run the pipeline, we need
- Demultiplexed R1/R2 read pairs (either \*.fasta or \*.fasta.gz)

## How to run
- Create a project directory, e.g. `/scr1/users/tuv/ITS_Run1`
- Copy the files from this repository into that directory
- Edit `config.yml` so that it suits your project. In particular,
  - **all: project_dir**: path to the project directory, e.g. `"/scr1/users/tuv/ITS_Run1"`
  - **all: ncbi_db**: path to a local ncbi nt database, e.g. `"/path/to/nt"`
  - **all: demux_dir**: the directory containing demultiplexed R1/R2 read pairs, which does not have to be in the project directory, e.g. `"/path/to/demux_files"` 
  - **all: threads**: number of threads to use
  - **all: ITS_subregion**: can leave blank or one of `ITS1` or `ITS2` for ITS subregion extraction
- To run the pipeline, activate the environment by entering `conda activate ITS_PIPITS_BROCC`, `cd` into the project directory and execute:
```bash
snakemake \
    --configfile path/to/config.yml \
    --keep-going \
    --latency-wait 90 \
    --notemp
```
- When submitting jobs using `qsub`, you may run `qsub run_snakemake.bash config.yml`
  
## Notes on BROCC
`create_local_taxonomy_db.py` may be used to install a local taxonomy db for faster processing
