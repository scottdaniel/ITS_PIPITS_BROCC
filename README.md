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
  - `ITS_primer_trim` (https://github.com/PennChopMicrobiomeProgram/Primer_trim)
  - To install (dnabc as example):
  ```bash
  git clone https://github.com/PennChopMicrobiomeProgram/dnabc
  cd dnabc
  conda activate ITS_PIPITS_BROCC
  pip install -e ./
  ```

## Required input files for the pipeline
To run the pipeline, we need
- Multiplexed R1/R2 read pairs

## How to run
- Create a project directory, e.g. `/scr1/users/tuv/ITS_Run1`
- Copy the files from this repository into that directory
- Edit `config.yml` so that it suits your project. In particular,
  - **all: project_dir**: path to the project directory, e.g. `"/scr1/users/tuv/ITS_Run1"`
  - **all: ncbi_db**: path to a local ncbi nt database, e.g. `"/path/to/nt"`
  - **all: mux_dir**: the directory containing multiplexed Illumina sequencing reads, which does not have to be in the project directory, e.g. `"/path/to/mux_files"` 
  - **all: demux_dir**: the directory to contain the demultiplexed R1/R2 read pairs
  - **all: threads**: number of threads to use
  - **all: ITS_subregion**: can leave blank or one of `ITS1` or `ITS2` for ITS subregion extraction
  - **all: mapping_file**: Mapping file of samples with barcode information for demultiplexing
  - **demux: mismatch**: Number of allowable basepair mismatches on barcode sequence for demultiplexing
  - **demux: revcomp**: If `TRUE`, reverse complement barcode sequence before demultiplexing
  - **trim: f_primer**: Sequence of forward primer used for ITS PCR
  - **trim: r_primer**: Sequence of reverse primer used for ITS PCR
  - **trim: mismatch**: Number of allowable basepair mismatches on ITS PCR primers for trimming
  - **trim: min_length**: Minimum length of primer to trim from reads
- To run the pipeline, activate the environment by entering `conda activate ITS_PIPITS_BROCC`, `cd` into the project directory and execute:
```bash
snakemake \
    --configfile path/to/config.yml \
    --keep-going \
    --latency-wait 90 \
    --notemp
```
- When submitting jobs using `qsub`, you may run `qsub run_snakemake.bash config.yml`
- You can use the [skeleton.Rmd](Rmd/skeleton.Rmd) to create a basic bioinformatic report from the results
  
## Notes on BROCC
`create_local_taxonomy_db.py` may be used to install a local taxonomy db for faster processing

## Rules
### Demultiplexing
Input: Multiplexed Illumina sequencing files  
Output: manifest.csv, total_read_counts.tsv, demultiplexed fastq files
### Trimming
#### Primer trim
Removes ITS forward and reverse primer sequences from reads
#### Read cleaning
Filter out reads with less than 50 bps after trimming
#### Vsearch filter
Remove reads that still retain ITS forward and reverse primer sequences in the reads. These reads are low quality due to the formation of primer dimers, and usually have the reverse complement of the primer sequence at the beginning of the read.
#### Pair reads
Only use reads with a forward and reverse sequence
### PIPITS
Run PIPITS pipeline to process reads and form OTUs clusters with or witout extracting the ITS subregion
### BROCC
Determine the taxonomic assignments of the OTUs by through a consensus based BLAST result
