import configparser
import yaml

from scripts import util_functions

PROJECT_DIR = config["all"]["project_dir"]
DEMUX_DIR = PROJECT_DIR + "/" + config["all"]["demux_dir"]
PRIMER_TRIM_FP = PROJECT_DIR + "/primer_trim"
V_TRIM_FP = PROJECT_DIR + "/vsearch_trim"
TRIM_DIR = PROJECT_DIR + "/final_trim"
READ_DIR = PROJECT_DIR + "/reads"
PIPITS_DIR = PROJECT_DIR + "/PIPITS_output"
BL_BR_DIR = PROJECT_DIR + "/BLAST_BROCC_output"
MAPPING_FP = PROJECT_DIR + "/" + config["all"]["mapping_file"]
SAMPLE_IDS = util_functions.get_sample(MAPPING_FP)
VSEARCH_DIR = PROJECT_DIR + "/VSEARCH_output"

include: "rules/targets.rules"
include: "rules/PIPITS.rules"
include: "rules/blast_n_brocc.rules"
include: "rules/demux.rules"
include: "rules/trim.rules"
include: "rules/vsearch.rules"

workdir: PROJECT_DIR

rule all:
  input: TARGET_ALL
