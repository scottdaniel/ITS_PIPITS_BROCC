import configparser
import yaml

from scripts import util_functions

PROJECT_DIR = config["all"]["project_dir"]
DEMUX_DIR = PROJECT_DIR + "/" + config["all"]["demux_dir"]
PIPITS_DIR = PROJECT_DIR + "/PIPITS_output"
BL_BR_DIR = PROJECT_DIR + "/BLAST_BROCC_output"
MAPPING_FP = PROJECT_DIR + "/" + config["all"]["mapping_file"]
SAMPLE_IDS = util_functions.get_sample(MAPPING_FP)

include: "rules/targets.rules"
include: "rules/PIPITS.rules"
include: "rules/blast_n_brocc.rules"
include: "rules/demux.rules"

workdir: PROJECT_DIR

rule all:
  input: TARGET_ALL
