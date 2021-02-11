import configparser
import yaml

PROJECT_DIR = config["all"]["project_dir"]
PIPITS_DIR = PROJECT_DIR + "/PIPITS_output"
BL_BR_DIR = PROJECT_DIR + "/BLAST_BROCC_output"

include: "rules/targets.rules"
include: "rules/PIPITS.rules"
include: "rules/blast_n_brocc.rules"

workdir: PROJECT_DIR

rule all:
  input: TARGET_ALL
