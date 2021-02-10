import configparser
import yaml

PROJECT_DIR = config["all"]["project_dir"]
PIPITS_DIR = PROJECT_DIR + "/PIPITS_output"

include: "rules/targets.rules"
include: "rules/PIPITS.rules"

workdir: PROJECT_DIR

rule all:
  input: TARGET_ALL
