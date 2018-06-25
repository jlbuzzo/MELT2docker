############################## HEADER #########################################
#
# This is an template configuration file.
#
###############################################################################



############################## DATA ###########################################

# General pipeline data: name and file suffixes.
PIPELINE := MELT
SUFFIXES := .bam
OUTER_BASE_DIR := /home/users/jlbuzzo/projects_old/melt


# Essential locales, inside the container.
#BASE_DIR := $(OUTER_BASE_DIR)
BASE_DIR := $(PWD)
CONFIG := $(BASE_DIR)/config
INPUTS := $(BASE_DIR)/inputs
OUTPUTS := $(BASE_DIR)/outputs
ASSETS := $(BASE_DIR)/assets
REFERENCE := $(BASE_DIR)/reference
ANNOTATION := $(BASE_DIR)/annotation

# Results.
MELT_PLACE := $(BASE_DIR)/MELTv2.1.4
MEI_PLACE := $(MELT_PLACE)/me_refs/Hg38/
HERVK_ZIP := $(MEI_PLACE)/HERVK_MELT.zip
LINE1_ZIP := $(MEI_PLACE)/LINE1_MELT.zip
ALU_ZIP := $(MEI_PLACE)/ALU_MELT.zip
SVA_ZIP := $(MEI_PLACE)/SVA_MELT.zip

BASE_DISCOVERY := $(BASE_DIR)/outputs/result
HERVK_DISCOVERY_DIR := $(BASE_DISCOVERY)/HERVK
LINE1_DISCOVERY_DIR := $(BASE_DISCOVERY)/LINE1
ALU_DISCOVERY_DIR := $(BASE_DISCOVERY)/ALU
SVA_DISCOVERY_DIR := $(BASE_DISCOVERY)/SVA


# Some commands.
JAR := $(MELT_PLACE)/MELT.jar
PREPROCESS := java -Xmx2G -jar $(JAR) Preprocess
INDIV_ANALYSIS := java -Xmx4G -jar $(JAR) IndivAnalysis
GROUP_ANALYSIS := java -Xmx6G -jar $(JAR) GroupAnalysis
GENOTYPE := java -Xmx2G -jar $(JAR) Genotype
MAKE_VCF := java -Xmx2G -jar $(JAR) MakeVCF


# Pay attention here!
INPUT := $(INPUTS)
OUTPUT_DIR := $(OUTPUTS)


# Reference files and annotations.
GENOMES_DIR := $(REFERENCE)/genomes
ANNOTATION_DIR := $(REFERENCE)/annotation
REF_ANNOTATION := $(ANNOTATION_DIR)/gencode.v26.annotation.gtf
REP_ANNOTATION := $(ANNOTATION_DIR)/rep.hg38.converted.bed
REP_ANNOTATION_manual := $(ANNOTATION_DIR)/rep.hg38.converted.manual.bed
REFERENCE_GENOME_FASTA := $(GENOMES_DIR)/hg38/hg38.fa
REFERENCE_BED := $(MELT_PLACE)/add_bed_files/Hg38/Hg38.genes.bed
BAM_COVERAGE := 40


# Temporary directories for runtime processing.
TEMP_PROCESS_DIR := $(OUTPUT_DIR)/temp/
DUMP_DIR := $(OUTPUT_DIR)/result/dump/


# Some other variables.
DISTANCE := 750000
MIN_DIST := 1000000
SEARICH_CRIT := $(BASE_DIR)/scripts/search_crit.awk
