############################## HEADER #########################################
#
# This is an template configuration file.
#
###############################################################################



############################## DATA ###########################################

# General pipeline data: name and file suffixes.
PIPELINE := MELT
SUFFIXES := .bam


# Essential locales.
BASE_DIR := $(PWD)
CONFIG := $(BASE_DIR)/config
INPUTS := $(BASE_DIR)/inputs
OUTPUTS := $(BASE_DIR)/outputs
ASSETS := $(BASE_DIR)/assets
REFERENCE := $(BASE_DIR)/reference

# Pay attention here!
INPUT := $(INPUTS)
OUTPUT_DIR := $(OUTPUTS)

# Some locales and commands.
MELT_PLACE=$(BASE_DIR)/MELTv2.1.4
JAR := $(MELT_PLACE)/MELT.jar
PREPROCESS := java -Xmx2G -jar $(JAR) Preprocess
ANALYSE := java -Xmx2G -jar $(JAR) IndivAnalysis


# Some additional informations.
AUX_DIR := /home/users/jlbuzzo/projects_old/RTC/ancillaries
GENOMES_DIR := $(AUX_DIR)/genomes
ANNOTATION_DIR := $(AUX_DIR)/annotation
REF_ANNOTATION := $(ANNOTATION_DIR)/gencode.v26.annotation.gtf
REP_ANNOTATION := $(ANNOTATION_DIR)/rep.hg38.converted.bed
REP_ANNOTATION_manual := $(ANNOTATION_DIR)/rep.hg38.converted.manual.bed
REFERENCE_GENOME_FASTA := $(GENOMES_DIR)/hg38/hg38.fa
REFERENCE_BED := $(MELT_PLACE)/add_bed_files/Hg38/Hg38.genes.bed


SEARCH_CRIT := $(BASE_DIR)/scripts/search_crit.awk

# Temporary directories for runtime processing.
TEMP_PROCESS_DIR := $(OUTPUT_DIR)/temp/
DUMP_DIR := $(OUTPUT_DIR)/result/dump/


# Some other variables.
DISTANCE := 750000
MIN_DIST := 1000000
