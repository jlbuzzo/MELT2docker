############################## HEADER ##########################################
#
# This is an template configuration file.
#
################################################################################



############################## DATA ############################################

# NOTE: Varable names can end with '_d', '_f', '_r', '_b', '_t' and 'l', to
# represent, respectively, directories, files, references, booleans, temporaries
# and lists values in them.
# Only names ending with '_l' can have non-single string values, but you can
# compose'em: OUTPUTS_dlt (for a list of temporary directories). 


# General pipeline data: name and file suffixes.
PIPELINE 		:= MELT
SUFFIXES 		:= .bam


# Main command line arguments.
#ARGS 			:=


# Host's base directory, outside the container: $(PWD).
H_BASE_d		:= $(PWD)
H_CONFIG_d		:= $(H_BASE_d)/config
H_INPUTS_d		:= $(H_BASE_d)/inputs
H_OUTPUTS_d		:= $(H_BASE_d)/outputs
H_ASSETS_d		:= $(H_BASE_d)/assets
H_REFERENCE_d	:= $(H_BASE_d)/reference
H_ANNOTATION_d	:= $(H_BASE_d)/annotation
H_EXTRA_d		:= $(H_BASE_d)/extra
H_TMP_d			:= $(H_BASE_d)/tmp


# Container's base directory: /home.
C_BASE_d		:= /home
C_CONFIG_d		:= $(C_BASE_d)/config
C_INPUTS_d		:= $(C_BASE_d)/inputs
C_OUTPUTS_d		:= $(C_BASE_d)/outputs
C_ASSETS_d		:= $(C_BASE_d)/assets
C_REFERENCE_d	:= $(C_BASE_d)/reference
C_ANNOTATION_d	:= $(C_BASE_d)/annotation
C_EXTRA_d		:= $(C_BASE_d)/extra
C_TMP_d			:= $(C_BASE_d)/tmp


# A common internal representation.
BASE_d			:= $(if $(SWITCH),$(C_BASE_d),$(H_BASE_d))
CONFIG_d		:= $(if $(SWITCH),$(C_CONFIG_d),$(H_CONFIG_d))
INPUTS_d		:= $(if $(SWITCH),$(C_INPUTS_d),$(H_INPUTS_d))
OUTPUTS_d		:= $(if $(SWITCH),$(C_OUTPUTS_d),$(H_OUTPUTS_d))
ASSETS_d		:= $(if $(SWITCH),$(C_ASSETS_d),$(H_ASSETS_d))
REFERENCE_d		:= $(if $(SWITCH),$(C_REFERENCE_d),$(H_REFERENCE_d))
ANNOTATION_d	:= $(if $(SWITCH),$(C_ANNOTATION_d),$(H_ANNOTATION_d))
EXTRA_d			:= $(if $(SWITCH),$(C_EXTRA_d),$(H_EXTRA_d))
TMP_d			:= $(if $(SWITCH),$(C_TMP_d),$(H_TMP_d))




############################## ANCILLARY VARIABLES ############################

# Results.
MELT_PLACE		:= $(BASE_d)/MELTv2.1.4
MEI_PLACE		:= $(MELT_PLACE)/me_refs/Hg38/
HERVK_ZIP		:= $(MEI_PLACE)/HERVK_MELT.zip
LINE1_ZIP		:= $(MEI_PLACE)/LINE1_MELT.zip
ALU_ZIP			:= $(MEI_PLACE)/ALU_MELT.zip
SVA_ZIP			:= $(MEI_PLACE)/SVA_MELT.zip

BASE_DISCOVERY_d	:= $(BASE_d)/outputs/result
HERVK_DISCOVERY_d	:= $(BASE_DISCOVERY_d)/HERVK
LINE1_DISCOVERY_d	:= $(BASE_DISCOVERY_d)/LINE1
ALU_DISCOVERY_d		:= $(BASE_DISCOVERY_d)/ALU
SVA_DISCOVERY_d		:= $(BASE_DISCOVERY_d)/SVA


# Some commands.
JAR				:= $(MELT_PLACE)/MELT.jar
PREPROCESS		:= java -Xmx2G -jar $(JAR) Preprocess
INDIV_ANALYSIS	:= java -Xmx4G -jar $(JAR) IndivAnalysis
GROUP_ANALYSIS	:= java -Xmx6G -jar $(JAR) GroupAnalysis
GENOTYPE		:= java -Xmx2G -jar $(JAR) Genotype
MAKE_VCF		:= java -Xmx2G -jar $(JAR) MakeVCF


# Temporary directories for runtime processing.
TEMP_PROCESS_DIR	:= $(OUTPUTS_d)/temp/
DUMP_DIR			:= $(OUTPUTS_d)/result/dump/




############################## EDIT HERE! #####################################

# Reference files and annotations.
GENOMES_d				:= $(REFERENCE_d)/genomes
REP_ANNOTATION			:= $(ANNOTATION_d)/rep.hg38.converted.bed
REP_ANNOTATION_manual	:= $(ANNOTATION_d)/rep.hg38.converted.manual.bed
REFERENCE_GTF			:= $(ANNOTATION_d)/gencode.v26.annotation.gtf
REFERENCE_GENOME		:= $(GENOMES_d)/hg38/hg38.fa
REFERENCE_BED			:= $(MELT_PLACE)/add_bed_files/Hg38/Hg38.genes.bed

###############################################################################



# Some other variables.
BAM_COVERAGE	:= 40
SEARCH_CRIT		:= $(BASE_d)/scripts/search_crit.awk

# Pay attention here!
INPUTS			:= /home/scratch60/lbuzzo/MELT/inputs/bams.txt
OUTPUTS_d		:= /home/scratch60/lbuzzo/MELT/outputs
