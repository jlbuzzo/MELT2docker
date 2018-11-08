############################## HEADER ##########################################
#
# This is a template configuration file.
#
################################################################################



############################## INFRASTRUCTURE ##################################
# This section must contain only general infrastructure variables.


# NOTE: Varable names can end with '_d', '_f', '_r', '_b', '_t' and 'l'.
# They represent, respectively, directories, files, references, booleans,
# temporaries and lists values inside them.
# Only names ending with '_l' can have non-single string values, but you can
# compose'em: OUTPUTS_dlt (for a list of temporary directories). 


# There could be an user for the container (ex. lion).
# So, C_BASE_d must be set accordingly:
#DEFAULT_USER			:= $(if $(SWITCH),lion,)


# Host's base directory, outside the container: $(PWD).
H_BASE_d				:= $(if $(prefix),$(prefix),$(PWD))
H_CONFIG_d				:= $(H_BASE_d)/config
H_INPUTS_d				:= $(H_BASE_d)/inputs
H_OUTPUTS_d				:= $(H_BASE_d)/outputs
H_ASSETS_d				:= $(H_BASE_d)/assets
H_REFERENCE_d			:= $(H_BASE_d)/reference
H_ANNOTATION_d			:= $(H_BASE_d)/annotation
H_EXTRA_d				:= $(H_BASE_d)/extra
H_TMP_d					:= $(H_BASE_d)/tmp


# Container's base directory: /home.
C_BASE_d				:= /home/$(DEFAULT_USER)
C_CONFIG_d				:= $(C_BASE_d)/config
C_INPUTS_d				:= $(C_BASE_d)/inputs
C_OUTPUTS_d				:= $(C_BASE_d)/outputs
C_ASSETS_d				:= $(C_BASE_d)/assets
C_REFERENCE_d			:= $(C_BASE_d)/reference
C_ANNOTATION_d			:= $(C_BASE_d)/annotation
C_EXTRA_d				:= $(C_BASE_d)/extra
C_TMP_d					:= $(C_BASE_d)/tmp


# Overwrite values with caution!
H_OUTPUTS_d				:= /home/scratch60/lbuzzo/MELT/rett



############################## APP SPECIFIC VARIABLES ##########################
# This section must contain variables specific to the application.


# General pipeline data: name and file suffixes.
PIPELINE				:= MELT
SUFFIXES				:= .bam


# A common internal representation. The best place to overwrite variables.
BASE_d					:= $(if $(SWITCH),$(C_BASE_d),$(H_BASE_d))
CONFIG_d				:= $(if $(SWITCH),$(C_CONFIG_d),$(H_CONFIG_d))
INPUTS_d				:= $(if $(SWITCH),$(C_INPUTS_d),$(H_INPUTS_d))
OUTPUTS_d				:= $(if $(SWITCH),$(C_OUTPUTS_d),$(H_OUTPUTS_d))
ASSETS_d				:= $(if $(SWITCH),$(C_ASSETS_d),$(H_ASSETS_d))
REFERENCE_d				:= $(if $(SWITCH),$(C_REFERENCE_d),$(H_REFERENCE_d))
ANNOTATION_d			:= $(if $(SWITCH),$(C_ANNOTATION_d),$(H_ANNOTATION_d))
EXTRA_d					:= $(if $(SWITCH),$(C_EXTRA_d),$(H_EXTRA_d))
TMP_d					:= $(if $(SWITCH),$(C_TMP_d),$(H_TMP_d))


# Mount some important files in their respective folders. Must use absolute paths!
MP_CONFIG_l				:= $(PWD)/Makefile $(PWD)/config.mk
MP_INPUTS_l				:= $(INPUTS_d)/ponga.txt
MP_OUTPUTS_l			:=
MP_ASSETS_l				:= /home/scratch60/lbuzzo/RTC/neopipe/assets/ref.perldb
MP_REFERENCE_l			:= /home/genomes/Homo_sapiens/hg38/hg38.fa
MP_ANNOTATION_l			:= /home/projects2/databases/gencode/release29/gencode.v29.annotation.gff3.gz
MP_EXTRA_l				:=
MP_TMP_l				:=


# Some important input files.
MELT_PLACE				:= $(BASE_d)/MELTv2.1.4
MEI_PLACE				:= $(MELT_PLACE)/me_refs/Hg38/
HERVK_ZIP				:= $(MEI_PLACE)/HERVK_MELT.zip
LINE1_ZIP				:= $(MEI_PLACE)/LINE1_MELT.zip
ALU_ZIP					:= $(MEI_PLACE)/ALU_MELT.zip
SVA_ZIP					:= $(MEI_PLACE)/SVA_MELT.zip

# Dicovery directories for each result.
BASE_DISCOVERY_d		:= $(OUTPUTS_d)/results
HERVK_DISCOVERY_d		:= $(BASE_DISCOVERY_d)/HERVK
LINE1_DISCOVERY_d		:= $(BASE_DISCOVERY_d)/LINE1
ALU_DISCOVERY_d			:= $(BASE_DISCOVERY_d)/ALU
SVA_DISCOVERY_d			:= $(BASE_DISCOVERY_d)/SVA


# Some commands.
JAR						:= $(MELT_PLACE)/MELT.jar
PREPROCESS				:= java -Xmx2G -jar $(JAR) Preprocess
INDIV_ANALYSIS			:= java -Xmx4G -jar $(JAR) IndivAnalysis
GROUP_ANALYSIS			:= java -Xmx6G -jar $(JAR) GroupAnalysis
GENOTYPE				:= java -Xmx2G -jar $(JAR) Genotype
MAKE_VCF				:= java -Xmx2G -jar $(JAR) MakeVCF




############################## MISCELLANEOUS ###################################
# This section must contain some common variables, like reference genomes and
# annotations files.


# Reference files and annotations.
GENOMES_d				:= $(REFERENCE_d)/genomes
REP_ANNOTATION			:= $(ANNOTATION_d)/rep.hg38.converted.bed
REP_ANNOTATION_manual	:= $(ANNOTATION_d)/rep.hg38.converted.manual.bed
REFERENCE_GTF			:= $(ANNOTATION_d)/gencode.v26.annotation.gtf
REFERENCE_GENOME		:= $(MP_REFERENCE_l)
REFERENCE_BED			:= $(MELT_PLACE)/add_bed_files/Hg38/Hg38.genes.bed
BAM_COVERAGE			:= 40



################################################################################
# This is the superfulous variables section.



# Some other variables.
TEMP_PROCESS_DIR		:= $(OUTPUTS_d)/temp/
DUMP_DIR				:= $(OUTPUTS_d)/result/dump/
SEARCH_CRIT				:= $(BASE_d)/scripts/search_crit.awk
