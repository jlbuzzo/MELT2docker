############################## HEADER ##########################################
#
# This is a template configuration file.
#
################################################################################





############################## APP GENERAL INFO ################################

#ROOT_DIR				:= /home/users/jlbuzzo/projects_old/MELT2docker
#RUNNING_DIR			:= $(prefix)

APP						:= MELT2docker
VERSION					:= v0.1
ID						:= npp_123456
KEY						:= $(ROOT_DIR)/keys/key.pub
LICENSE					:= $(ROOT_DIR)/LICENSE
MANUAL					:= $(ROOT_DIR)/doc/neopipe.pod
SITE					:= https://www.docker.com/
REPO					:= https://download.docker.com/linux/static/stable/x86_64/docker-18.06.1-ce.tgz
CLOUD					:= https://www.aws.com/galantelab/neopipe

CONFIGFILE				:= $(RUNNING_DIR)/config.mk
ENCAPSULATE				:=
CLOUDERIZE				:=
CLOUDFILE				:= $(RUNNING_DIR)/aws.yml
DOCKERIZE				:=
DOCKERFILE				:= $(RUNNING_DIR)/docker/Dockerfile
IMAGE					:= galantelab/MELT2docker
MAKERIZE				:=
MAKEFILE				:= $(RUNNING_DIR)/Meltfile
TARBALL					:=

CMD						:=



############################## APP SPECIFIC VARIABLES ##########################

# General pipeline's data files: suffixes.
SFX						:= .bam
INPUT					?= list.txt

# Mount some important files in their respective folders.
# Must use only absolute paths!
MP_CONFIG_l				:= $(MAKEFILE_LIST)
MP_INPUTS_l				:= $(INPUT)
MP_OUTPUTS_l			:=
MP_ASSETS_l				:= /home/scratch60/lbuzzo_13_may/data/hg38_hash.perldb
MP_REFERENCE_l			:= /home/scratch60/lbuzzo_13_may/data/hg38.fa
MP_ANNOTATION_l			:= /home/scratch60/lbuzzo_13_may/data/gencode.v29.annotation.gtf
MP_EXTRA_l				:=
MP_TMP_l				:= /home/scratch60/lbuzzo_13_may/tmp


# App runtime Arguments.
ARGS					:=


# Dicovery directories for each result.
BASE_DISCOVERY_d		:= $(OUTPUTS_d)/results
HERVK_DISCOVERY_d		:= $(BASE_DISCOVERY_d)/HERVK
LINE1_DISCOVERY_d		:= $(BASE_DISCOVERY_d)/LINE1
ALU_DISCOVERY_d			:= $(BASE_DISCOVERY_d)/ALU
SVA_DISCOVERY_d			:= $(BASE_DISCOVERY_d)/SVA



############################## MISCELLANEOUS ###################################
# This section must contain some common variables, like reference genomes and
# annotations files.


# Reference files and annotations.
REFERENCE_GENOME		:= $(MP_REFERENCE_l)
REFERENCE_GTF			:= $(MP_ANNOTATION_l)
REFERENCE_BED			:= $(MELT_PLACE)/add_bed_files/Hg38/Hg38.genes.bed
REP_ANNOTATION			:= $(ANNOTATION_d)/rep.hg38.converted.bed
REP_ANNOTATION_manual	:= $(ANNOTATION_d)/rep.hg38.converted.manual.bed

BAM_COVERAGE			:= 40
