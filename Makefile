############################## HEADER #########################################

# Makefile data:
#
# 		Default goal:	analyse_all
# 		File:			Makefile
# 		Module:			none
#		Project:		melt
# 		Author:			Jose L. L. Buzzo
# 		Organization:	RetroTeam
# 		Year:			2018


# Usage:
#
# 		Inputs: INPUT variable (string list).
# 					The string list must contain, at least, a path to a folder
# 					or a file in the system.
#
# 		Outputs: PROCESSED_INPUT variable (string list).
# 					The returned string list caontains all the existent files
# 					specified in the input args with canonical paths and
# 					filtered by the SUFFIXES auxilliary variable.
#
#		Options:





############################## PREAMBLE #######################################

# Include user's configuration file (only the first one).
CONFIG_FILE ?= ./config.mk
aux := $(word 1, $(abspath $(strip $(wildcard $(CONFIG_FILE)))))
include $(if $(aux), $(aux), $(error Configuration file not found))


# Presentation.
$(info ###############################################################################)
$(info .                             $(PIPELINE))
$(info ###############################################################################)


# Define a CALL message for the log.
CALL = [$(PIPELINE): $(shell date "+%Y-%m-%d(%H:%M:%S)")]


# ===> VALIDATIONS <===

# Validate INPUT against emptyness.
INPUT_PROCESSED := $(strip $(INPUT))
$(if $(INPUT_PROCESSED),, $(error Variabe INPUT is empty))

# Search and validate INPUT from a given string list, directory, or file, by
# the SUFFIXES list criteria.
#INPUT_PROCESSED := $(shell $(SCRIPTS)/ufinder.sh $(INPUT_PROCESSED) $(SUFFIXES))

# A limited alternative validation process (without 'ufinder' script).
INPUT_PROCESSED := $(if $(shell [ -f $(INPUT_PROCESSED) ] && echo 1), $(wildcard $(filter %$(SUFFIXES), $(shell cat $(INPUT_PROCESSED)))), $(wildcard $(INPUT_PROCESSED)/*$(SUFFIXES)))
INPUT_PROCESSED := $(abspath $(strip $(INPUT_PROCESSED)))

# Verifying 'SUFFIXES' terminated files' remaining after the filtering of the
# INPUT_PROCESSED list.
$(if $(INPUT_PROCESSED),, $(error No valid $(SUFFIXES) file specified))

# Verifying OUTPUT_DIR previous existence.
ifneq ($(wildcard $(OUTPUT_DIR)),)
$(info )
$(info $(CALL) Directory '$(OUTPUT_DIR)' will be overwritten!)
endif



# ===> PREPROCESSING <===

# ATTENTION: Variables that are string lists has an '_l' suffix appendedd to
# the end of their names. Their behavior differently in pattern rules.

# Carrefully set the INPUT_l variable after success validations.
INPUT_l := $(INPUT_PROCESSED)
INPUT_FILENAME_l := $(strip $(notdir $(INPUT_l)))
INPUT_BASENAME_l := $(strip $(basename $(notdir $(INPUT_l))))
INPUT_DIR_l := $(strip $(dir $(INPUT_l)))
#INPUT_DIR_BASENAME_l :=


#OUTPUT_l :=
#OUTPUT_FILENAME_l := $(INPUT_FILENAME_l)
#OUTPUT_BASENAME_l := $(INPUT_BASENAME_l)
OUTPUT_DIR := $(abspath $(strip $(OUTPUT_DIR)))# This is not a string list!

OUTPUT_DIR_BASENAME_l := $(addprefix $(OUTPUT_DIR)/, $(INPUT_BASENAME_l))
OUTPUT_l := $(strip $(join $(OUTPUT_DIR_BASENAME_l), $(addprefix /, $(INPUT_FILENAME_l))))

# List of '.abnormal', '.disc' '.disc.bai' and '.disc.fq' files, without their paths.
OUTPUT_ABNORMAL_l := $(addsuffix .abnormal, $(OUTPUT_l))
OUTPUT_DISC_l := $(addsuffix .disc, $(OUTPUT_l))
OUTPUT_DISC_BAI_l := $(addsuffix .disc.bai, $(OUTPUT_l))
OUTPUT_DISC_FQ_l:= $(addsuffix .disc.fq, $(OUTPUT_l))
OUTPUT_VCF_l := $(addsuffix .vcf, $(OUTPUT_l))


# Export all variables.
export



# ===> DEBUG CODE <===

# Enable some prints, if variable DBG="yes".
ifeq ($(DBG),yes)
$(info )
$(info INPUT_l: >$(INPUT_l).)
$(info INPUT_FILENAME_l: >$(INPUT_FILENAME_l).)
$(info INPUT_BASENAME_l: >$(INPUT_BASENAME_l).)
$(info INPUT_DIR: >$(INPUT_DIR_l).)
$(info )
$(info OUTPUT_DIR: >$(OUTPUT_DIR).)
$(info OUTPUT_DIR_BASENAME_l: >$(OUTPUT_DIR_BASENAME_l).)
$(info OUTPUT_l: >$(OUTPUT_l).)
$(info OUTPUT_ABNORMAL_l: >$(OUTPUT_ABNORMAL_l).)
$(info OUTPUT_DISC_BAI_l: >$(OUTPUT_DISC_BAI_l).)
## The pattern must be extended to all, must be constant, not an array of
## different values. See:
$(info )
$(info 1>$(patsubst /home/leonel/%$(SUFFIX), %.bam.o, $(INPUT)).)
$(info 2>$(patsubst $(OUTPUT_DIR)%.abnormal,%.abnormal.o, $(OUTPUT_ABNORMAL_l)).)
$(info )
endif

# Enable debug stop, if variable STP is non-empty.
ifdef STP
$(error Emergency stop)
endif



############################## TARGETS ########################################

all: process
	$(info )
	$(info $(CALL) All done!)

process: preprocess
	$(info )
	$(info $(CALL) Running MELT.)
	$(info $(ANALYSE))


preprocess: $(OUTPUT_DISC_l)
	$(info )
	$(info $(CALL) Target 'preprocess' complete!)


$(OUTPUT_DISC_l): %.bam.disc: %.bam %.bam.bai
	$(info )
	$(info $(CALL) Preprocessing file: $(word 1, $^).)
	$(info $(PREPROCESS) -bamfile $< -h $(REFERENCE_GENOME_FASTA).)


%.bam.bai: %.bam
	$(info )
	$(info $(CALL) Indexing file $<.)
	$(info >$(REQ_BAI)<)
	if [ -s "$(BAI)" ]; then \
		ln -sf "$(BAI)" $@; \
	else \
		samtools index -b -@ 8 $< $@; \
	fi


%.bam: $(INPUT_l)
	$(info )
	$(info $(CALL) Creating link for input file: $(REQ).)
	$(info >$(REQ)<)
	mkdir -p $(*D)
	if [ "$(CMD)" = "SO:coordinate" ]; then \
		ln -sf "$(BAM)" $@; \
	elif [ -s "$(STBAM)" ]; then \
		ln -sf "$(STBAM)" $@; \
	else \
		samtools sort -O BAM -m 8G -@ 8 $(REQ) -o $@; \
	fi



############################## LASY EVALUATION VARIABLES ######################

# Look for a *.sorted.bam.bai file, according to the link to 'OUTPUT_l'.
REQ_BAI = $(shell readlink -f $(filter %$(*F)$(SUFFIXES), $(OUTPUT_l)))
aux2 = $(strip $(filter %.sorted.bam, $(REQ_BAI)))
STBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*sorted*.bai')
NSTBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*.bai' -not -regex '.*sorted.*')
BAI = $(if $(aux2),$(STBAI),$(NSTBAI))

# Ancilary variables to this target.
REQ = $(filter %$(*F)$(SUFFIXES), $^)
CMD = $(shell samtools view -H $(REQ) 2> /dev/null | head -n1 | cut -f3)
BAM = $(shell readlink -f $(REQ))
STBAM = $(shell find $(dir $(BAM)) -type f -name '*$(*F)*.sorted.bam')



# Other independent targets.
dockerize: Dockerfile
	$(info )
	$(info Make a docker image.)
	docker build -t melt:latest .

dockerRun:
	$(info )
	$(info Run docker image with MELT command.)
	docker run \
		--rm \
		-u $$(id -u):$$(id -g) \
		-v \
		melt make


# .PHONY targets.
#.PHONY:
