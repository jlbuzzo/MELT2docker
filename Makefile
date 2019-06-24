#!/usr/bin/env bash
############################## HEADER ##########################################
#
# Makefile data:
#
# 		Default goal:	analyse_all
# 		File:			Makefile
# 		Module:			none
#		Project:		melt
#		Pipeline: 
# 		Author:			Jose L. L. Buzzo
# 		Organization:	RetroTeam
# 		Year:			2018
#
#
# Usage:
#
# 		Inputs: INPUTS variable (string list).
# 					The string list must contain, at least, a path to a folder
# 					or a file in the system.
#
# 		Outputs: PROCESSED_INPUTS variable (string list).
# 					The returned string list caontains all the existent files
# 					specified in the input args with canonical paths and
# 					filtered by the SUFFIXES auxilliary variable.
#
#		Options:
#
################################################################################





############################## PREAMBLE ########################################

# Command macros.
SHELL				:= bash
MODULE_NAME			:= main
#CALL				:= [$(MODULE_NAME): $(shell date --utc)]
ECHO				:= echo -e
MKDIR				:= mkdir -p
PING				:= ping -c1
DOCKER_RUN			:= docker run
DOCKER_FLAGS		:= --rm -ti



############################## INFRASTRUCTURE ##################################
# This section must contain only general infrastructure variables.


# NOTE: Varable names can end with '_d', '_f', '_r', '_t' and '_l'.
# They represent, respectively, directories, files, references, temporaries and
# lists values inside them.
# Only names ending with '_l' can have non-single string values, but you can
# compose'em: OUTPUTS_dlt (for a list of temporary directories). 

# Running mode.
DOCKERIZE			:= no
CLOUDRIZE			:= no

# Environment determination (container or host).
SWITCH_HAS			:= [ -f ".ponga_switch" ] && echo "1"
SWITCH				:= $(if $(shell $(SWITCH_HAS)),1,)


# Include user's configuration file (only the first one).
CONFIG_HAS			:= [ -f "config.mk" ] && echo "1"
CONFIG				:= $(if $(shell $(CONFIG_HAS)),config.mk,)

CONFIG_f 			?= $(if $(SWITCH),/home/config/config.mk,./config.mk)
aux 				:= $(word 1, $(abspath $(strip $(wildcard $(CONFIG_f)))))
include $(if $(aux), $(aux), $(error Configuration file not found))



# Infrastructure directories.
INFRASTRUCTURE_dl	:= $(CONFIG_d) $(INPUTS_d) $(OUTPUTS_d) $(ASSETS_d) $(REFERENCE_d) $(ANNOTATION_d) $(EXTRA_d) $(TMP_d)



############################## VALIDATIONS #####################################

# Presentation.
$(info ********************************************************************************)
$(info .                               $(PIPELINE))
$(info ********************************************************************************)


# Validate INPUTS against emptyness.
#$(info >>$(MP_INPUTS_l)<<2)
INPUTS_PROCESSED	:= $(word 1, $(abspath $(strip $(wildcard $(MP_INPUTS_l)))))
$(if $(INPUTS_PROCESSED),, $(error Variabe INPUTS is empty))

# Search and validate INPUTS from a given string list, directory, or file, by
# the SUFFIXES list criteria.
#INPUTS_PROCESSED := $(shell $(SCRIPTS)/ufinder.sh $(INPUTS_PROCESSED) $(SUFFIXES))

# A limited alternative validation process (without 'ufinder' script).
INPUTS_PROCESSED := $(if $(shell [ -f "$(INPUTS_PROCESSED)" ] && echo "1"),$(filter %$(SUFFIXES), $(shell cat $(INPUTS_PROCESSED))),$(wildcard $(INPUTS_PROCESSED)/*$(SUFFIXES)))
#INPUTS_PROCESSED := $(if $(shell [ -f "$(INPUTS_PROCESSED)" ] && echo "1"),$(info foi),$(info nao foi))

#$(info >>$(INPUTS_PROCESSED)<<2)
INPUTS_PROCESSED := $(abspath $(strip $(wildcard $(INPUTS_PROCESSED))))

# Verifying 'SUFFIXES' terminated files' remaining after the filtering of the
# INPUTS_PROCESSED list.
$(if $(INPUTS_PROCESSED),, $(error No valid $(SUFFIXES) file specified))


############################## PREPROCESSING ###################################

# ATTENTION: Variables that are string lists has an '_l' suffix appendedd to
# the end of their names. Their behavior differently in pattern rules.

OUTPUTS_d				:= $(abspath $(strip $(OUTPUTS_d)))# This is not a string list!

# Folder for timestamps.
TMSTP_d					:= $(OUTPUTS_d)/tmstp


# Carrefully set the INPUTS_l variable after success validations.
INPUTS_l				:= $(INPUTS_PROCESSED)
INPUTS_FILENAME_l		:= $(strip $(notdir $(INPUTS_l)))
INPUTS_BASENAME_l		:= $(strip $(basename $(notdir $(INPUTS_l))))
INPUTS_dl				:= $(strip $(dir $(INPUTS_l)))
#INPUTS_BASENAME_dl		:=


#OUTPUTS_l				:=
#OUTPUTS_FILENAME_l		:=
#OUTPUTS_BASENAME_l		:=
#OUTPUTS_dl				:=
OUTPUTS_BASENAME_dl		:= $(addprefix $(OUTPUTS_d)/, $(INPUTS_BASENAME_l))
OUTPUTS_l				:= $(strip $(join $(OUTPUTS_BASENAME_dl), $(addprefix /, $(INPUTS_FILENAME_l))))
#OUTPUTS_l				:= $(OUTPUTS_l:.bam=.sorted.bam)


# Derivated lists: '.bai', '.abnormal', '.disc', '.disc.bai' and others.
OUTPUTS_BAI_l			:= $(addsuffix .bai, $(OUTPUTS_l))
OUTPUTS_ABNORMAL_l		:= $(addsuffix .abnormal, $(OUTPUTS_l))
OUTPUTS_DISC_l			:= $(addsuffix .disc, $(OUTPUTS_l))
OUTPUTS_DISC_BAI_l		:= $(addsuffix .disc.bai, $(OUTPUTS_l))
OUTPUTS_DISC_FQ_l		:= $(addsuffix .disc.fq, $(OUTPUTS_l))

# Timestamp lists, terminated in '_lt'.
# Lists for individual analysis.
OUTPUTS_HERVK_INDIV_lt	:= $(addsuffix .hervk.indiv.tmstp, $(OUTPUTS_l))
OUTPUTS_LINE1_INDIV_lt	:= $(addsuffix .line1.indiv.tmstp, $(OUTPUTS_l))
OUTPUTS_ALU_INDIV_lt	:= $(addsuffix .alu.indiv.tmstp, $(OUTPUTS_l))
OUTPUTS_SVA_INDIV_lt	:= $(addsuffix .sva.indiv.tmstp, $(OUTPUTS_l))

# Simple timestamps files for group analysis.
OUTPUTS_HERVK_GRP_t		:= $(TMSTP_d)/common.bam.hervk.group.tmstp
OUTPUTS_LINE1_GRP_t		:= $(TMSTP_d)/common.bam.line1.group.tmstp
OUTPUTS_ALU_GRP_t		:= $(TMSTP_d)/common.bam.alu.group.tmstp
OUTPUTS_SVA_GRP_t		:= $(TMSTP_d)/common.bam.sva.group.tmstp

# Lists for genotype analysis.
OUTPUTS_HERVK_GEN_lt	:= $(addsuffix .hervk.gen.tmstp, $(OUTPUTS_l))
OUTPUTS_LINE1_GEN_lt	:= $(addsuffix .line1.gen.tmstp, $(OUTPUTS_l))
OUTPUTS_ALU_GEN_lt		:= $(addsuffix .alu.gen.tmstp, $(OUTPUTS_l))
OUTPUTS_SVA_GEN_lt		:= $(addsuffix .sva.gen.tmstp, $(OUTPUTS_l))

# Simple VCF files.
OUTPUTS_HERVK_VCF		:= $(HERVK_DISCOVERY_d)/HERVK.final_comp.vcf
OUTPUTS_LINE1_VCF		:= $(LINE1_DISCOVERY_d)/LINE1.final_comp.vcf
OUTPUTS_ALU_VCF			:= $(ALU_DISCOVERY_d)/ALU.final_comp.vcf
OUTPUTS_SVA_VCF			:= $(SVA_DISCOVERY_d)/SVA.final_comp.vcf


# Export all variables.
export



############################## DEBUG ###########################################

# Enable some prints, if variable DBG="yes".
ifeq ($(DBG),yes)
$(info )
$(info ****************************** DEBUG *******************************************)
$(info General information.)
$(info PIPELINE:$(PIPELINE).)
$(info SWITCH:$(SWITCH).)
$(info INPUT:$(INPUT).)
$(info INPUTS_PROCESSED:$(INPUTS_PROCESSED).)
$(info OUTPUT_DIR:$(OUTPUT_DIR).)
$(info OUTPUTS_d:$(OUTPUTS_d).)
$(info )

$(info Infrastructure.)
$(info DOCKERIZE:$(DOCKERIZE).)
$(info IMAGE:$(IMAGE).)
$(info ENCAPSULATE:$(ENCAPSULATE).)
$(info CAPSULE:$(POCKET).)
$(info INFRASTRUCTURE_dl:$(INFRASTRUCTURE_dl).)
$(info )

$(info Input files and locals.)
$(info INPUTS_l:$(INPUTS_l).)
$(info INPUTS_FILENAME_l:$(INPUTS_FILENAME_l).)
$(info INPUTS_BASENAME_l:$(INPUTS_BASENAME_l).)
$(info INPUTS_dl:$(INPUTS_dl).)
$(info INPUTS_BASENAME_dl:$(INPUTS_BASENAME_dl).)
$(info )

$(info Output files and locals.)
$(info OUTPUTS_l:$(OUTPUTS_l).)
$(info OUTPUTS_FILENAME_l:$(OUTPUTS_FILENAME_l).)
$(info OUTPUTS_BASENAME_l:$(OUTPUTS_BASENAME_l).)
$(info OUTPUTS_dl:$(OUTPUTS_dl).)
$(info OUTPUTS_BASENAME_dl:$(OUTPUTS_BASENAME_dl).)
$(info )

$(info OUTPUTS_BAI_l:$(OUTPUTS_BAI_l).)
$(info OUTPUTS_ABNORMAL_l:$(OUTPUTS_ABNORMAL_l).)
$(info OUTPUTS_HERVK_INDIV_lt:$(OUTPUTS_HERVK_INDIV_lt).)
$(info OUTPUTS_ALU_GEN_lt:$(OUTPUTS_ALU_GEN_lt).)
$(info ********************************************************************************)
$(info )
endif


# Enable debug stop, if variable STP is non-empty.
ifeq ($(STP),yes)
$(error Emergency stop)
endif



############################## MELT MAIN TARGETS ##############################

all: all_vcf
	$(info )
	$(info $(CALL) All done!)



############################## MAKE VCF #######################################

# Note: VCF are made for all BAMs together, so it can't run in parallel.

# Token target:
all_vcf: hervk_vcf line1_vcf alu_vcf sva_vcf
hervk_vcf: $(OUTPUTS_HERVK_VCF)
line1_vcf: $(OUTPUTS_LINE1_VCF)
alu_vcf: $(OUTPUTS_ALU_VCF)
sva_vcf: $(OUTPUTS_SVA_VCF)

$(OUTPUTS_HERVK_VCF): %HERVK.final_comp.vcf: $(OUTPUTS_HERVK_GEN_lt)
	$(info )
	$(info $(CALL) Create VCF for HERVK from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(HERVK_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_d) \
		-p $(HERVK_DISCOVERY_d) \
		-o $(@D)

$(OUTPUTS_LINE1_VCF): %LINE1.final_comp.vcf: $(OUTPUTS_LINE1_GEN_lt)
	$(info )
	$(info $(CALL) Create VCF for LINE1 from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(LINE1_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_d) \
		-p $(LINE1_DISCOVERY_d) \
		-o $(@D)

$(OUTPUTS_ALU_VCF): %ALU.final_comp.vcf: $(OUTPUTS_ALU_GEN_lt)
	$(info )
	$(info $(CALL) Create VCF for ALU from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(ALU_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_d) \
		-p $(ALU_DISCOVERY_d) \
		-o $(@D)

$(OUTPUTS_SVA_VCF): %SVA.final_comp.vcf: $(OUTPUTS_SVA_GEN_lt)
	$(info )
	$(info $(CALL) Create VCF for SVA from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(SVA_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_d) \
		-p $(SVA_DISCOVERY_d) \
		-o $(@D)



############################## GENOTYPING #####################################

# Token target:
all_gen: hervk_gen line1_gen alu_gen sva_gen
hervk_gen: $(OUTPUTS_HERVK_GEN_lt)
line1_gen: $(OUTPUTS_LINE1_GEN_lt)
alu_gen: $(OUTPUTS_ALU_GEN_lt)
sva_gen: $(OUTPUTS_SVA_GEN_lt)

$(OUTPUTS_HERVK_GEN_lt): %.bam.hervk.gen.tmstp: %.bam $(OUTPUTS_HERVK_GRP_t)
	$(info )
	$(info $(CALL) Genotyping HERVK in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_d) \
		-p $(HERVK_DISCOVERY_d)
	touch $@

$(OUTPUTS_LINE1_GEN_lt): %.bam.line1.gen.tmstp: %.bam $(OUTPUTS_LINE1_GRP_t)
	$(info )
	$(info $(CALL) Genotyping LINE1 in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_d) \
		-p $(LINE1_DISCOVERY_d)
	touch $@

$(OUTPUTS_ALU_GEN_lt): %.bam.alu.gen.tmstp: %.bam $(OUTPUTS_ALU_GRP_t)
	$(info )
	$(info $(CALL) Genotyping ALU in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_d) \
		-p $(ALU_DISCOVERY_d)
	touch $@

$(OUTPUTS_SVA_GEN_lt): %.bam.sva.gen.tmstp: %.bam $(OUTPUTS_SVA_GRP_t)
	$(info )
	$(info $(CALL) Genotyping SVA in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_d) \
		-p $(SVA_DISCOVERY_d)
	touch $@



############################## GROUP ANALYSES #################################

# Note: Group analysis are made for all together, so it can't run in parallel.

# Token target:
all_group: hervk_group line1_group alu_group sva_group
hervk_group: $(OUTPUTS_HERVK_GRP_t)
line1_group: $(OUTPUTS_LINE1_GRP_t)
alu_group: $(OUTPUTS_ALU_GRP_t)
sva_group: $(OUTPUTS_SVA_GRP_t)

$(OUTPUTS_HERVK_GRP_t): $(OUTPUTS_HERVK_INDIV_lt) | $(TMSTP_d)
	$(info )
	$(info $(CALL) Group analysis of HERVK on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(HERVK_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_d) \
		-n $(REFERENCE_BED)
	touch $@

$(OUTPUTS_LINE1_GRP_t): $(OUTPUTS_LINE1_INDIV_lt) | $(TMSTP_d)
	$(info )
	$(info $(CALL) Group analysis of LINE1 on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(LINE1_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_d) \
		-n $(REFERENCE_BED)
	touch $@

$(OUTPUTS_ALU_GRP_t): $(OUTPUTS_ALU_INDIV_lt) | $(TMSTP_d)
	$(info )
	$(info $(CALL) Group analysis of ALU on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(ALU_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_d) \
		-n $(REFERENCE_BED)
	touch $@

$(OUTPUTS_SVA_GRP_t): $(OUTPUTS_SVA_INDIV_lt) | $(TMSTP_d)
	$(info )
	$(info $(CALL) Group analysis of SVA on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(SVA_DISCOVERY_d) \
		-h $(REFERENCE_GENOME) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_d) \
		-n $(REFERENCE_BED)
	touch $@



############################## INDIVIDUAL ANALYSES ############################

# Token target:
all_indiv: hervk_indiv line1_indiv alu_indiv sva_indiv
hervk_indiv: $(OUTPUTS_HERVK_INDIV_lt)
line1_indiv: $(OUTPUTS_LINE1_INDIV_lt)
alu_indiv: $(OUTPUTS_ALU_INDIV_lt)
sva_indiv: $(OUTPUTS_SVA_INDIV_lt)

# Setting timestamps for several MEI types processes.
$(OUTPUTS_HERVK_INDIV_lt): %.bam.hervk.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of HERVK in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_d) \
		-c $(BAM_COVERAGE)
	touch $@

$(OUTPUTS_LINE1_INDIV_lt): %.bam.line1.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of LINE1 in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_d) \
		-c $(BAM_COVERAGE)
	touch $@

$(OUTPUTS_ALU_INDIV_lt): %.bam.alu.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of ALU in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_d) \
		-c $(BAM_COVERAGE)
	touch $@

$(OUTPUTS_SVA_INDIV_lt): %.bam.sva.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of SVA in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_d) \
		-c $(BAM_COVERAGE)
	touch $@



############################## PREPROCESSING ##################################

# Token target:
preprocess: $(OUTPUTS_DISC_l)

# Preprocessing BAMs.
$(OUTPUTS_DISC_l): %.bam.disc: %.bam %.bam.bai $(REFERENCE_GENOME)
	$(info )
	$(info $(CALL) Preprocessing file $<.)
	$(PREPROCESS) -bamfile $< -h $(REFERENCE_GENOME)


# Projects' complementary infastructure directories.
$(HERVK_d) $(LINE1_d) $(ALU_d) $(SVA_d) $(TMSTP_d):
	$(info )
	$(info $(CALL) Create directory: $@.)
	mkdir -p $@



############################## COMMON TARGETS #################################

# Lazy evaluation variable for this part.
REQ = $(filter %$(*F)$(SUFFIXES), $(INPUTS_l))

# Token target:
index: $(OUTPUTS_BAI_l) 

# Make .bai files list.
$(OUTPUTS_BAI_l): %.bam.bai: %.bam
	$(info )
	$(info $(CALL) Indexing file $^.)
	if [ -s "$$(readlink -f $<).bai" ]; then \
		ln -sf "$$(readlink -f $<).bai" $@; \
	else \
		samtools index -b -@ 8 $< $@; \
	fi


sort: $(OUTPUTS_l)

# Make .bam files list.
.SECONDEXPANSION:
$(OUTPUTS_l): %.bam: $$(REQ) | validation $(OUTPUTS_d)
	$(info )
	$(info $(CALL) Creating link for input file: $(REQ) -- $^.)
	mkdir -p $(*D)
	if [ "$$(samtools view -H $< 2> /dev/null | head -n1 | cut -f3)" = "SO:coordinate" ]; then \
		ln -sf "$$(readlink -f $<)" $@; \
	elif [ -s "$$(find $$(dirname $$(readlink -f $<)) -type f -name '*$(*F)*.sorted.bam')" ]; then \
		ln -sf "$$(find $$(dirname $$(readlink -f $<)) -type f -name '*$(*F)*.sorted.bam')" $@; \
	else \
		samtools sort -O BAM -m 4G -@ 8 $< -o $@; \
	fi


# Projects' complementary infastructure directories.
$(INFRASTRUCTURE_dl):
	$(info )
	$(info $(CALL) Create directory: $@.)
	mkdir -p $@

validation:
	$(info )
	$(info $(CALL) Making $@.)

# Simple multipurpose test.
simple_test:
	echo "It's a simple test for arguments:$(ARGS)."



############################## EXTRA TARGETS ###################################

# Search for docker.
docker_have: Dockerfile
	# Search for docker image.

# Other independent targets.
dockerize: Dockerfile
	$(info )
	$(info Make a docker image.)
	docker build -f docker/Dockerfile -t melt:latest

docker_run: Makefile $(CONFIG_f) | $(INFRASTRUCTURE_dl)
	$(info )
	$(info Run docker image with MELT command.)
	cp $(CONFIG_f) $(CONFIG_d)/$(notdir $(CONFIG_f))
	cp Makefile $(CONFIG_d)/Makefile
	$(DOCKER_RUN) \
		$(DOCKER_FLAGS) \
		-u 1541:1000 \
		-w $(C_BASE_d) \
		-v $(H_CONFIG_d):$(C_CONFIG_d) \
		-v $(H_INPUTS_d):$(C_INPUTS_d) \
		-v $(H_OUTPUTS_d):$(C_OUTPUTS_d) \
		-v $(H_ASSETS_d):$(C_ASSETS_d) \
		-v $(H_REFERENCE_d):$(C_REFERENCE_d) \
		-v $(H_ANNOTATION_d):$(C_ANNOTATION_d) \
		-v $(H_EXTRA_d):$(C_EXTRA_d) \
		-v $(H_TMP_d):$(C_TMP_d) \
		melt $(ARGS)
