############################## HEADER #########################################
#
# Makefile data:
#
# 		Default goal:	analyse_all
# 		File:			Makefile
# 		Module:			none
#		Project:		melt
# 		Author:			Jose L. L. Buzzo
# 		Organization:	RetroTeam
# 		Year:			2018
#
#
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
#
###############################################################################





############################## PREAMBLE #######################################

# Include user's configuration file (only the first one).
SWITCH			:= $(shell [ -f .ponga_switch ] && echo 1)
CONFIG_FILE 	?= $(if $(SWITCH),/home/config/config.mk,config.mk)
aux 			:= $(word 1, $(abspath $(strip $(wildcard $(CONFIG_FILE)))))
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

OUTPUT_DIR := $(abspath $(strip $(OUTPUT_DIR)))# This is not a string list!

# Folder for timestamps.
TMSTP := $(OUTPUT_DIR)/tmstp


# Carrefully set the INPUT_l variable after success validations.
INPUT_l := $(INPUT_PROCESSED)
INPUT_FILENAME_l := $(strip $(notdir $(INPUT_l)))
INPUT_BASENAME_l := $(strip $(basename $(notdir $(INPUT_l))))
INPUT_DIR_l := $(strip $(dir $(INPUT_l)))
#INPUT_DIR_BASENAME_l :=


#OUTPUT_l :=
#OUTPUT_FILENAME_l :=
#OUTPUT_BASENAME_l :=
#OUTOUT_DIR_l :=
OUTPUT_DIR_BASENAME_l := $(addprefix $(OUTPUT_DIR)/, $(INPUT_BASENAME_l))
OUTPUT_l := $(strip $(join $(OUTPUT_DIR_BASENAME_l), $(addprefix /, $(INPUT_FILENAME_l))))



############################## LISTS ##########################################

# Derivated lists: '.bai', '.abnormal', '.disc', '.disc.bai' and others.
OUTPUT_BAI_l := $(addsuffix .bai, $(OUTPUT_l))
OUTPUT_ABNORMAL_l := $(addsuffix .abnormal, $(OUTPUT_l))
OUTPUT_DISC_l := $(addsuffix .disc, $(OUTPUT_l))
OUTPUT_DISC_BAI_l := $(addsuffix .disc.bai, $(OUTPUT_l))
OUTPUT_DISC_FQ_l:= $(addsuffix .disc.fq, $(OUTPUT_l))

# Timestamp lists, terminated in '_tml'.
# Lists for individual analysis.
OUTPUT_HERVK_INDIV_tml := $(addsuffix .hervk.indiv.tmstp, $(OUTPUT_l))
OUTPUT_LINE1_INDIV_tml := $(addsuffix .line1.indiv.tmstp, $(OUTPUT_l))
OUTPUT_ALU_INDIV_tml := $(addsuffix .alu.indiv.tmstp, $(OUTPUT_l))
OUTPUT_SVA_INDIV_tml := $(addsuffix .sva.indiv.tmstp, $(OUTPUT_l))

# Simple timestamps files for group analysis.
OUTPUT_HERVK_GRP := $(TMSTP)/common.bam.hervk.group.tmstp
OUTPUT_LINE1_GRP := $(TMSTP)/common.bam.line1.group.tmstp
OUTPUT_ALU_GRP := $(TMSTP)/common.bam.alu.group.tmstp
OUTPUT_SVA_GRP := $(TMSTP)/common.bam.sva.group.tmstp

# Lists for genotype analysis.
OUTPUT_HERVK_GEN_tml := $(addsuffix .hervk.gen.tmstp, $(OUTPUT_l))
OUTPUT_LINE1_GEN_tml := $(addsuffix .line1.gen.tmstp, $(OUTPUT_l))
OUTPUT_ALU_GEN_tml := $(addsuffix .alu.gen.tmstp, $(OUTPUT_l))
OUTPUT_SVA_GEN_tml := $(addsuffix .sva.gen.tmstp, $(OUTPUT_l))

# Simple VCF files.
OUTPUT_HERVK_VCF := $(HERVK_DISCOVERY_DIR)/HERVK.final_comp.vcf
OUTPUT_LINE1_VCF := $(LINE1_DISCOVERY_DIR)/LINE1.final_comp.vcf
OUTPUT_ALU_VCF := $(ALU_DISCOVERY_DIR)/ALU.final_comp.vcf
OUTPUT_SVA_VCF := $(SVA_DISCOVERY_DIR)/SVA.final_comp.vcf


# Export all variables.
export



# ===> DEBUG CODE <===

# Enable some prints, if variable DBG="yes".
ifeq ($(DBG),yes)
$(info )
$(info OUTPUT_DIR: >$(OUTPUT_DIR).)
$(info )
$(info INPUT_l: >$(INPUT_l).)
$(info INPUT_FILENAME_l: >$(INPUT_FILENAME_l).)
$(info INPUT_BASENAME_l: >$(INPUT_BASENAME_l).)
$(info INPUT_DIR_l: >$(INPUT_DIR_l).)
$(info INPUT_DIR_BASENAME_l: >$(INPUT_DIR_BASENAME_l).)
$(info )
$(info OUTPUT_l: >$(OUTPUT_l).)
$(info OUTPUT_FILENAME_l: >$(OUTPUT_FILENAME_l).)
$(info OUTPUT_BASENAME_l: >$(OUTPUT_BASENAME_l).)
$(info OUTPUT_DIR_l: >$(OUTPUT_DIR_l).)
$(info OUTPUT_DIR_BASENAME_l: >$(OUTPUT_DIR_BASENAME_l).)
$(info )
$(info OUTPUT_ABNORMAL_l: >$(OUTPUT_ABNORMAL_l).)
$(info OUTPUT_HERVK_INDIV_tml: >$(OUTPUT_HERVK_INDIV_tml).)
$(info OUTPUT_ALU_GEN_tml: >$(OUTPUT_ALU_GEN_tml).)
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



############################## LASY EVALUATION VARIABLES ######################

# Look for a *.sorted.bam.bai file, according to the link to 'OUTPUT_l'.
#REQ_BAI = $(shell readlink -f $(filter %$(*F)$(SUFFIXES), $(OUTPUT_l)))
#aux2 = $(strip $(filter %.sorted.bam, $(REQ_BAI)))
#STBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*.bai' -regex '.*sorted.*')
#NSTBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*.bai' -not -regex '.*sorted.*')
#BAI = $(if $(aux2),$(STBAI),$(NSTBAI))

# Ancilary variables to this target.
REQ = $(filter %$(*F)$(SUFFIXES), $(INPUT_l))
#CMD = $(shell samtools view -H $(REQ) 2> /dev/null | head -n1 | cut -f3)
#BAM = $(shell readlink -f $(REQ))
#STBAM = $(shell find $(dir $(BAM)) -type f -name '*$(*F)*.sorted.bam')

############################## MELT MAIN TARGETS ##############################

all: all_vcf
	$(info )
	$(info $(CALL) All done!)




############################## MAKE VCF #######################################

# Note: VCF are made for all BAMs together, so it can't run in parallel.

# Token target:
all_vcf: hervk_vcf line1_vcf alu_vcf sva_vcf
hervk_vcf: $(OUTPUT_HERVK_VCF)
line1_vcf: $(OUTPUT_LINE1_VCF)
alu_vcf: $(OUTPUT_ALU_VCF)
sva_vcf: $(OUTPUT_SVA_VCF)

$(OUTPUT_HERVK_VCF): %HERVK.final_comp.vcf: $(OUTPUT_HERVK_GEN_tml)
	$(info )
	$(info $(CALL) Create VCF for HERVK from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(HERVK_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_DIR) \
		-p $(HERVK_DISCOVERY_DIR) \
		-o $(@D)

$(OUTPUT_LINE1_VCF): %LINE1.final_comp.vcf: $(OUTPUT_LINE1_GEN_tml)
	$(info )
	$(info $(CALL) Create VCF for LINE1 from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(LINE1_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_DIR) \
		-p $(LINE1_DISCOVERY_DIR) \
		-o $(@D)

$(OUTPUT_ALU_VCF): %ALU.final_comp.vcf: $(OUTPUT_ALU_GEN_tml)
	$(info )
	$(info $(CALL) Create VCF for ALU from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(ALU_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_DIR) \
		-p $(ALU_DISCOVERY_DIR) \
		-o $(@D)

$(OUTPUT_SVA_VCF): %SVA.final_comp.vcf: $(OUTPUT_SVA_GEN_tml)
	$(info )
	$(info $(CALL) Create VCF for SVA from files $^.)
	$(MAKE_VCF) \
		-genotypingdir $(SVA_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_DIR) \
		-p $(SVA_DISCOVERY_DIR) \
		-o $(@D)



############################## GENOTYPING #####################################

# Token target:
all_gen: hervk_gen line1_gen alu_gen sva_gen
hervk_gen: $(OUTPUT_HERVK_GEN_tml)
line1_gen: $(OUTPUT_LINE1_GEN_tml)
alu_gen: $(OUTPUT_ALU_GEN_tml)
sva_gen: $(OUTPUT_SVA_GEN_tml)

$(OUTPUT_HERVK_GEN_tml): %.bam.hervk.gen.tmstp: %.bam $(OUTPUT_HERVK_GRP)
	$(info )
	$(info $(CALL) Genotyping HERVK in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_DIR) \
		-p $(HERVK_DISCOVERY_DIR)
	touch $@

$(OUTPUT_LINE1_GEN_tml): %.bam.line1.gen.tmstp: %.bam $(OUTPUT_LINE1_GRP)
	$(info )
	$(info $(CALL) Genotyping LINE1 in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_DIR) \
		-p $(LINE1_DISCOVERY_DIR)
	touch $@

$(OUTPUT_ALU_GEN_tml): %.bam.alu.gen.tmstp: %.bam $(OUTPUT_ALU_GRP)
	$(info )
	$(info $(CALL) Genotyping ALU in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_DIR) \
		-p $(ALU_DISCOVERY_DIR)
	touch $@

$(OUTPUT_SVA_GEN_tml): %.bam.sva.gen.tmstp: %.bam $(OUTPUT_SVA_GRP)
	$(info )
	$(info $(CALL) Genotyping SVA in file $<.)
	$(GENOTYPE) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_DIR) \
		-p $(SVA_DISCOVERY_DIR)
	touch $@



############################## GROUP ANALYSES #################################

# Note: Group analysis are made for all together, so it can't run in parallel.

# Token target:
all_group: hervk_group line1_group alu_group sva_group
hervk_group: $(OUTPUT_HERVK_GRP)
line1_group: $(OUTPUT_LINE1_GRP)
alu_group: $(OUTPUT_ALU_GRP)
sva_group: $(OUTPUT_SVA_GRP)

$(OUTPUT_HERVK_GRP): $(OUTPUT_HERVK_INDIV_tml) | $(TMSTP)
	$(info )
	$(info $(CALL) Group analysis of HERVK on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(HERVK_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_DIR) \
		-n $(REFERENCE_BED)
	touch $@

$(OUTPUT_LINE1_GRP): $(OUTPUT_LINE1_INDIV_tml) | $(TMSTP)
	$(info )
	$(info $(CALL) Group analysis of LINE1 on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(LINE1_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_DIR) \
		-n $(REFERENCE_BED)
	touch $@

$(OUTPUT_ALU_GRP): $(OUTPUT_ALU_INDIV_tml) | $(TMSTP)
	$(info )
	$(info $(CALL) Group analysis of ALU on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(ALU_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_DIR) \
		-n $(REFERENCE_BED)
	touch $@

$(OUTPUT_SVA_GRP): $(OUTPUT_SVA_INDIV_tml) | $(TMSTP)
	$(info )
	$(info $(CALL) Group analysis of SVA on files $^.)
	$(GROUP_ANALYSIS) \
		-discoverydir $(SVA_DISCOVERY_DIR) \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_DIR) \
		-n $(REFERENCE_BED)
	touch $@



############################## INDIVIDUAL ANALYSES ############################

# Token target:
all_indiv: hervk_indiv line1_indiv alu_indiv sva_indiv
hervk_indiv: $(OUTPUT_HERVK_INDIV_tml)
line1_indiv: $(OUTPUT_LINE1_INDIV_tml)
alu_indiv: $(OUTPUT_ALU_INDIV_tml)
sva_indiv: $(OUTPUT_SVA_INDIV_tml)

# Setting timestamps for several MEI types processes.
$(OUTPUT_HERVK_INDIV_tml): %.bam.hervk.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of HERVK in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(HERVK_ZIP) \
		-w $(HERVK_DISCOVERY_DIR) \
		-c $(BAM_COVERAGE)
	touch $@

$(OUTPUT_LINE1_INDIV_tml): %.bam.line1.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of LINE1 in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(LINE1_ZIP) \
		-w $(LINE1_DISCOVERY_DIR) \
		-c $(BAM_COVERAGE)
	touch $@

$(OUTPUT_ALU_INDIV_tml): %.bam.alu.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of ALU in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(ALU_ZIP) \
		-w $(ALU_DISCOVERY_DIR) \
		-c $(BAM_COVERAGE)
	touch $@

$(OUTPUT_SVA_INDIV_tml): %.bam.sva.indiv.tmstp: %.bam %.bam.disc
	$(info )
	$(info $(CALL) Individual analysis of SVA in file $<.)
	$(INDIV_ANALYSIS) \
		-bamfile $< \
		-h $(REFERENCE_GENOME_FASTA) \
		-t $(SVA_ZIP) \
		-w $(SVA_DISCOVERY_DIR) \
		-c $(BAM_COVERAGE)
	touch $@



############################## PREPROCESSING ##################################

# Token target:
preprocess: $(OUTPUT_DISC_l)

# Preprocessing BAMs.
$(OUTPUT_DISC_l): %.bam.disc: %.bam %.bam.bai
	$(info )
	$(info $(CALL) Preprocessing file $<.)
	$(PREPROCESS) -bamfile $< -h $(REFERENCE_GENOME_FASTA)



############################## COMMON TARGETS #################################

# Token target:
common: $(OUTPUT_BAI_l) 

$(OUTPUT_BAI_l): %.bam.bai: %.bam
	$(info )
	$(info $(CALL) Indexing file $^.)
	if [ -s "$$(readlink -f $<).bai" ]; then \
		ln -sf "$$(readlink -f $<).bai" $@; \
	else \
		samtools index -b -@ 8 $< $@; \
	fi

.SECONDEXPANSION:
$(OUTPUT_l): %.bam: $$(REQ) | validation $(OUTPUT_DIR)
	$(info )
	$(info $(CALL) Creating link for input file: $(REQ).)
	$(info $(CALL) SEE: $@ | $< | $^.)
	mkdir -p $(*D)
	if [ "$$(samtools view -H $< 2> /dev/null | head -n1 | cut -f3)" = "SO:coordinate" ]; then \
		ln -sf "$$(readlink -f $<)" $@; \
	elif [ -s "$$(find $$(dirname $$(readlink -f $<)) -type f -name '*$(*F)*.sorted.bam')" ]; then \
		ln -sf "$$(find $$(dirname $$(readlink -f $<)) -type f -name '*$(*F)*.sorted.bam')" $@; \
	else \
		samtools sort -O BAM -m 8G -@ 8 $< -o $@; \
	fi


# Directories.
$(OUTPUT_DIR) $(HERVK_DIR) $(LINE1_DIR) $(ALU_DIR) $(SVA_DIR) $(TMSTP):
	$(info )
	$(info $(CALL) Create directory: $@.)
	mkdir -p $@


validation:
	$(info )
	$(info $(CALL) Making $@.)



############################## OTHER TARGETS ##################################

# Simple test.
simple_test:
	echo "It's a simple test for arguments: $(ARGS)."

# Infrastructure directories.
infra_dirs: $(HOST_CONFIG) $(HOST_INPUTS) $(HOST_OUTPUTS) $(HOST_ASSETS)
infra_dirs: $(HOST_REFERENCE) $(HOST_ANNOTATION) $(HOST_EXTRA) $(HOST_TMP)

# Search for docker.
docker_have: Dockerfile
	# Search for docker image.

# Other independent targets.
dockerize: Dockerfile
	$(info )
	$(info Make a docker image.)
	docker build -t melt:latest .

docker_run: Makefile $(CONFIG_FILE) | infra_dirs
	$(info )
	$(info Run docker image with MELT command.)
	cp $(CONFIG_FILE) $(CONFIG)/$(notdir $(CONFIG_FILE))
	docker run \
		--rm \
		-ti \
		-u 1541:1000 \
		-v $(HOST_CONFIG):$(CONFIG) \
		-v $(HOST_INPUTS):$(INPUTS) \
		-v $(HOST_OUTPUTS):$(OUTPUTS) \
		-v $(HOST_ASSETS):$(ASSETS) \
		-v $(HOST_REFERENCE):$(REFERENCE) \
		-v $(HOST_ANNOTATION):$(ANNOTATION) \
		-v $(HOST_EXTRA):$(EXTRA) \
		-v $(HOST_TMP):$(TMP) \
		melt $(ARGS)
