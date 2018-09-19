#On donne au programme l'information que les trois individus sont apparentés pour voir les différences de traitement

#!/bin/bash
WORK_DIR=/mnt/data/variant_calling/trio
cd ${WORK_DIR}

########################################################################################################################
# Requirements:
#   Java (version 8)
#   GATK (version 3.3)
########################################################################################################################

java -version
java -jar ${GATK} --help
java -jar ${PICARD}

############################
### Joint variant calling ##
############################

#Permet de grouper en un seul fichier tous les variants calling

REF_GENOME=Homo_sapiens.Chr20.fa
PEDIGREE=${WORK_DIR}/20130606_g1k.ped
FILE_NAME1=HG02024
FILE_NAME2=HG02025
FILE_NAME3=HG02026

# Perform joint variant calling
# Command: gatk GenotypeGVCFs
# Input : list of genomic variant calling files (.g.vcf) + reference genome (.fa)
# Output: Variant calling file (.vcf)
java -jar ${GATK} -T GenotypeGVCFs \
	-R ${REF_GENOME} \
	-V ${FILE_NAME1}.g.vcf \
	-V ${FILE_NAME2}.g.vcf \
	-V ${FILE_NAME3}.g.vcf \
	-o trio.vcf


##########################
### Recalibrate variants #
##########################

# Not possible here because exomes
# => Skip this part


########################################
### Post-processing: Analysis of trios #
########################################

#Indique que les trois individus sont apparentés

# Modify PEDIGREE to keep only first 6 columns and no header
cut -f 1-6 ${PEDIGREE} |sed '1,1d' > ${PEDIGREE}.txt

# Phase trios
# Command: gatk PhaseByTransmission
# Input: Variant calling file (.vcf) + reference genome (.fa) + PEDIGREE file (.ped)
# Output: Phased variant calling file (.vcf)
java -jar ${GATK} -T PhaseByTransmission \
	-R ${REF_GENOME} \
	--variant trio.vcf \
	-ped ${PEDIGREE}.txt \
	-o trio.phased.vcf

# Evaluate variants by computing control metrics
# Command: gatk VariantEval
# Input: List of variant calling files to evaluate (.vcf) + reference genome (.fa)
# Output: Variant evaluation file (.txt)
java -jar ${GATK} -T VariantEval \
	-R ${REF_GENOME} \
	--eval:set1 trio.vcf \
	--eval:set2 trio.phased.vcf \
	-o trio.phased.VE.txt

# Tabulate the number of sites which overlap and share alleles
# Command: gatk GenotypeConcordance
# Input: 2 variant callings files (.vcf) + reference genome (.fa)
# Output: Report file (.txt)
java -jar ${GATK} -T GenotypeConcordance \
	-R ${REF_GENOME} \
	--eval trio.vcf \
	--comp trio.phased.vcf \
	-o report.trio.txt
