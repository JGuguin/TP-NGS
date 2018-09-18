#!/bin/bash				#dit que le language est du bash
mkdir -p /mnt/data/variant_calling		#-p permet de créer les dossiers récursivement
cd /mnt/data/variant_calling

########################################################################################################################
# Requirements:
#	Java (version 8)
#	FastQC (version 0.11.7)
#	BWA-MEM (version 0.7.17-r1194-dirty)
#	SAMtools (version 1.9)
#	IGV (version 2.4.14)
#	GATK (version 3.3)
########################################################################################################################

java -version
fastqc -version
bwa		#fait l'alignement
samtools	#permet de manipuler les fichiers
java -jar ${GATK} --help	# GATK pour détection de variants
java -jar ${PICARD}			#PICARD pour duplicat de PCR

##########################################################
## Download, extract and index the reference chromosome ##
##########################################################

# Download the reference Human chromosome (chromosome 20) from Ensembl
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed reference sequence (.fa.gz)
wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz -O Homo_sapiens.Chr20.fa.gz

# Extract the reference chromosome
# Command: gunzip
# Input: compressed reference sequence (.fa.gz)
# Ouput: reference sequence (remove .gz file)
gunzip Homo_sapiens.Chr20.fa.gz

# Index the reference chromosome
# Command: bwa index
# Input: reference (.fa)
# Ouput: indexed reference (.fa.amb, .fa.ann, .fa.bwt, fa.pac, .fa.sa)
bwa index Homo_sapiens.Chr20.fa		#créer plein de fichiers qui servent d'accéder rapidement à dofférentes parties du chr

######################################################
## Mapping of a family trio to the reference genome ##
######################################################

# The sequences are from an East Asian (Kinh Vietnamese) family forming a trio : daughter/mother/father
# Data available at http://www.internationalgenome.org/data-portal/sample/HG02024
# Daughter:
#       StudyId: SRP004063
#       SampleName: HG02024
#       Library: Pond-206419
#       ExperimentID: SRX001595
#       RunId: SRR822251
#       PlatformUnit: C1E0PACXX121221.6.tagged_373
#       InstrumentModel: Illumina HiSeq 2000
#       InsertSize: 160
# Mother:
#       RUN_ID: SRP361100
#       SampleName: HG02025
# Father:
#       SampleName: HG02026

#############################
## Mapping of the daughter ##
#############################

# Download paired sequencing reads for the daughter (SampleName: HG02024, RunId: SRR822251)
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_1.filt.fastq.gz	#Première partie de la paire
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_2.filt.fastq.gz	#Seconde partie de la paire


# Map the paired sequencing reads against the reference Human chromosome 20
# Command: bwa mem
# Options: -M (Mark shorter split hits as secondary for GATK compatibility)
#          -t [number of CPU] (multi-threading)
# Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
# Ouput: alignment (.sam)
bwa mem -t 2 -M Homo_sapiens.Chr20.fa SRR822251_1.filt.fastq.gz SRR822251_2.filt.fastq.gz > HG02024_SRR822251.sam

# (Optional)
# Compute summary statistics of the alignment
# Command: samtools flagstats
# Input: alignment (.sam)
# Ouput: text file (human and computer readable)
samtools flagstat HG02024_SRR822251.sam > HG02024_SRR822251.sam.flagstats

head HG02024_SRR822251.sam	#visualiser les alignements

# Compress the alignment and filter unaligned reads		#enlève les séquences non alignées
# Command: samtools view
# Options: -@ [number of CPU] (multi-threading)
#	   -S (input format is auto-detected)
# 	   -b (output BAM)
#	   -h (include header in output)
#          -f [flag] (include reads with all  of the FLAGs in INT present)
# 	        flag=3 for read paired & read mapped in proper pair, see
#	        https://broadinstitute.github.io/picard/explain-flags.html
# Input: alignment (.sam)
# Ouput: compressed alignment (.bam)
samtools view -@ 2 -Sbh -f 3 HG02024_SRR822251.sam > HG02024_SRR822251.bam

# Sort the alignment	#trie selon l'ordre physique
# Command: samtools sort
# Input: compressed alignment (.bam)
# Ouput: sorted and compressed alignment (.bam)
samtools sort -@ 2 HG02024_SRR822251.bam > HG02024_SRR822251.sorted.bam

# Add Read group (cf https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)
# Command: picard AddOrReplaceReadGroups
# Input: alignment (.bam) and read group (Read group identifier, DNA preparation library identifier, Platform, Platform Unit, Sample)
# Ouput: annotated alignment (.bam)
java -jar ${PICARD} AddOrReplaceReadGroups I=HG02024_SRR822251.sorted.bam \
                                         O=daughter.bam \
                                         RGID=SRR822251 RGLB=Pond-206419 RGPL=illumina \
                                         RGPU=C1E0PACXX121221.6.tagged_373 RGSM=HG02024 RGPI=160

# (Optional)
# Compute statistics of the alignment
# Command: samtools-stats
# Input: alignment (.bam)
# Ouput: text file (human and computer readable)
samtools stats daughter.bam > daughter.bam.stats

# (Optional)
# Plot statistics of the alignment
# Command: plot-bamstats
# Input: statistics text file (output of samtools-stats)
# Ouput: plots (.png)
mkdir plots
plot-bamstats -p /mnt/data/plots/ daughter.bam.stats

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index daughter.bam	#revenir dans data, permet de faire des calculs rapidement, création d'un sommaire


###########################
## Mapping of the mother ##
###########################

#Trouver les infos de la mère
#Obtenir l'ID de la mère
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
grep "HG02024" 20130606_g1k.ped
head 20130606_g1k.ped

#ID father HG02026 = sample_name
#ID mother HG02025 =sample_name

#Trouver les infos complémentaires via le run_id de la mère
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.analysis.sequence.index
head n-1 20130502.phase3.analysis.sequence.index > entete.txt #crée un fichier avec juste la première ligne
grep "SRR361100" 20130502.phase3.analysis.sequence.index >> entete.txt #place ce grep à la suite de ce qui a été fait précédemment
cat entete.txt	#visualisation

#Transfert en local pour ouverture avec excel
Utilisation de filezilla

# Variables definition - on rentre tout ça dans la console et cela va permettre de définir des variables. Ainsi, quand on fera appelle aux variables plus tard, ça adaptera d'office
FTP_SEQ_FOLDER=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3 # Ftp folder from 1000Genomes project
RUN_ID=SRR361100 # Read group identifier
SAMPLE_NAME=HG02025 # Sample
INSTRUMENT_PLATFORM=Illumina # Platform/technology used to produce the read
LIBRARY_NAME=Catch-88584 # DNA preparation library identifier
RUN_NAME=BI.PE.110902_SL-HBC_0182_AFCD046MACXX.2.tagged_851.srf # Platform Unit
INSERT_SIZE=96 # Insert size

# Download paired sequencing reads for the mother
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_1.filt.fastq.gz
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_2.filt.fastq.gz

# Map, filter, and sort the paired sequencing reads of the mother against the reference genome
# Command: bwa mem && samtools view && samtools sort
# Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
# Ouput: sorted alignment (.bam)
bwa mem -t 2 -M Homo_sapiens.Chr20.fa ${RUN_ID}_1.filt.fastq.gz ${RUN_ID}_2.filt.fastq.gz | samtools view -@ 2 -Sbh -f 3 | samtools sort -@ 2 > ${SAMPLE_NAME}_${RUN_ID}.sorted.bam

# Add Read group
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group
# Ouput: alignment (.bam)
java -jar ${PICARD} AddOrReplaceReadGroups I=${SAMPLE_NAME}_${RUN_ID}.sorted.bam O=mother.bam \
                                         RGID=${RUN_ID} RGLB=${LIBRARY_NAME} RGPL=${INSTRUMENT_PLATFORM} \
                                         RGPU=${RUN_NAME} RGSM=${SAMPLE_NAME} RGPI=${INSERT_SIZE}

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index mother.bam

###########################
## Mapping of the father ##
###########################

# Variables definition
SAMPLE_NAME=HG02026 # Sample

# Download index file containing sequencing runs information
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: text file (.index)
wget ${FTP_SEQ_FOLDER}/20130502.phase3.analysis.sequence.index -O 20130502.phase3.index

# Filter paired exome sequencing runs related to father (HG02026)
# Command: grep && grep -v
# Input: tab-separated values file (.index)
# Ouput: filtered comma-separated values file (.index)
grep ${SAMPLE_NAME} 20130502.phase3.index | grep "exome" | grep 'PAIRED' | grep -v 'Catch-88526' | grep -v 'Solexa' | grep -v 'from blood' | grep -v '_1.filt.fastq.gz' | grep -v '_2.filt.fastq.gz' | sed 's/\t/,/g' > father.index

#sed remplace les tabulations par des virgules car certaines colonnes sont vides et cela permet de garder la séparation en colonne

# File containing the list of alignments (each line is a .bam file)
# This file is necessary to merge multiple alignments into a single alignment.
# Command: touch
# Input: file name
# Ouput: empty file (.bamlist)
touch father.bamlist

# for each sequencing run (the first 10), align to the reference, sort, add read group and index
head -6 father.index | while IFS="," read FASTQ_FILE MD5 RUN_ID STUDY_ID STUDY_NAME CENTER_NAME SUBMISSION_ID SUBMISSION_DATE SAMPLE_ID SAMPLE_NAME POPULATION EXPERIMENT_ID INSTRUMENT_PLATFORM INSTRUMENT_MODEL LIBRARY_NAME RUN_NAME RUN_BLOCK_NAME INSERT_SIZE LIBRARY_LAYOUT PAIRED_FASTQ WITHDRAWN WITHDRAWN_DATE COMMENT READ_COUNT BASE_COUNT ANALYSIS_GROUP
do
#head -6 peut être remplacé par cat et ça le fera pour tous les échantillons et pas seulement que pour les 6 premières lignes

    # Variables definition
    FASTQ_FILE_1=${FASTQ_FILE/.filt.fastq.gz/_1.filt.fastq.gz} # Path of the fasta file in the FTP folder
    FASTQ_FILE_2=${FASTQ_FILE/.filt.fastq.gz/_2.filt.fastq.gz} # Path of the fasta file in the FTP folder (pairing file)

    # Download paired sequencing reads for the father
    # Command: wget
    # Input: url (http:// or ftp://)
    # Ouput: compressed sequencing reads (.fastq.gz)
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_1} -O ${SAMPLE_NAME}_${RUN_ID}_1.filt.fastq.gz
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_2} -O ${SAMPLE_NAME}_${RUN_ID}_2.filt.fastq.gz

    # Map, filter, and sort the paired reads of the sequencing run against the reference genome
    # Command: bwa mem && samtools view && samtools sort
    # Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
    # Ouput: sorted alignment (.bam)
    bwa mem -t 2 -M Homo_sapiens.Chr20.fa ${SAMPLE_NAME}_${RUN_ID}_1.filt.fastq.gz ${SAMPLE_NAME}_${RUN_ID}_2.filt.fastq.gz | samtools view -@ 2 -Sbh -f 3 | samtools sort -@ 2 > ${SAMPLE_NAME}_${RUN_ID}.sorted.bam

    # Add Read group
    # Command: gatk AddOrReplaceReadGroups
    # Input: alignment (.bam) and read group
    # Ouput: alignment (.bam)
    java -jar ${PICARD} AddOrReplaceReadGroups I=${SAMPLE_NAME}_${RUN_ID}.sorted.bam O=${SAMPLE_NAME}_${RUN_ID}.sorted.RG.bam \
                                         RGID=${RUN_ID} RGLB=${LIBRARY_NAME} RGPL=${INSTRUMENT_PLATFORM} \
                                         RGPU=${RUN_NAME} RGSM=${SAMPLE_NAME} RGPI=${INSERT_SIZE}

    # Index the alignment
    # Command: samtools index
    # Input: alignment (.bam)
    # Ouput: indexed alignment (.bam.bai)
    samtools index ${SAMPLE_NAME}_${RUN_ID}.sorted.RG.bam

    # Append the file name (.bam) to the list of alignments that will be merged
    echo ${SAMPLE_NAME}_${RUN_ID}.sorted.RG.bam >> father.bamlist
done

# Merge the list of alignments into a single file
# Command: samtools merge
# Ouput: alignment (.bam)
# Input: file containing the list of alignments (each line is a .bam file)
samtools merge father.bam -b father.bamlist

# Index the alignment
# Command: samtools index
# Input: alignment (.sam or .bam)
# Ouput: indexed alignment (.sam.bai or .bam.bai)
samtools index father.bam

#Alignement des fichiers .bam de la fille, de la mère et du père sur IGV.
#Comparaison avec le génome de référence Human hg38
#Ne pas oublier de transférer avec les fichiers .bai nécessaires à la lecture des .bam
