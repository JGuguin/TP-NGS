#!/bin/bash

# Oupppsss, I'm so empty ... I want to be filled to actually perform the analysis

cd /mnt/data

#Polymorphism dataset

#Téléchargement du dossier avec tout
wget https://tinyurl.com/2018-tp-polymorphism
unzip 2018-tp-polymorphism

cd 2018-tp-polymorphism-part2/
mkdir 1000genomes
cd 1000genomes/

#Téléchargement des données .vcf
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr20_GRCh38.genotypes.20170504.vcf.gz
gunzip ALL.chr20_GRCh38.genotypes.20170504.vcf


#Filtering out non-coding polymorphism

#Téléchargement de la référence
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.94.chr.gtf.gz

#Filtre les intervalles bien annotés avec début et fin d'exons
gtf_to_bed.py -g Homo_sapiens.GRCh38.94.chr.gtf     #[-h] optionnel

#Garder les exons en un seul exemplaire
bedtools sort -i Homo_sapiens.GRCh38.94.chr.gtf > Homo_sorted.bed
bedtools merge -c 4 -o distinct -i Homo_sorted.bed > Homo_merge.bed
bedtools intersect -header -wb -a ALL.chr20_GRCh38.genotypes.20170504.vcf -b Homo_merge.bed > Chr20_intersect.vcf
#-wb garde la localisation pour quel transcrit on a fait l'intersection



