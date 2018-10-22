#!/bin/bash

# Oupppsss, I'm so empty ... I want to be filled to actually perform the analysis

cd /mnt/data

#Téléchargement du dossier avec tout
wget https://tinyurl.com/2018-tp-polymorphism
unzip 2018-tp-polymorphism

cd 2018-tp-polymorphism-part2/
mkdir 1000genomes
cd 1000genomes/

#Téléchargement des données .vcf
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
gunzip ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz



