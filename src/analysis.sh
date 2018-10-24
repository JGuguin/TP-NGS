#!/bin/bash
cd /mnt/data/
mkdir 1000genomes
cd 1000genomes/

# Téléchargement des données

# Génotype des individus (.vcf)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr20_GRCh38.genotypes.20170504.vcf.gz
gunzip ALL.chr20_GRCh38.genotypes.20170504.vcf.gz

# Téléchargement de la référence
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.94.chr.gtf.gz

# Filtre les intervalles bien annotés avec début et fin d'exons
gtf_to_bed.py -g Homo_sapiens.GRCh38.94.chr.gtf
# [-h] optionnel

# Garder les exons en un seul exemplaire
bedtools sort -i Homo_sapiens.GRCh38.94.chr.bed > Homo_sorted.bed
bedtools merge -c 4 -o distinct -i Homo_sorted.bed > Homo_merge.bed
bedtools intersect -header -wb -a ALL.chr20_GRCh38.genotypes.20170504.vcf -b Homo_merge.bed > Chr20_intersect.vcf
# -wb garde la localisation pour quel transcrit on a fait l'intersection

# Classification des mutations synonymes,non-synonymes ou stop
#Téléchargement des transcrits
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gunzip Homo_sapiens.GRCh38.cds.all.fa.gz
vcf_coding_polymorphism.py -f Homo_sapiens.GRCh38.cds.all.fa -g Homo_sapiens.GRCh38.94.chr.gtf -v Chr20_intersect.vcf

# Filtrage par population + analyse de la variance
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
# Fichier qui contient la population pour chaque individu
for CLASS in "Stop" "Syn" "NonSyn"
do
    for POP in "EUR" "AFR" "EAS" "AMR" "SAS" "GBR" "FIN" "CHS" "PUR" "CDX" "JPT" "CLM" "IBS" "PEL" "PJL" "KHV" "LWK" "ACB" "GWD" "ESN" "BEB" "MSL" "MXL" "STU" "ITU" "CEU" "YRI" "CHB" "ASW" "TSI" "GIH"
    do
    extract_pop.py -p integrated_call_samples_v3.20130502.ALL.panel -k ${POP}
    vcftools --vcf Chr20_intersect.${CLASS}.vcf --keep ${POP}.txt --recode --out ${CLASS}_${POP}
    done
done

for POP in "EUR" "AFR" "EAS" "AMR" "SAS" "GBR" "FIN" "CHS" "PUR" "CDX" "JPT" "CLM" "IBS" "PEL" "PJL" "KHV" "LWK" "ACB" "GWD" "ESN" "BEB" "MSL" "MXL" "STU" "ITU" "CEU" "YRI" "CHB" "ASW" "TSI" "GIH"
    do
    vcf_meta_analysis.py -o Stop_${POP}.recode.vcf -s Syn_${POP}.recode.vcf -n NonSyn_${POP}.recode.vcf -c 200
done