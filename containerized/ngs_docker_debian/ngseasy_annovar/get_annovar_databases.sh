#!/bin/bash

#Get ANNOVAR DATA BASES
for annovardatabase in \
  ljb26_phylop100way_vertebrate \
  ljb26_all \
  cosmic70 \
  esp6500si_all \
  1000g2014sep \
  snp138 \
  nci60 \
  clinvar_20140929 \
  gerp++elem \
  targetScanS \
  tfbsConsSites \
  wgRna;
do 
  ./annotate_variation.pl --buildver hg19 --downdb --webfrom annovar ${annovardatabase} humandb/ ; 
done

# [ljb26_all] whole-exome SIFT scores, PolyPhen2 HDIV scores, PolyPhen2 HVAR scores, LRT scores, MutationTaster scores, MutationAssessor score, FATHMM scores, MetaSVM scores, MetaLR scores, VEST scores, CADD scores, GERP++ scores, PhyloP scores and SiPhy scores from dbsnp version 2.6

