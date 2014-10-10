#!/bin/bash

#ANNOVAR DATA BASES
./annotate_variation.pl --buildver hg19 --downdb seq humandb/hg19_seq

./annotate_variation.pl --buildver hg19 --downdb --webfrom annovar refGene  humandb/

./annotate_variation.pl --buildver hg19 --downdb --webfrom annovar knownGene  humandb/

./annotate_variation.pl --buildver hg19 --downdb --webfrom annovar ensGene  humandb/

./retrieve_seq_from_fasta.pl humandb/hg19_refGene.txt -seqdir humandb/hg19_seq -format refGene -outfile humandb/hg19_refGeneMrna.fa

./retrieve_seq_from_fasta.pl humandb/hg19_knownGene.txt -seqdir humandb/hg19_seq -format knownGene -outfile humandb/hg19_knownGeneMrna.fa

./retrieve_seq_from_fasta.pl humandb/hg19_ensGene.txt -seqdir humandb/hg19_seq -format ensGene -outfile humandb/hg19_ensGeneMrna.fa
