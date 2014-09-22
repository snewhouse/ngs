INSTALL_TARGET=${1}  # /usr/local/pi
PATH_FILE=${2}  # .bashrc

# + Trimmomatic
wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip \
    && unzip /tmp/Trimmomatic-0.32.zip -d $INSTALL_TARGET/ \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:$INSTALL_TARGET/Trimmomatic-0.32/trimmomatic-0.32.jar' $PATH_FILE \
    && sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/Trimmomatic-0.32' $PATH_FILE

# + FastQC
wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip \
    && unzip /tmp/fastqc_v0.11.2.zip -d $INSTALL_TARGET/ \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:$INSTALL_TARGET/FastQC/jbzip2-0.9.jar:$INSTALL_TARGET/FastQC/sam-1.103.jar' $PATH_FILE \
    && sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/FastQC' $PATH_FILE

# + bowtie
wget -O /tmp/bowtie2-2.2.3-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip \
    && unzip /tmp/bowtie2-2.2.3-linux-x86_64.zip -d $INSTALL_TARGET/ \
    && sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/bowtie2-2.2.3:$INSTALL_TARGET/bowtie2-2.2.3/scripts' $PATH_FILE \
    && echo "alias ngsBowtie2='$INSTALL_TARGET/bowtie2-2.2.3'" >>  $PATH_FILE

# + bwa
wget -O /tmp/bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2 \
    && tar xjvf /tmp/bwa-0.7.10.tar.bz2 -C $INSTALL_TARGET/ \
    && cd $INSTALL_TARGET/bwa-0.7.10 && make \
    && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/bwa-0.7.10' $PATH_FILE \
    && echo "alias ngsBWA='$INSTALL_TARGET/bwa-0.7.10'" >>  $PATH_FILE

# + stampy (registration required, get compressed file and put in context dir for the build)
wget -O /tmp/stampy-latest.tgz  http://www.well.ox.ac.uk/~gerton/software/Stampy/stampy-latest.tgz \
    && tar xvf /tmp/stampy-latest.tgz -C $INSTALL_TARGET/ \
    && sed -i 's/-Wl//' $INSTALL_TARGET/stampy-1.0.23/makefile \
    && chmod -R 755 $INSTALL_TARGET/stampy-1.0.23 \
    && cd $INSTALL_TARGET/stampy-1.0.23 && make \
    && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/stampy-1.0.23' $PATH_FILE \
    && echo "alias ngsStampy='$INSTALL_TARGET/stampy-1.0.23'" >>  $PATH_FILE

# + VCFtools: http://vcftools.sourceforge.net/index.html
wget -O /tmp/vcftools_0.1.12a.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12a.tar.gz/download \
    && tar xzvf /tmp/vcftools_0.1.12a.tar.gz -C $INSTALL_TARGET/  \
    && cd $INSTALL_TARGET/vcftools_0.1.12a/ && make \
    && sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/vcftools_0.1.12a/bin' $PATH_FILE

# + GATK (see licence, and registration)
wget --no-check-certificate -O /tmp/GenomeAnalysisTK-3.2-2.tar.bz2 https://www.dropbox.com/s/wey1edv9tqdc7so/GenomeAnalysisTK-3.2-2.tar.bz2 \
  && mkdir $INSTALL_TARGET/GenomeAnalysisTK-3.2-2 \
  && tar xjvf /tmp/GenomeAnalysisTK-3.2-2.tar.bz2 -C $INSTALL_TARGET/GenomeAnalysisTK-3.2-2 \
  && chmod -R 755 $INSTALL_TARGET/GenomeAnalysisTK-3.2-2 \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:$INSTALL_TARGET/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar' $PATH_FILE \
  && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/GenomeAnalysisTK-3.2-2' $PATH_FILE \
  && echo "alias ngsGATK='$INSTALL_TARGET/GenomeAnalysisTK-3.2-2'" >>  $PATH_FILE

# + GATK resources
wget --no-check-certificate -O /tmp/gatk_resources.tar.gz https://www.dropbox.com/s/ff1b5l2xbhisidx/gatk_resources.tar.gz \
 && tar xvf /tmp/gatk_resources.tar.gz -C $INSTALL_TARGET/ \
 && gunzip $INSTALL_TARGET/gatk_resources/*.gz \
 && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/gatk_resources' $PATH_FILE \
 && echo "alias ngsGATKresources='$INSTALL_TARGET/gatk_resources'" >>  $PATH_FILE

# + samtools, htslib and bcftools
cd $INSTALL_TARGET \
	&& git clone --branch=develop git://github.com/samtools/htslib.git \
	&& git clone --branch=develop git://github.com/samtools/bcftools.git \
	&& git clone --branch=develop git://github.com/samtools/samtools.git \
	&& cd $INSTALL_TARGET/bcftools \
	&& make \
	&& cd $INSTALL_TARGET/samtools \
	&& make \
	&& cd $INSTALL_TARGET/htslib \
	&& make \
	&& sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/samtools' $PATH_FILE \
  && sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/bcftools' $PATH_FILE \
  && sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/htslib' $PATH_FILE \

# + SAMtools
wget -O /tmp/samtools-0.1.19.tar.bz2 http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download \
    && tar xjvf /tmp/samtools-0.1.19.tar.bz2 -C $INSTALL_TARGET/ \
    && cd $INSTALL_TARGET/samtools-0.1.19 && make \
    && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/samtools-0.1.19/bin' $PATH_FILE

# + Picard
wget -O /tmp/picard-tools-1.115.zip http://sourceforge.net/projects/picard/files/picard-tools/1.115/picard-tools-1.115.zip \
    && mkdir $INSTALL_TARGET/picardtools \
    && unzip /tmp/picard-tools-1.115.zip -d $INSTALL_TARGET/picardtools/ \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:$INSTALL_TARGET/picardtools/picard-tools-1.115/snappy-java-1.0.3-rc3.jar' $PATH_FILE \
    && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/picardtools/picard-tools-1.115' $PATH_FILE \
    && echo "alias ngsPicard='$INSTALL_TARGET/picardtools/picard-tools-1.115'" >>  $PATH_FILE

# + BEDtools
cd $INSTALL_TARGET \
    && git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2 && make clean && make all \
    && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/bedtools2/bin' $PATH_FILE \
    && echo "alias ngsBedtools='$INSTALL_TARGET/bedtools2/bin'" >> $PATH_FILE

# + VCFtools: http://vcftools.sourceforge.net/index.html
wget -O /tmp/vcftools_0.1.12a.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12a.tar.gz \
    && tar xzvf /tmp/vcftools_0.1.12a.tar.gz -C $INSTALL_TARGET/  \
    && cd $INSTALL_TARGET/vcftools_0.1.12a/ && make \
    && sed  -i '$aPATH=${PATH}:$INSTALL_TARGET/vcftools_0.1.12a/bin' $PATH_FILE \
    && echo "alias ngsVCFtools='$INSTALL_TARGET/vcftools_0.1.12a/bin/'" >>  $PATH_FILE

# + SnpEff
wget -O /tmp/snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    && unzip /tmp/snpEff_latest_core.zip -d $INSTALL_TARGET/ \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:$INSTALL_TARGET/snpEff/snpEff.jar' $PATH_FILE \
    && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/snpEff' $PATH_FILE \
    && echo "alias ngsSNPeff='$INSTALL_TARGET/snpEff'" >> $PATH_FILE

# + ANNOVAR (see licence, and registration)
wget -O /tmp/annovar.latest.tar.gz http://www.openbioinformatics.org/annovar/download/mP628pfL21/annovar.latest.tar.gz \
  && tar xzvf /tmp/annovar.latest.tar.gz -C $INSTALL_TARGET/ \
  && cd $INSTALL_TARGET/annovar \
  && cp -v $INSTALL_TARGET/annovar/table_annovar.pl table_annovar.pl.original \
  && wget http://www.openbioinformatics.org/annovar/download/table_annovar.pl \
  && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/annovar' $PATH_FILE \
   && echo "alias ngsSNPeff='$INSTALL_TARGET/annovar'" >> $PATH_FILE

# + ANNOVAR databases
mkdir $INSTALL_TARGET/annovar \
  && cd $INSTALL_TARGET/annovar \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar snp138 humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar ljb23_all humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar cosmic70 humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar clinvar_20140902 humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar targetScanS humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar wgRna humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar tfbsConsSites humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar knownGene humandb \
  && $INSTALL_TARGET/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar ensGene humandb \

# + freebayes
cd $INSTALL_TARGET \
  && git clone --recursive git://github.com/ekg/freebayes.git \
  && cd $INSTALL_TARGET/freebayes \
  && make \
  && chmod -R 755 $INSTALL_TARGET/freebayes \
  && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/freebayes/bin' $PATH_FILE \
  && echo "alias ngsFreebayes='$INSTALL_TARGET/freebayes/bin'" >>  $PATH_FILE

# + Platypus
wget -O /tmp/Platypus-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz \
  && tar xvf /tmp/Platypus-latest.tgz -C $INSTALL_TARGET/ \
  && cd $INSTALL_TARGET/Platypus_0.7.9.1 \
  && sh ./buildPlatypus.sh \
  && chmod -R 755 $INSTALL_TARGET/Platypus_0.7.9.1 \
  && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/Platypus_0.7.9.1' $PATH_FILE \
  && echo "alias ngsPlatypus='$INSTALL_TARGET/Platypus_0.7.9.1'" >>  $PATH_FILE

# + SAVANT
wget -O /tmp/Savant-full-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Savant-full-latest.tgz \
  && tar xvf /tmp/Savant-full-latest.tgz -C $INSTALL_TARGET/ \
  && cd $INSTALL_TARGET/savant-v1.1.0-complete \
  && chmod -R 755 $INSTALL_TARGET/savant-v1.1.0-complete \
  && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/savant-v1.1.0-complete' $PATH_FILE \
  && echo "alias ngsSavant='$INSTALL_TARGET/savant-v1.1.0-complete'" >>  $PATH_FILE


