


#Good Run on Google Cloud

ngseasy_initiate_project -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects

ngseasy_volumes_container -d /media/container-vol/ngs_projects

ngseasy_initiate_fastq -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects

# get annovar databases
sudo docker run \
	--rm=true \
	--name annovar_db_get \
	--volumes-from volumes_container \
	-i \
	-t \
	compbio/ngseasy-annovar:v0.9.2 \
	/bin/bash \
	/usr/local/pipeline/annovar/get_annovar_gene_databases.sh

sudo docker run \
	--rm=true \
	--volumes-from volumes_container \
	-i \
	-t \
	compbio/ngseasy-annovar:v0.9.2 \
	/bin/bash \
	/usr/local/pipeline/annovar/get_annovar_databases.sh

## To Test

# PIPELINE SET TO ngseasy_fastq in config file
# this creates new sample confog file from project confif file and feeds as input to pipeline ngseasy_fastq
ngseasy -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects





