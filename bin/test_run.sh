


#Good Run on Google Cloud

ngseasy_initiate_project -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects

ngseasy_volumes_container-d /media/container-vol/ngs_projects

ngseasy_initiate_fastq -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects

# get 
sudo docker run --name annovar_db_get --volumes_from data_volumes -i -t compbio/ngseasy-annovar:v0.9 /bin/bash 


## To Test



