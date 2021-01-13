Bootstrap: docker
From: ubuntu:18.04
%setup


%files
    /nfs/home/software/packages/fastp_v0.20.1/ /nfs/home/software/packages/fastp_v0.20.1
    /nfs/home/software/packages/FastQC_v0.11.8/ /nfs/home/software/packages/FastQC_v0.11.8
    /home/pnatividad/BRC1-ONC1-0834-S-A2_S7_L001_R1_001.fastq.gz /data/UMPSING_BRC1/rawdata/
    /home/pnatividad/BRC1-ONC1-0834-S-A2_S7_L001_R2_001.fastq.gz /data/UMPSING_BRC1/rawdata/
    /nfs/home/panel_designs/BRC1/capture.bed /nfs/home/panel_designs/BRC1/
    /nfs/home/databases/adapters/adapters.fa /nfs/home/databases/adapters/
    /media/scratch3/docker/docker-fastp/inc_pipe.conf /data/UMPSING_BRC1/scripts/
    /media/scratch3/docker/docker-fastp/inc_pipe_singularity.sh /data/UMPSING_BRC1/scripts/    
    /nfs/home/software/packages/bedtoolsv2.28.0-33/bin/bedtools /nfs/home/software/packages/bedtoolsv2.28.0-33/bin/
    /nfs/home/databases/GATK/small_exac_common_3.hg38.vcf.gz /nfs/home/databases/GATK/small_exac_common_3.hg38.vcf.gz
    /home/pnatividad/UMPSING_BRC1/analysis_BRC1_250x/QC/rawdata/ /data/UMPSING_BRC1/analysis_BRC1_250x/QC/rawdata/

%post
    apt-get update && apt-get install -y tabix perl zip default-jre tree
    mkdir -p /data/UMPSING_BRC1/analysis_BRC1_250x/QC/rawdata
    chmod -R 777 /data/UMPSING_BRC1/

%test
    grep -q NAME=\"Ubuntu\" /etc/os-release
    if [ $? -eq 0 ]; then
        echo "Container base is Ubuntu as expected."
    else
        echo "Container base is not Ubuntu."
    fi

%runscript
    exec "$@"

%labels
    Author pnatividad@incliva.es
    Version v0.0.1

%help
    Este contenedor es una prueba
