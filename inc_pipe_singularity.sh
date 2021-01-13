#!/bin/bash

function checkDirectories(){
        for dir in ${@}
        do
            if [ ! -d "${dir}" ]
            then
                mkdir -p ${dir}
            fi
        done
}


function fastq_preprocessing(){
    echo "Running fastq preprocessing"
    date
    echo "${fastqc} -o ${5}/QC/rawdata -t ${6} ${1} ${2}"
    ${fastqc} -o ${5}/QC/rawdata -t ${6} ${1} ${2}
    date
    if [ "${7}" == "yes" ]
    then
        echo "Processing UMIs"
        date
        echo "$umi_tools extract --bc-pattern=NNNNNNNNN --stdin=${8} --read2-in=${1} --stdout=${4}/QC/${sample}_umi_R1.fastq.gz --read2-stdout"
	echo "$umi_tools extract --bc-pattern=NNNNNNNNN --stdin=${8} --read2-in=${2} --stdout=${4}/QC/${sample}_umi_R2.fastq.gz --read2-stdout"
	$umi_tools extract --bc-pattern=NNNNNNNNN --stdin=${8} --read2-in=${1} --stdout=${4}/QC/${sample}_umi_R1.fastq.gz --read2-stdout
        $umi_tools extract --bc-pattern=NNNNNNNNN --stdin=${8} --read2-in=${2} --stdout=${4}/QC/${sample}_umi_R2.fastq.gz --read2-stdout
        R1=`echo ${4}/QC/${sample}_umi_R1.fastq.gz`
        R2=`echo ${4}/QC/${sample}_umi_R2.fastq.gz`
    else
        echo "Sample does not have UMIs"
        R1=${1}
        R2=${2}
    fi

    date
    $fastp -i ${R1} -I ${R2} -o ${4}/QC/filtered_${3}_1.fastq.gz -O ${4}/QC/filtered_${3}_2.fastq.gz --adapter_fasta $adapterList -g -x --average_qual 20 --length_required $min_read_length -p -h ${4}/QC/filtered_${3}_report.html -w ${threads} --cut_right --cut_mean_quality 30 --poly_x_min_len 5 --low_complexity_filter --poly_g_min_len 10 --poly_x_min_len 10 -j ${4}/QC/filtered_${3}.json
    ${fastqc} -o ${5}/QC/rawdata -t ${6} ${4}/QC/filtered_${3}_1.fastq.gz ${4}/QC/filtered_${3}_2.fastq.gz
    date
}


function usage() {
    echo "Usage: $0"
    echo "This script runs the pipeline for target seq data"
    echo "Mandatory parameters:"
    echo "-f configuration file"
    1>&2; exit 1;
}

while getopts "f:" opt; do
    case ${opt} in
        f)
                f=${OPTARG} ;;
        *)
                usage ;;
    esac
done

if [ -z "$f" ]
then
        echo "ERROR: -f is a mandatory parameter"
        usage
        exit
fi

date
echo "Your command: "$@

source ${f}

checkDirectories ${tmp_dir}/mapping
checkDirectories ${tmp_dir}/variant_calling
checkDirectories ${tmp_dir}/QC
checkDirectories ${analysis_dir}/mapping
checkDirectories ${analysis_dir}/variant_calling
checkDirectories ${analysis_dir}/QC/mapping
checkDirectories ${analysis_dir}/QC/rawdata
mkdir -p ${tmp_dir}/variant_calling
mkdir -p ${tmp_dir}/QC
mkdir -p ${analysis_dir}/mapping
mkdir -p ${analysis_dir}/variant_calling
mkdir -p ${analysis_dir}/QC/mapping
mkdir -p ${analysis_dir}/QC/rawdata

#ln -s -f ${INC_somatic} ${panel}/
#ln -s -f ${INC_somatic}.tbi ${panel}/
#ln -s -f ${vep_db} ${panel}/
#ln -s -f ${dbSNP} ${panel}/
#ln -s -f ${dbSNP}.tbi ${panel}/
#ln -s -f ${GATK_indels} ${panel}/
#ln -s -f ${GATK_indels}.tbi ${panel}/
#ln -s -f ${COSMIC} ${panel}/
#ln -s -f ${COSMIC}.tbi ${panel}/
#ln -s -f ${HOTSPOTS} ${panel}/
#ln -s -f ${HOTSPOTS}.tbi ${panel}/
#ln -s -f ${INC_somatic_unzip} ${panel}/
#ln -s -f ${MIXED} ${panel}/
#ln -s -f ${MIXED}.tbi ${panel}/
${bedtools} intersect -a ${GATK_small} -b ${panel}/capture.bed -header > ${panel}/exac_roche.vcf
bgzip -f ${panel}/exac_roche.vcf && tabix -p vcf -f ${panel}/exac_roche.vcf.gz


echo "Sample,Rawdata,%LowQReads,%MappedReads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_100x,Nt_250x,Nt_500x,Mean_insert,Insert_SD,NumVars" > ${analysis_dir}/QC/mapping/${name}_global_mapping_stats.csv
if [ ${panel} != "/nfs/home/panel_designs/TP53" ]
    then
        echo "------NO AMPLICONS------"
        echo "NO ES TP53"
        p="no"
    else
        echo "Sample,Rawdata,%LowQReads,%MappedReads,%DuplicateReads,%OnTargetNoDupReads,Kit_specificity,Cov_1stQ,Cov_3rdQ,Cov_Median,Nt_100x,Nt_300x,Nt_500x,Mean_insert,Insert_SD,NumVars" > ${analysis_dir}/QC/mapping/${name}_global_mapping_stats.csv
        echo "ES TP53"
        p=`echo ${samples} | awk '{split($0,a,","); print a[1]}'`
    fi
if [ ${BRCA} == "yes" ]
    then
        p=`echo ${samples} | awk '{split($0,a,","); print a[1]}'`
    fi
sample_list=`echo ${samples} | sed 's/,/ /gi'`
cwd=`pwd`

for sample in ${sample_list[@]}
do
    echo "Running analysis pipeline on sample ${sample}"
    R1=`ls ${rawdata}/$sample*_R1*fast*`
    R2=`ls ${rawdata}/$sample*_R2*fast*`
    if [ ${has_umi} == "no" ]
    then
	echo "------$sample DOES NOT HAVE UMIs------"
    else
	echo "------$sample HAS UMIs------"
	R1=`ls ${rawdata}/$sample*_R1*fast*`
    	R2=`ls ${rawdata}/$sample*_R3*fast*`
	umi=`ls ${rawdata}/$sample*_R2*fast*`
    fi
    if [ ${amplicon} == "yes" ]
    then
        echo "------DESIGN FOR SAMPLE $sample IS AMPLICON-BASED------"
        #amplicon="yes"
    else
        echo "------DESIGN FOR SAMPLE $sample IS NOT AMPLICON-BASED------"
        #amplicon="no"
    fi

    fastq_preprocessing ${R1} ${R2} ${sample} ${tmp_dir} ${analysis_dir} ${threads} ${has_umi} ${umi}

done
