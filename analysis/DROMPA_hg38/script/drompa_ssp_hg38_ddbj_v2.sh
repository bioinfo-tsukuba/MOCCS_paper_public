#!/bin/bash
#$ -S /bin/bash

PROJECT_PATH=/home/s-tahara/DROMPA
srxfile=/home/s-tahara/DROMPA/code_hg38_2/SRX_TO_SRR_hg19.txt
DROMPA_PATH=/home/s-tahara/DROMPA
CHIP_PATH=/home/s-tahara/chipseq

# srx作成
declare -a srx_array=()
row_N=`cat ${srxfile} | wc -l`

for row_num in `seq 1 ${row_N}`
do
 line=`head -n ${row_num} ${srxfile} | tail -n 1`
 SRX=$(echo ${line} | awk '{print $1}')
 srx_array+=(${SRX})
done

declare -a srx_array_uniq=()
srx_array_uniq=(`echo ${srx_array[@]} | tr ' ' '\n' | uniq | xargs` )


# SRXに対応するSRRを取得
SRX_TO_SRR_target=`grep ${srx_array_uniq[((${SGE_TASK_ID}-1))]} ${srxfile}`
mkdir -p ./num_txt
echo ${SRX_TO_SRR_target}  | tr ' ' '\n' > ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num.txt
grep -v ${srx_array_uniq[((${SGE_TASK_ID}-1))]} ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num.txt > ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num2.txt 
 cat ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num2.txt
SRR_N=`wc -l ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num2.txt | awk '{print $1}'`

############################################################################################
### ここから先は、DROMPA Plusの出力がある場合はスキップする (jobを途中で止めた時のため) ####
###########################################################################################
DROMPA_out=${DROMPA_PATH}/DROMPA_PLUS_2/hg38/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${srx_array_uniq[((${SGE_TASK_ID}-1))]}*ReadCountDist.tsv

if [ -f $DROMPA_out ]; then
	echo ${srx_array_uniq[((${SGE_TASK_ID}-1))]} >> /home/s-tahara/DROMPA/code_hg38_2/skip_ID_hg19.txt
else	
	# singularity execでfastq-dumpを実行
	#target_SRR=`echo ${SRX_TO_SRR_target}  | tr ' ' '\n' | sed -e "/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/d"`
	#python -m venv cwltoolenv
	. cwltoolenv/bin/activate
	
	for SRR_num in `seq 1 ${SRR_N}`
	do
		echo ${target_SRR}
 		target_SRR=`head -n ${SRR_num} ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num2.txt | tail -n 1`
		singularity -d exec sra-tools_2.11.0.sif fasterq-dump ${target_SRR} --threads 8 --skip-technical --split-files --split-spot --outdir /home/s-tahara/DROMPA/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}
	
	done
	
	rm ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num.txt
	rm ./num_txt/${srx_array_uniq[((${SGE_TASK_ID}-1))]}_num2.txt

	### cwl workflow で得られない場合、ddbjのfastqを直接取ってくる
	mkdir -p /home/s-tahara/DROMPA/FASTQ_CON
	concatenated_fastq=${PROJECT_PATH}/FASTQ_CON/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.fastq
	fastq_file=/home/s-tahara/DROMPA/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/*.fastq
	echo ls /home/s-tahara/DROMPA/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}
	ls /home/s-tahara/DROMPA/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/

	# fastq fileが複数ある場合に備えて、一つ目のファイルをfastq_file_path2 に格納して、ifの条件文にいれる
	fastq_file_first=`ls -l /home/s-tahara/DROMPA/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]} | head -n 2 | awk '{print $9}'`
	fastq_file_first_path=/home/s-tahara/DROMPA/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/${fastq_file_first}
	fastq_file_first_path2=`echo ${fastq_file_first_path} | tr -d ' '`
	echo ${fastq_file_first_path2}

	if [ -f $fastq_file_first_path2 ] ; then

		### SRXごとにfastqを連結 ###
		echo downloaded throuth cwltool
 		cat ${fastq_file} >> ${concatenated_fastq}
 		gzip ${concatenated_fastq}

	else

 		### ddbjから直接fastqをダウンロード ###
 		echo donwloaded from ddbj database as fastq format
 		target_ID_list_path=/home/s-tahara/DROMPA/code_hg38_2/SRX_SRA_hg19_3.tab
 		line=`grep ${srx_array_uniq[((${SGE_TASK_ID}-1))]} ${target_ID_list_path}`
 
  		index2=`echo ${line} | awk '{print $2}'`
  		index1=`echo ${index2:0:-3}`
  		index3=`echo ${line} | awk '{print $1}'`

  		#prepare directory per ID
  		wget -P ${DROMPA_PATH}/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]} ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/${index1}/${index2}/${index3}/*fastq*

 		### SRXごとにfastqを連結 ###
 		bzip2 -d ${DROMPA_PATH}/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/*.fastq.bz2
	 	fastq_file=${PROJECT_PATH}/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/*.fastq
 		cat ${fastq_file} >> ${concatenated_fastq}
 		gzip ${concatenated_fastq}
 		rm ${DROMPA_PATH}/FASTQ/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/*.fastq.bz2
	fi

	rm ${fastq_file}

	### mapping bowtie2 (fastq -> SAM) ###
	echo bowtie2 start
	mkdir -p ${DROMPA_PATH}/SAM/
	bowtie2 -p 8 --seed 0 -x /home/s-tahara/bowtie2_index/GRCh38_noalt_as/GRCh38_noalt_as -U ${concatenated_fastq}.gz > ${DROMPA_PATH}/SAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.trim.sam
	rm ${concatenated_fastq}.gz
	rm ${concatenated_fastq}

### samtools (SAM -> BAM) ### 
	echo samtools start
	mkdir -p ${DROMPA_PATH}/BAM
	samtools view -bhS -F 0x4 -q 42 ${DROMPA_PATH}/SAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.trim.sam | samtools sort -T ${DROMPA_PATH}/BAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.trim > ${DROMPA_PATH}/BAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.trim.uniq.bam
	rm ${DROMPA_PATH}/SAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.trim.sam
	rm ${DROMPA_PATH}/SAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.trim
	
	### DROMPA SSP ###
	echo SSP start
	mkdir -p ${DROMPA_PATH}/SSP_2/
	mkdir -p ${DROMPA_PATH}/SSP_2/hg38/
	mkdir -p ${DROMPA_PATH}/SSP_2/hg38/${srx_array_uniq[((${SGE_TASK_ID}-1))]}
	singularity run lolcow.sif 
	singularity exec ssp_drompa.img ssp -i ${DROMPA_PATH}/BAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}*.bam -o ${srx_array_uniq[((${SGE_TASK_ID}-1))]} --gt ${CHIP_PATH}/genometable_hg38.txt --verbose -p 4 --odir ${DROMPA_PATH}/SSP_2/hg38/${srx_array_uniq[((${SGE_TASK_ID}-1))]} 

	### DROMPA Plus (Library complexity & GC content) ###
	echo DROMPA Plus start
	mkdir -p ${DROMPA_PATH}/DROMPA_PLUS_2
	mkdir -p ${DROMPA_PATH}/DROMPA_PLUS_2/hg38
	mkdir -p ${DROMPA_PATH}/DROMPA_PLUS_2/hg38/${srx_array_uniq[((${SGE_TASK_ID}-1))]}
	singularity exec ssp_drompa.img parse2wig+ -i ${DROMPA_PATH}/BAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}*.bam -o ${srx_array_uniq[((${SGE_TASK_ID}-1))]} --gt ${CHIP_PATH}/genometable_hg38.txt --verbose -p 4 --odir ${DROMPA_PATH}/DROMPA_PLUS_2/hg38/${srx_array_uniq[((${SGE_TASK_ID}-1))]} --chrdir /home/s-tahara/chipseq/chrdir/hg38 

	rm ${DROMPA_PATH}/BAM/${srx_array_uniq[((${SGE_TASK_ID}-1))]}.trim.uniq.bam
	
	# DROMPA outputの使わないファイルを消す
	rm ${DROMPA_PATH}/DROMPA_PLUS_2/hg38/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/*.bw
	rm ${DROMPA_PATH}/DROMPA_PLUS_2/hg38/${srx_array_uniq[((${SGE_TASK_ID}-1))]}/*.R

fi # skip 用 ifの終わり
