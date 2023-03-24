#RNA-seq Reads QC, Trimming, and Mapping via Kallisto
#Arguments: ./rnaseq_pipe.sh [input reads directory] [reference genome location]

#!/bin/bash

#load necessary modules
module load fastqc
module load trimmomatic
module load kallisto
module load sambamba

#load args
HOMEPATH="$1"
REF="$2"

#go to data
cd $HOMEPATH
mkdir $HOMEPATH/abundance
mkdir $HOMEPATH/sorted_bam/

#set up output directories in read directory
mkdir fastqc_out
mkdir kal_out
mkdir kal_out/counts
mkdir kal_out/bam

#QC and Trim fastq files
for SUBDIR in *; do
	mkdir $HOMEPATH/$SUBDIR/trimmed
	cd $SUBDIR
	for FILE in *; do
		if [ "${FILE: -8}" == "fastq.gz" ]; then
			fastqc $FILE -o $HOMEPATH/fastqc_out
			java -jar /ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar SE $FILE $HOMEPATH/$SUBDIR/trimmed/$FILE CROP:50 HEADCROP:2 MINLEN:40 AVGQUAL:30
		fi
	done
	cd ..
done

#Map Reads via Kallisto
kallisto index -i Lp_ref $REF
for SUBDIR in *; do
	cd $SUBDIR/trimmed
	mkdir kal_map	
	kallisto quant -i $HOMEPATH/Lp_cds_ref -o kal_map *R1_001.fastq.gz *R2_001.fastq.gz --pseudobam
	cd kal_map
	sambamba sort *.bam
	find * -maxdepth 0 -exec mv {} "$SUBDIR"_{} \;
	mv *abundance* $HOMEPATH/abundance_lvbR/
	mv *sorted* $HOMEPATH/sorted_bam_lvbR/
	cd $HOMEPATH
done

rm Lp_ref
rm trimmed
rm kal_map
