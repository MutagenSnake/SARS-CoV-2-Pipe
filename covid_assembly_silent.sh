#!/bin/bash

#Adding path to the input folder
input_folder="$1"
#Refering to previously created reference genome file
ref_genome="reference_files/Reference_EPI_ISL_402124.fasta"
#Refering to preciously created reference genome file
bed_file="reference_files/primers.bed"
#Refering to previously created nextclade reference folder
nextclade_ref_folder="reference_files/nextclade_reference"
# Get the input folder name without the path
input_folder_name=$(basename "$input_folder")
#Creating output folder
mkdir -p "${input_folder_name}_output"
#Refering to output folder
output_folder="${input_folder_name}_output"
#Creating temp folder, for storing intermediate pipeline data
mkdir -p temp

find $input_folder -type f -exec mv {} $input_folder \;

toilet "NVSPL SARS-CoV-2 Version: 0.11" -F gay -f smblock

# Creating a while loop to iterate through files (including the nested folders)
for file in "$input_folder"/*; do
    # Checking if the reads are paired-end Illumina reads
    if [[ $file == *"_R1_"* ]]; then
        # Assigning the R1 file
        R1_file="$file"
        # Assigning the R2 file
        R2_file="${R1_file/_R1_/_R2_}"
        # Assigning sample name
        sample_name="${R1_file##*/}"    # Remove the path
        sample_name="${sample_name/_R1_*/}"  # Remove _R1_ and everything after
		toilet ${sample_name} -f term -F border --gay
        # Mapping both pairs to the reference genome
        bwa mem "$ref_genome" "$R1_file" "$R2_file" -o "temp/${sample_name}.sam" > /dev/null 2>&1
		echo "bwa mem complete"
        # Generating the .bam file and sorting it
        samtools view -b "temp/${sample_name}.sam" | samtools sort -o "temp/${sample_name}.bam" > /dev/null 2>&1
		echo "samtools view complete"
        # Trimming primers from sorted .bam file
        ivar trim -e -i "temp/${sample_name}.bam" -b "$bed_file" -m 30 -q 20 -p "temp/${sample_name}.trimmed.bam" > /dev/null 2>&1
		echo "ivar trim complete"
        # Sorting the trimmed .bam file
        samtools sort "temp/${sample_name}.trimmed.bam" -o "temp/${sample_name}.trimmed.sorted.bam" > /dev/null 2>&1
		echo "samtools sort complete"
        # Generating a .vcf file with gene variations between reads and reference genome
        bcftools mpileup -B -d 250 --max-idepth 1000 --annotate INFO/AD,FORMAT/AD -Q 30 -f $ref_genome temp/"${sample_name}.trimmed.sorted.bam" | 
		bcftools call -Ou -m | 
		bcftools +fill-tags -- -t FORMAT/VAF | 
		bcftools +setGT -- -t q -i 'GT="1/1" && FORMAT/VAF < 0.8' -n 'c:0/1' |
		bcftools +setGT -- -t q -i 'GT="0/1" && FORMAT/VAF >= 0.8' -n 'c:1/1' |
		bcftools filter -o temp/"${sample_name}.vcf" -e 'INFO/IMF < 0.5' -- > /dev/null 2>&1
		echo "mpileup complete"
        # Compressing the variant calling file
        bgzip "temp/${sample_name}.vcf" > /dev/null 2>&1 # Automatically generates a .vcf.gz file
		echo "bgzip complete"
        # Indexing the compressed variant calling file
        tabix -p vcf "temp/${sample_name}.vcf.gz" > /dev/null 2>&1 # Automatically generates a .vcf.gz.tbi file
		echo "Tabix complete"
        # Assemble the genome between indexed variant calling and the reference genome
        bcftools consensus -a 'N' -p "$sample_name" -f "$ref_genome" -H I -i 'INFO/DP >= 10' "temp/${sample_name}.vcf.gz" > "temp/${sample_name}.fasta"
		echo "bcftools consensus complete"
		#Removing temporary files
		rm temp/${sample_name}.sam
		rm temp/${sample_name}.bam
		rm temp/${sample_name}.bam.bai
		rm temp/${sample_name}.trimmed.bam
		rm temp/${sample_name}.trimmed.sorted.bam
		rm temp/${sample_name}.vcf.gz
		rm temp/${sample_name}.vcf.gz.tbi
    fi
done


#Conacatenating all the fasta files to one
cat temp/*.fasta > ${output_folder}/${input_folder_name}_consensus.fasta
#Running nextclade to find a lineage
nextclade run -D $nextclade_ref_folder $output_folder/${input_folder_name}_consensus.fasta --output-tsv $output_folder/${input_folder_name}_lineage.tsv
#Running pangolin to find a lineage
conda activate pangolin
pangolin $nextclade_ref_folder $output_folder/${input_folder_name}_consensus.fasta -o $output_folder/ --outfile ${input_folder_name}_pango_lineage.csv

toilet "Complete" -F gay -f smblock
