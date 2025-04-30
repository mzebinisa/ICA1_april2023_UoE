#!/bin/bash

# Asking the user for the input directory
echo "Please enter the path to the directory containing the input files, e.g. /localdisk/data/BPSM/ICA1/fastq/"
read INPUT_DIR

# Telling the user about copying the input files to the current directory
echo "Dear USER, this directory contains some other files, to which we don't have access. Therefore we will copy the ones we need into the current directory. It will be a bit messy now, but it will help us reduce other messiness further. Thanks for understanding and collaboration."

# Copying the input files to the current directory
cp ${INPUT_DIR}/Tco-{5053..6950}_{1,2}.fq.gz .

# Prompting the user for the output directory
echo "Dear user, please enter the directory where you would like to store the output files:"
read OUTPUT_DIR

# Going through the input files and run fastqc
echo "Running FastQC on input files..."
for file in Tco-*.fq.gz; do
  fastqc --extract "${file}" -o "$OUTPUT_DIR"
done
echo "FastQC completed."

# Generating a summary report for all FastQC results
echo "Generating a summary report for FastQC results..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR"
echo "MultiQC completed."

# Informing the user that FastQC reports are saved in the output directory
echo "FastQC reports are saved in the "$OUTPUT_DIR" directory. Please analyze them before further processing."

# Asking the user for the path to the reference genome FASTA file
echo "Please enter the path to the reference genome FASTA file, e.g. /localdisk/home/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz"
read REF_GENOME

# Again asking the user for the desired name for the Bowtie2 index
echo "Please enter a name for the Bowtie2 index, e.g. TcongolenseIL3000_index:"
read INDEX_NAME

# Generating Bowtie2 index using the reference genome FASTA file
echo "Building Bowtie2 index using reference genome FASTA file..."
bowtie2-build $REF_GENOME $INDEX_NAME
echo "Reference index is created."

# Looping through all paired-end read files in the directory and align them to the reference genome
echo "Aligning paired-end read files to the reference genome..."
for file in Tco-*_1.fq.gz; do
  # Extract the common part of the filename
  base=$(basename $file _1.fq.gz)

  # Align the paired-end reads to the reference genome using Bowtie2
  bowtie2 -x $INDEX_NAME -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz -S ${base}.sam
  
  # Turnign the SAM file to BAM file and sort it
  samtools view -bS ${base}.sam | samtools sort -o ${base}.sorted.bam
  
  # Indexing the sorted BAM file
  samtools index ${base}.sorted.bam
  
  # Removing the intermediate SAM file
  rm ${base}.sam
done
echo "Alignment completed."

# Informing the user that the alignment is complete
echo "Alignment of fastq files to ${INDEX_NAME} is complete."
echo "*.sam files have been generated for each sam file."

# Prompting the user for the path to the bedfile
echo "Please enter the full path and filename to the bedfile, e.g. /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed:"
read BED_DIR

# Creating an array of sorted BAM file paths in sampleID order
bam_files=($(ls -1v *.sorted.bam))

# Looping over the sorted BAM files in sampleID order
for file in "${bam_files[@]}"; do
  # Eliminating  files that end in ".sorted.bam.bai"
  if [[ $file == *.sorted.bam.bai ]]; then
     continue
   fi

  # Getting the sample name
   sample=$(basename "$file" .sorted.bam)

  # bedtools intersect to identify overlaps with the genomic regions
  echo "Finding overlaps between ${sample}.sorted.bam and ${BED_DIR}..."
  bedtools intersect -a "$BED_DIR" -b "$file" -bed > "$sample.overlaps.bed"

  # and bedtools multicov to generate counts data
  echo "Generating counts data for ${sample}.sorted.bam..."
  bedtools multicov -bams "$file" -bed "$BED_DIR" > "$sample.counts.txt"
done

# Concatenating all  count data  into a common count data file with all the samples ordered
echo "Concatenating all the count data files into a common count data file with all the samples ordered..."
cat *.counts.txt | sort -k1,1 > all_counts.txt
echo "All the bam files are processed with bedtools intersect first and then with multicov. Results are saved into individual *.bed and *.count.txt files as well as one common all_counts.txt file."

# Copying the file containing the fqfile information to the current directory
echo "Copying the file containing the fqfile information to the current directory..."
cp /localdisk/data/BPSM/ICA1/fastq/Tco.fqfiles .

# Modifying the fqfile to include a hyphen in the SampleName field
echo "Modifying the fqfile to include a hyphen in the SampleName field..."
sed -i 's/^\(Tco\|\tTco\)/&-/' Tco.fqfiles

# Printing the first 10 lines of the fqfile
echo "The first 10 lines of the modified fqfile are:"
head -n 10 Tco.fqfiles

# Printing the the details we have
echo "The content of fastq details file is as follows:"
cat Tco.fqfiles

# Setting the variable fqfile to the fastq details file path
fqfile="Tco.fqfiles"

# Counting the number of replicates for each group and save the results to a file
echo "Counting the number of replicates for each group and saving the results to reps.groups file..."
while read SampleName SampleType Replicate Time Treatment End1 End2; do
    echo "${SampleType}.${Time}.${Treatment}"
done < "${fqfile}" | tail -n +2 | sort | uniq -c > reps.groups
cat reps.groups

# Scanning each line in the fqfile and extract the group name and fqfile names
echo "Looping through each line in the fqfile and extracting the group name and fqfile names..."
while read SampleName SampleType Replicate Time Treatment End1 End2; do
    # Constructing the group name using the SampleType, Time, and Treatment fields
    group_name="${SampleType}.${Time}.${Treatment}"
    # Creating a directory for the group (if it doesn't already exist)
    mkdir -p "${group_name}"

    # Finding the corresponding count file for this fqfile
    count_file="${SampleName}.counts.txt"

    # Copying the count file to the group directory
    cp "${count_file}" "${group_name}/"
    
    # Combining the count files for this group into a single file
    paste "${group_name}"/*.counts.txt > "${group_name}/${group_name}.genecounts" 
done < "${fqfile}"
echo "Group name and fqfile names extracted, count files copied to corresponding group directories, and genecounts files created."

# Looping through each group and calculate the statistical means
echo "Looping through each group and calculating the statistical means..."
while read reps group; do
    # Calculating the statistical means for each row in the input file
    awk -v reps="${reps}" '{
        sum = 0;
        # Start at the relevant column and sum to the end of each row
        for (i=NF-reps+1; i<=NF; i++) sum+=$i;
        # Print the gene name, chromosome, start position, end position, and mean expression
        print $5,$(5+1),$(5+2),$NF, sum/reps;
    }' "${group}/${group}.genecounts_cleaned" > "${group}/${group}_stmean.txt"
    echo "Statistical means calculated for ${group}."
done < reps.groups

# Printing a message to indicate that the statistical means have been calculated and saved
echo "The statistical means for each group are generated and saved in corresponding files: ${group}_stmean.txt."