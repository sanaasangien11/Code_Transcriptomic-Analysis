import pandas as pd
from tqdm import tqdm
import os


#if not os.path.exists(f"arabidopsis_data/hisat2_index/arabidopsis.8.ht2"):
# chroms = ",".join([f"arabidopsis_data/CHROM/{x}" for x in os.listdir(f"arabidopsis_data/CHROM")])
# os.system(f"hisat2-build -p 7 {chroms} arabidopsis_data/hisat2_index/arabidopsis")



raw_data_path = "/nam-99/ablage/nam/ss_rnaseq/wheat_mapping/raw_data"
outputs_folder_path = "/nam-99/ablage/nam/ss_rnaseq/wheat_mapping/output_folder"
cleanup_folder_path = "/nam-99/ablage/nam/ss_rnaseq/wheat_mapping/cleanup_folder"
processed_folder_path = "/nam-99/ablage/nam/ss_rnaseq/wheat_mapping/processed_data"
uncompressed_folder_path = "/nam-99/ablage/nam/ss_rnaseq/wheat_mapping/uncompressed_data"


if not os.path.exists(uncompressed_folder_path):
    os.mkdir(uncompressed_folder_path)


if not os.path.exists(processed_folder_path):
    os.mkdir(processed_folder_path)


if not os.path.exists(outputs_folder_path):
    os.mkdir(outputs_folder_path)


if not os.path.exists(cleanup_folder_path):
    os.mkdir(cleanup_folder_path)


# List all files in the raw data path
files = os.listdir(raw_data_path)

# Sort files to ensure they are processed in order


# Iterate over pairs of files (assuming the files are named in pairs as per your example)
for i in tqdm(range(0, len(files), 2)):

    # Construct the full paths to the input files
    file_path_forward = os.path.join(raw_data_path, files[i])
    file_path_reverse = os.path.join(raw_data_path, files[i + 1])

    print(f"file_path_forward is: {file_path_forward}")
    print(f"file_path_reverse is: {file_path_reverse}")


    # Extract the base names without extension for both forward and reverse reads
    file_name_forward, _ = os.path.splitext(files[i])
    file_name_reverse, _ = os.path.splitext(files[i + 1])

    print(f"file_name_forward is: {file_name_forward}")
    print(f"file_name_reverse is: {file_name_reverse}")

    #Construct the output file paths in the uncompressed folder
    output_file_forward = os.path.join(uncompressed_folder_path, f"{file_name_forward}")
    output_file_reverse = os.path.join(uncompressed_folder_path, f"{file_name_reverse}")

# Unzip the input files
# os.system(f"gzip -d {file_path_forward}")
# os.system(f"gzip -d {file_path_reverse}")

# Update file paths to point to the unzipped files
# file_path_forward = os.path.join(raw_data_path, f"{file_name_forward}")
# file_path_reverse = os.path.join(raw_data_path, f"{file_name_reverse}")

    # Unzip the input files directly to the uncompressed folder
    os.system(f"gzip -d {file_path_forward} -c > {output_file_forward}")
    os.system(f"gzip -d {file_path_reverse} -c > {output_file_reverse}")

    # Update file paths to point to the uncompressed files
    file_path_forward = output_file_forward
    file_path_reverse = output_file_reverse

    print(f"file_path_forward updated to: {file_path_forward}")
    print(f"file_path_reverse updated to: {file_path_reverse}")

    #Remove the .fastq from the file paths
    file_name_forward = os.path.splitext(file_name_forward)[0]
    file_name_reverse = os.path.splitext(file_name_reverse)[0]
    print(f"file_name_forward updated to: {file_name_forward}")
    print(f"file_name_reverse updated to: {file_name_reverse}")

    # Construct the output file paths with '_trim_R1' and '_trim_R2' in the names
    output_file_forward = os.path.join(processed_folder_path, f"{file_name_forward}_trim.fastq")
    output_file_reverse = os.path.join(processed_folder_path, f"{file_name_reverse}_trim.fastq")
    print(f"output_file_forward is: {output_file_forward}")
    print(f"output_file_reverse is: {output_file_reverse}")

    #Do some string work to obtain output name
    common_name = f"{file_name_forward.split('_')[0]}_001"
    print(f"common name is: {common_name}")

    # Construct the output SAM file path
    output_sam = os.path.join(outputs_folder_path, f"{common_name}_aligned.sam")

    output_bam_unsorted = os.path.join(outputs_folder_path, f"{common_name}.unsorted.bam")
    output_bam_sorted = os.path.join(outputs_folder_path, f"{common_name}.bam")
    output_counts = os.path.join(outputs_folder_path, f"{common_name}_counts.txt")

    if not os.path.exists(output_counts):

        # Run trimmomatic PE for the pair
        os.system(f"trimmomatic PE -phred33 {file_path_forward} {file_path_reverse} {output_file_forward} /dev/null {output_file_reverse} /dev/null ILLUMINACLIP:TruSeq3-PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
        print(f"Checking outputs of trimmmomatic: output_file_forward: {output_file_forward}, output_file_reverse: {output_file_reverse}")
        os.system(f"fastqc {output_file_forward}")
        os.system(f"fastqc {output_file_reverse}")

        #Align with hisat2
        os.system(f"hisat2 -p 20 -x /nam-99/ablage/nam/ss_rnaseq/index_new -1 {output_file_forward} -2 {output_file_reverse} -S {output_sam}")

        #Convert sam to bam
        os.system("module load samtools")
        os.system(f"samtools view -b {output_sam} > {output_bam_unsorted}")
        os.system(f"samtools sort {output_bam_unsorted} > {output_bam_sorted}")

        #Extract counts
        os.system(f"featureCounts -t mRNA -g Parent -a /nam-99/ablage/nam/ss_rnaseq/Triticum_aestivum.IWGSC.58.gff3.gz -o {output_counts} -p {output_bam_sorted}")

        #Clean up
        os.system(f"rm -rf {output_bam_unsorted}")
        os.system(f"rm -rf {output_sam}")
        os.system(f"mv {processed_folder_path}/*.zip {cleanup_folder_path}")
        os.system(f"mv {processed_folder_path}/*.html {cleanup_folder_path}")

