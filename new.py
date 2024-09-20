import os
import subprocess

# Define the paths
index_file = "/nam-99/ablage/nam/ss_rnaseq/tauschii_data/hisat2_index/tauschi.8.ht2"
chrom_dir = "/nam-99/ablage/nam/ss_rnaseq/tauschii_data/CHROM"
index_output_prefix = "/nam-99/ablage/nam/ss_rnaseq/tauschii_data/hisat2_index/tauschi"

# Check if the HISAT2 index already exists
if not os.path.exists(index_file):
    # Verify the chromosome directory exists
    if not os.path.exists(chrom_dir):
        print(f"Error: Chromosome directory {chrom_dir} does not exist.")
    else:
        # Get the list of chromosome files with the correct extensions
        chrom_files = [os.path.join(chrom_dir, x) for x in os.listdir(chrom_dir) if x.endswith(('.fa', '.fasta'))]
        
        # Ensure there are chromosome files to process
        if not chrom_files:
            print(f"Error: No chromosome files found in {chrom_dir}.")
        else:
            # Join the chromosome files with commas
            chroms = ",".join(chrom_files)
            
            # Build the HISAT2 index
            print(f"Building HISAT2 index with the following chromosome files: {chroms}")
            try:
                result = subprocess.run(f"hisat2-build -p 7 {chroms} {index_output_prefix}", shell=True, check=True)
                print("HISAT2 index built successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while building HISAT2 index: {e}")

