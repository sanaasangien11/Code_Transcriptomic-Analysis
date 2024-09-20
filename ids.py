# Read gene names from t_genenames.txt
gene_names = set()
with open('t_genenames_trimmed.txt', 'r') as gene_file:
    for line in gene_file:
        gene_names.add(line.strip())

# Create a new file to write filtered data
with open('filtered_gene_cds_mapping.csv', 'w') as output_file:
    # Read gene_cds_mapping.txt, filter rows based on gene names, and write to output file
    with open('gene_cds_mapping.csv', 'r') as input_file:
        for line in input_file:
            gene_name = line.split()[0]  # Get the gene name from the first column
            if gene_name in gene_names:
                output_file.write(line)
