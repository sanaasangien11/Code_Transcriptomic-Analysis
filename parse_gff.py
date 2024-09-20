from Bio import SeqIO

# Function to parse GFF file and extract gene IDs and CDS IDs
def parse_gff(gff_file):
    gene_to_cds_mapping = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    seq_id = parts[0]
                    feature_type = parts[2]
                    attributes = parts[8]
                    attributes_dict = dict(item.split('=') for item in attributes.split(';') if '=' in item)
                    if feature_type == 'gene' and 'ID' in attributes_dict:
                        gene_id = attributes_dict['ID']
                        if gene_id.startswith('ID=gene-'):
                            cds_id = attributes_dict.get('Parent', '')
                            if cds_id.startswith('ID=cds-'):
                                gene_to_cds_mapping[gene_id] = cds_id
    return gene_to_cds_mapping

# Parse GFF file to get gene IDs and their corresponding CDS IDs
gff_file = 'genomic.gff'
gene_to_cds_mapping = parse_gff(gff_file)

# Print or process the mapping as needed
for gene_id, cds_id in gene_to_cds_mapping.items():
    print(f"Gene ID: {gene_id}, CDS ID: {cds_id}")







