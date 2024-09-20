def parse_gff(gff_file):
    gene_to_cds_mapping = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    seq_id = parts[0]
                    feature_type = parts[2]
                    if feature_type == 'CDS':
                        attributes = parts[8]
                        attributes_dict = dict(item.split('=') for item in attributes.strip().split(';') if '=' in item)
                        cds_id = attributes_dict.get('ID', '')
                        gene_id = attributes_dict.get('gene', '')
                        if cds_id.startswith('cds-') and gene_id.startswith('LOC'):
                            gene_to_cds_mapping[gene_id] = cds_id
    return gene_to_cds_mapping

# Example usage
gff_file = 'genomic.gff'
gene_to_cds_mapping = parse_gff(gff_file)

# Print the mapping
for gene_id, cds_id in gene_to_cds_mapping.items():
    print(f"Gene ID: {gene_id}, CDS ID: {cds_id}")

