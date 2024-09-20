def parse_gff_and_write_to_file(gff_file, output_file):
    with open(output_file, 'w') as out_file:
        with open(gff_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        feature_type = parts[2]
                        if feature_type == 'CDS':
                            attributes = parts[8]
                            attributes_dict = dict(item.split('=') for item in attributes.strip().split(';') if '=' in item)
                            cds_id = attributes_dict.get('ID', '')
                            gene_id = attributes_dict.get('gene', '')
                            if cds_id.startswith('cds-') and gene_id.startswith('LOC'):
                                out_file.write(f"Gene ID: {gene_id}, CDS ID: {cds_id}\n")

# Example usage
gff_file = 'genomic.gff'
output_file = 'gene_cds_mapping.txt'
parse_gff_and_write_to_file(gff_file, output_file)

