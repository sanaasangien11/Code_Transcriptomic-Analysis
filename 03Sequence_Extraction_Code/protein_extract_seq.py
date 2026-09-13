def extract_protein_sequences(protein_fasta, gene_names_file):
  """
  Extracts protein sequences for genes in gene_names_file from a protein FASTA file.

  Args:
      protein_fasta (str): Path to the protein FASTA file (assumed format: >protein_id|other_info).
      gene_names_file (str): Path to the file containing gene names (one per line).

  Returns:
      dict: A dictionary where keys are gene names and values are corresponding protein sequences.
  """

  # Read gene names
  with open(gene_names_file, 'r') as f:
    gene_names = [line.strip() for line in f]

  # Extract protein sequences
  final_protein_seqs = {}
  current_protein_name = None
  current_protein_seq = ""
  with open(protein_fasta, 'r') as f:
    for line in f:
      if line.startswith('>'):
        # New protein entry
        if current_protein_name:
          # Check previous protein for gene match and add sequence if applicable
          if current_protein_name in gene_names:
            final_protein_seqs[current_protein_name] = current_protein_seq.strip()
        current_protein_name = line.strip()[1:]
        current_protein_seq = ""
      else:
        # Append sequence line
        current_protein_seq += line.strip()
    # Handle the last protein entry
    if current_protein_name:
      if current_protein_name in gene_names:
        final_protein_seqs[current_protein_name] = current_protein_seq.strip()

  return final_protein_seqs

# Example usage
protein_fasta = "Triticum_aestivum.IWGSC.pep.all.fa"
gene_names_file = "gnn_wheat.txt"

protein_sequences = extract_protein_sequences(protein_fasta, gene_names_file)

# Access extracted protein sequences
for gene, protein_seq in protein_sequences.items():
  print(f"Gene: {gene}, Protein Sequence: {protein_seq}")

