import pandas as pd
from Bio import SeqIO

textfile = pd.read_csv("gene_names_speltoides.csv", header=None, names=['genes'])
textfile_list = textfile['genes'].tolist()

fasta_sequences = SeqIO.parse("mercator4_awk_final_merged_output.fasta.transdecoder.pep", "fasta")

sequence_ids = []

for seq_record in fasta_sequences:
    sequence_id = seq_record.id
    for gene_part in textfile_list:
        if gene_part in sequence_id:
            sequence_ids.append(sequence_id)

output_file = "sequence_ids.txt"

# Write sequence IDs to the output file
with open(output_file, 'w') as file: [file.write(sequence_id + '\n') for sequence_id in sequence_ids]

print("Results saved to:", output_file)

