import pandas as pd
from Bio import SeqIO

textfile=pd.read_csv("gnn_wheat_edited.csv", header=None, names=['genes'])
textfile_list=textfile['genes'].tolist()

fasta_sequences = SeqIO.parse("Triticum_aestivum.IWGSC.pep.all.fa", "fasta")

sequence_ids=[]

for seq_record in fasta_sequences:
     sequence_id=seq_record.id
     for gene_part in textfile_list:
             if gene_part in sequence_id:
                 sequence_ids.append(sequence_id)    
                   
print(sequence_ids)

