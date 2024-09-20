# Define path to your query and subject files 
query_file = 'speltoides_genes.txt'
subject_file = 'tauschii_genes.txt'

# Read in the query and subject files
with open(query_file, 'r') as query:
    query_data = query.readlines()

with open(subject_file, 'r') as subject_file:
    subject_data = subject.readlines()


# Iterate over each gene in the query file and each gene in the subject file (one file as rows and other file as columns)
for query_line in query_data:
    for subject_line in subject_data:

#write the blast command 




