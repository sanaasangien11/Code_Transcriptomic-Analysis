from Bio import SeqIO
import os

#if not os.path.exists('Triticum_aestivum.IWGSC.dna.toplevel.fa.gz'):
    #os.system("wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz")


#if not os.path.exists('Triticum_aestivum.IWGSC.58.gff3.gz'):
    #os.system("wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/gff3/triticum_aestivum/Triticum_aestivum.IWGSC.58.gff3.gz")


if not os.path.exists('tauschii_data'):
    os.mkdir('tauschii_data')


if not os.path.exists('tauschii_data/CHROM'):
    os.mkdir('tauschii_data/CHROM')
    

if not os.path.exists('tauschii_data/hisat2_index'):
    os.mkdir('tauschii_data/hisat2_index')


output_path = "tauschii_data/CHROM"

# Decompress the file using gzip command line tool
#os.system("gzip -d Triticum_aestivum.IWGSC.dna.toplevel.fa.gz")

with open('Aet.pgsb.Mar2019.all.chr_scaffolds.cds.fa', 'rt') as handle:
    for rec in SeqIO.parse(handle, format='fasta'):
        if not os.path.exists(f"{output_path}/{rec.id}.fa"):
            SeqIO.write(SeqIO.SeqRecord(seq=rec.seq, id=rec.id), handle=f"{output_path}/{rec.id}.fa", format="fasta")
