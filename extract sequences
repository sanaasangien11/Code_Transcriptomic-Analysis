awk 'BEGIN {
    # Load gene names into an array
    while ((getline line < "genenames.txt") > 0) {
        genes[line] = 1
    }
    RS = ">"
    FS = "\n"
}

NR > 1 {
    header = $1
    seq = ""
    for (i = 2; i <= NF; i++) {
        seq = seq $i "\n"
    }
    # Extract the gene ID from the header
    split(header, arr, " ")
    gene_id = arr[1]
    # Check if gene ID is in the array
    if (gene_id in genes) {
        print ">" header "\n" seq
    }
}' Aet.pgsb.Mar2019.all.chr_scaffolds.aa.fa

