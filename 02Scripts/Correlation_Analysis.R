setwd("/run/user/1001/gvfs/smb-share:server=filer-5.ipk-gatersleben.de,share=agruppen/CSF/Sanaa/")
#Load the required packages
library(tximport)
library(dplyr)
library(edgeR)
library(tibble)
library(WGCNA)
library(tidyverse)
library(DESeq2)

#Read in the tsv file
speltoides <- read.delim("/run/user/1001/gvfs/smb-share:server=filer-5.ipk-gatersleben.de,share=agruppen/CSF/Sanaa/coexp_aegilops_counts/combined_abundance.tsv", header = TRUE, sep = "\t")

#Rename the column 'target_id'
aeg_speltoides <- speltoides %>%
    rename(transcript_id = target_id)

#import gtf file
library(rtracklayer)
gtf_file <- rtracklayer::import("Final_transcriptome_stringTie.gtf")
aesAnno <- as_tibble(gtf_file)
aes_anno <- aesAnno[,c("transcript_id", "gene_id")] 

# Replace transcript IDs with gene IDs directly
aeg_speltoides$transcript_id <- aes_anno$gene_id[match(aeg_speltoides$transcript_id, aes_anno$transcript_id)]
# Perform a 'left join'to replace transcript IDs with gene IDs
#aeg_speltoides <- left_join(aeg_speltoides, aes_anno, by = "transcript_id")
#Remove the duplicate rows
#aeg_speltoides <- distinct(aeg_speltoides)
#Migrate 'gene_id' as the first column
#aeg_speltoides <- select(aeg_speltoides, gene_id, everything())
#Rename the column 'transcript_id'
aeg_speltoides <- rename(aeg_speltoides, gene_id = transcript_id)
#Use 'aggregate function to get the gene level counts
aeg_speltoides <- aggregate(. ~ gene_id, data = aeg_speltoides, FUN = sum)
#Remove the columns that end with length or tpm
aeg_speltoides <- select(aeg_speltoides, -ends_with("length"), -ends_with("tpm"))
#edit column names
names(aeg_speltoides) <- sub("_est_counts", "", names(aeg_speltoides))

#Read in the study design from previous analysis
targets <- read.csv("sample_design.csv")
path <- file.path(targets$sample, "abundance.tsv")
#Now import the aegilops counts from the expression data generated from our own study
all(file.exists(path)) 

# import Kallisto transcript counts into R using Tximport
txi <- tximport(path, 
                type = "kallisto", 
                tx2gene = aes_anno,
                txIn = T,
                txOut = F, #How does the result change if this =FALSE vs =TRUE?
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = F)

#Extract the counts from txi object and extract the sample names from the target object
myCounts <- txi$counts
samplelables <- targets$sample

#Make the data frame of the counts extracted 
aeg_counts_df <- as_tibble(myCounts, rownames = "gene_id") #data frame does not consider rownames so we force it to consider it 
# Add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(aeg_counts_df) <- c("gene_id", samplelables)
#Remove unwanted columns
aeg_counts_df <- select(aeg_counts_df, -starts_with("LCM"))
aeg_counts_df <- select(aeg_counts_df, -B0_Leaf_1, -B0_Leaf_2, -B0_Leaf_3, -B_Leaf_1, -B_Leaf_2, -B_Leaf_3)

#Merge the two data frames (one from public expression data of aegilops speltoids and 30 samples of aegilops speltoides expression data generated from our own study)
aeg_speltoides_final <- merge(aeg_speltoides, aeg_counts_df, by = "gene_id")

#Normalize your data
#cpm normalization
cpm_matrix <- as.matrix(aeg_speltoides_final[, -1])
cpm_matrix <- DGEList(cpm_matrix)
cpm_values <- cpm(cpm_matrix, log = TRUE)
cpm_values <- cbind(aeg_speltoides_final[, 1, drop = FALSE], cpm_values)

#Filter the data
#adjust the cutoff for the number of samples in the smallest group of comparison.
keepers <- rowSums(cpm_values>1)>=3
# now use base R's simple subsetting method to filter your DGEList based on the logical produced above
cpm_filtered <- cpm_values[keepers,]

#QC
#Detect the outliers genes
# Store gene_id column separately
gene_ids <- cpm_filtered$gene_id
# Remove gene_id column before applying goodSamplesGenes
data_without_gene_id <- cpm_filtered[, -1]
# Apply goodSamplesGenes function
outlier_genes <- goodSamplesGenes(t(data_without_gene_id))
# Reattach gene_id column to the result
outlier_genes$gene_id <- gene_ids
#remove genes detected as outliers
cpm_filtered <- cpm_filtered[outlier_genes$goodGenes == TRUE,]

#get sample names
samplenames <- colnames(cpm_filtered)
# Remove 'gene_id'
samplenames <- samplenames[-1]

#detect outlier samples- PCA
data_without_gene_id2 <- cpm_filtered[, -1]
pca.res <- prcomp(t(data_without_gene_id2), scale.=F, retx=T)
pca.data <- pca.res$x
pca.data <- as.data.frame(pca.data)

pca.var <- pca.res$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits=2)

ggplot(pca.data,aes(PC1,PC2, label=samplenames, color=samplenames)) +
  geom_point() +
  #geom_label() +
  #stat_ellipse(level = 0.95) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2], '%'))
write.csv(cpm_filtered, file = "gene_matrix.csv")

#NAM-99 server
#Co-expression analyisis
#aegilops speltoides

#read the gene matrix into server
gene_matrix <- read.csv("gene_matrix.csv", header = TRUE, row.names = 1)
#create the correlation matrix
correlation_matrix <- cor(t(gene_matrix))
# square the correlation matrix
squared_cor_matrix <- correlation_matrix^2
#rows to keep
rows <- c("MSTRG.36904","MSTRG.29143","MSTRG.56948")
new_matrix <- squared_cor_matrix[rows,]
# Find columns where all values are less than 0.7
cols <- which(all(new_matrix < 0.7))
subset_matrix <- new_matrix[, -cols]
#from the subset matrix, assign NA to values below than 0.70
subset_matrix_na <- subset_matrix
subset_matrix_na[subset_matrix_na < 0.70] <- NA
#save both the objects for further referance
write.csv(subset_matrix, file="subset_matrix.csv") #similar for 'NA' matrix
#get the column names from the subset matrix to create a gene vector
colnames <- colnames(subset_matrix)
#create gene vector
gene_vector <- as.vector(colnames)
#subset the gene vector from the squared correlation matrix previously created
subset_cor_matrix <- squared_cor_matrix[gene_vector, gene_vector] #create the subset matrix for correlation_matrix as well without square
#save the matrix in csv format to load into R for visualization
write.csv(subset_cor_matrix, file="subset_cor_matrix.csv")
#load the subset cor matrix in R and visualize using 'ggcorplot'.
subset_cor_matrix <- read.csv("subset_cor_matrix.csv", header = TRUE, row.names = 1)

#visualize using the multi dimensional scaling 
#perform it on both subset matrixes without square and with square
#create the dissimilarity matrix and perform classical MDS on that
mds.corr_speltoides <- (1 - subset_corr_matrix) %>%
     cmdscale() %>%
     as_tibble()
colnames(mds.corr_speltoides) <- c("Dim.1", "Dim.2")
#create the scatter plot of the MDS result 
ggscatter(mds.corr_speltoides, x = "Dim.1", y = "Dim.2", 
          size = 1,
          color = ifelse(colnames(subset_corr_matrix) %in% candidate_genes, "red", "black"),
          label =  ifelse(colnames(subset_corr_matrix) %in% candidate_genes, colnames(subset_corr_matrix), ""), 
          repel = TRUE) +
  labs(title = "MDS Plot of Aegilops Speltoides Genes")

#visualize the subsetted correlation matrixes using the heatmap as well
ggcorrplot(subset_corr_matrix, hc.order = TRUE, outline.col = "white") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +  # Remove axis text
  ggtitle("Correlation Heatmap of Gene Expression in Aegilops Speltoides (squared)")  # Add title

#continue the same code with the gene expression data of wheat and tauschii
#wheat
#read in the gene count data
wheat <- read.csv("/run/user/1001/gvfs/smb-share:server=filer-5.ipk-gatersleben.de,share=agruppen/CSF/Sanaa/wheat_counts.tsv", header = TRUE, row.names = 1)

#normalize the counts
cpm_matrix_wheat <- DGEList(wheat)
cpm_values_wheat <- cpm(cpm_matrix_wheat, log = TRUE)
cpm_values_wheat <- cbind(wheat[, 1, drop = FALSE], cpm_values_wheat)
 
#Filter the data
#adjust the cutoff for the number of samples in the smallest group of comparison.
keepers_wheat <- rowSums(cpm_values_wheat>1)>=3
cpm_filtered_wheat <- cpm_values_wheat[keepers_wheat,]

QC
#Detect the outliers genes
# Apply goodSamplesGenes function
outlier_genes_wheat <- goodSamplesGenes(t(cpm_filtered_wheat))
#remove genes detected as outliers
cpm_filtered_wheat <- cpm_filtered_wheat[outlier_genes_wheat$goodGenes == TRUE,]

#get sample names
samplenames_wheat <- colnames(cpm_filtered_wheat)

#detect outlier samples- PCA
pca.res <- prcomp(t(cpm_filtered_wheat), scale.=F, retx=T)
pca.data <- pca.res$x
pca.data <- as.data.frame(pca.data)

pca.var <- pca.res$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits=2)

ggplot(pca.data,aes(PC1,PC2, label=samplenames_wheat, color=samplenames_wheat)) +
  geom_point() +
  #geom_label() +
  #stat_ellipse(level = 0.95) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2], '%'))
#save the normalized gene count matrix
write.csv(cpm_filtered_wheat, file = "gene_matrix_wheat.csv")

#NAM-99 server
#Co-expression analyisis
#wheat

#read the gene matrix into server
gene_matrix_wheat <- read.csv("gene_matrix_wheat.csv", header = TRUE, row.names = 1)
#create the correlation matrix
correlation_matrix_wheat <- cor(t(gene_matrix_wheat))
# square the correlation matrix
squared_cor_matrix_wheat <- correlation_matrix_wheat^2
#rows to keep
rows_wheat <- c("gene:TraesCS6B02G126900","gene:TraesCS2B02G364900","gene:TraesCS4B02G196000")
new_matrix_wheat <- squared_cor_matrix_wheat[rows_wheat,]
# Find columns where all values are less than 0.7
cols_wheat <- which(all(new_matrix_wheat < 0.7))
subset_matrix_wheat <- new_matrix_wheat[, -cols_wheat]
#from the subset matrix, assign NA to values below than 0.70
subset_matrix_na <- subset_matrix
subset_matrix_na[subset_matrix_na < 0.70] <- NA
#save both the objects for further referance
write.csv(subset_matrix_wheat, file="subset_matrix_wheat.csv") #similar for 'NA' matrix
#get the column names from the subset matrix to create a gene vector
colnames_wheat <- colnames(subset_matrix_wheat)
#create gene vector
gene_vector_wheat <- as.vector(colnames_wheat)
#subset the gene vector from the squared correlation matrix previously created
subset_cor_matrix_wheat <- squared_cor_matrix_wheat[gene_vector_wheat, gene_vector_wheat] #create the subset matrix for correlation_matrix as well without square
#save the matrix in csv format to load into R for visualization
write.csv(subset_cor_matrix_wheat, file="subset_cor_matrix_wheat.csv")
#load the subset cor matrix in R and visualize using 'ggcorplot'.
subset_cor_matrix_wheat <- read.csv("subset_cor_matrix_wheat.csv", header = TRUE, row.names = 1)

#visualize using the multi dimensional scaling 
#perform it on both subset matrixes without square and with square
#create the dissimilarity matrix and perform classical MDS on that
mds.corr_wheat <- (1 - subset_corr_matrix_wheat) %>%
  cmdscale() %>%
  as_tibble()
colnames(mds.corr_wheat) <- c("Dim.1", "Dim.2")
#create the scatter plot of the MDS result 
ggscatter(mds.corr_wheat, x = "Dim.1", y = "Dim.2", 
          size = 1,
          color = ifelse(colnames(subset_corr_matrix_wheat) %in% candidate_genes, "red", "black"),
          label =  ifelse(colnames(subset_corr_matrix_wheat) %in% candidate_genes, colnames(subset_corr_matrix), ""), 
          repel = TRUE) +
  labs(title = "MDS Plot of Aegilops Speltoides Genes")

#visualize the subsetted correlation matrixes using the heatmap as well
ggcorrplot(subset_cor_matrix, hc.order = FALSE, outline.col = "white") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +  # Remove axis text
  ggtitle("Correlation Heatmap of Gene Expression in Aegilops Speltoides (squared)")  # Add title

#Aegilops tauschii 
#read in the gene count data
tauschii <- read.csv("/run/user/1001/gvfs/smb-share:server=filer-5.ipk-gatersleben.de,share=agruppen/CSF/Sanaa/tauschii_counts.csv", header = TRUE, row.names = 1)

#normalize the counts
cpm_matrix_tauschii <- DGEList(tauschii)
cpm_values_tauschii <- cpm(cpm_matrix_tauschii, log = TRUE)
cpm_values_tauschii <- cbind(wheat[, 1, drop = FALSE], cpm_values_wheat)

#Filter the data
#adjust the cutoff for the number of samples in the smallest group of comparison.
keepers_tauschii <- rowSums(cpm_values_tauschii>1)>=3
cpm_filtered_tauschii <- cpm_values_tauschii[keepers_tauschii,]

#QC
#Detect the outliers genes
# Apply goodSamplesGenes function
outlier_genes_tauschii <- goodSamplesGenes(t(cpm_filtered_tauschii))
#remove genes detected as outliers
cpm_filtered_tauschii <- cpm_filtered_tauschii[outlier_genes_tauschii$goodGenes == TRUE,]

#get sample names
samplenames_tauschii <- colnames(cpm_filtered_tauschii)

#detect outlier samples- PCA
pca.res <- prcomp(t(cpm_filtered_tauschii), scale.=F, retx=T)
pca.data <- pca.res$x
pca.data <- as.data.frame(pca.data)

pca.var <- pca.res$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits=2)

ggplot(pca.data,aes(PC1,PC2, label=samplenames_tauschii, color=samplenames_tauschii)) +
  geom_point() +
  #geom_label() +
  #stat_ellipse(level = 0.95) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2], '%'))
#save the normalized gene count matrix
write.csv(cpm_filtered_tauschii, file = "gene_matrix_tauschii.csv")

#NAM-99 server
#Co-expression analyisis
#wheat

#read the gene matrix into server
gene_matrix_tauschii <- read.csv("gene_matrix_tauschii.csv", header = TRUE, row.names = 1)
#create the correlation matrix
correlation_matrix_tauschii <- cor(t(gene_matrix_tauschii))
# square the correlation matrix
squared_cor_matrix_tauschii <- correlation_matrix_tauschii^2
#rows to keep
rows_tauschii <- c("LOC109760873","LOC109746448","LOC109786100")
new_matrix_tauschii <- squared_cor_matrix_tauschii[rows_tauschii,]
# Find columns where all values are less than 0.7
cols_tauschii <- which(all(new_matrix_tauschii < 0.7))
subset_matrix_tauschii <- new_matrix_tauschii[, -cols_tauschii]
#from the subset matrix, assign NA to values below than 0.70
subset_matrix_na_tauschii <- subset_matrix_tauschii
subset_matrix_na_tauschii[subset_matrix_na_tauschii < 0.70] <- NA
#save both the objects for further referance
write.csv(subset_matrix_tauschii, file="subset_matrix_tauschii.csv") #similar for 'NA' matrix










