# --- 1. SETUP ---
# Load all necessary libraries for the entire workflow
library(minfi)
library(limma)
library(ggpubr)
library(DMRcate)
library(reshape2)
library(missMethyl)
library(doParallel)
library(TCGAbiolinks)
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(sesame)
library(sesameData)
library(DMRcatedata)
library(GenomicRanges)
library(rtracklayer) # Added for liftOver
library(readr)
library(dplyr)

# !!! SET YOUR PROJECT AND DIRECTORIES HERE !!!
project <- "TCGA-LIHC"
setwd("C:/Users/konst/OneDrive/Documents/MethylationAnalyis")

# --- 2. DATA DOWNLOAD AND PREPARATION ---
# Query, download, and prepare data using the project variable
query_met <- GDCquery(project = project,
                      data.category = "DNA Methylation",
                      data.type = "Methylation Beta Value",
                      platform = "Illumina Human Methylation 450")


GDCdownload(query_met, directory = "GDCdata_beta_values/", files.per.chunk = 20)
data.met <- GDCprepare(query_met, directory = "GDCdata_beta_values/")

# --- 3. DATA FILTERING ---
# Filter for paired samples and remove unwanted probes
paired.sample.ids <- data.met[, duplicated(data.met$patient)]$patient
data.met <- data.met[, !((data.met$patient %in% paired.sample.ids) & (data.met$definition == "Primary solid Tumor"))]
data.met <- data.met[, !data.met$definition == "Recurrent Solid Tumor"]
data.met <- data.met[!data.met@rowRanges$probeType == "rs", ]
data.met <- data.met[!data.met@rowRanges@seqnames %in% c("chrX", "chrY"),]
data.met <- data.met[rowSums(is.na(assay(data.met))) == 0, ]

# Clean up group names
data.met$definition <- gsub("Primary solid Tumor", "Primary_solid_Tumor", data.met$definition)
data.met$definition <- gsub("Solid Tissue Normal", "Solid_Tissue_Normal", data.met$definition)


# --- 4. VISUALIZATION PLOTS ---
# All plot filenames now include the project name automatically
df <- data.frame("Sample.mean" = colMeans(assay(data.met), na.rm = TRUE), "Groups" = data.met$definition)

ggboxplot(data = df, y = "Sample.mean", x = "Groups", color = "Groups", add = "jitter",
          ylab = "Mean DNA methylation (beta-values)", xlab = "",
          title = paste("Mean DNA Methylation Levels:", project)) +
  stat_compare_means()
ggsave(filename = sprintf("%s_mean_methylation.png", project), bg = "white", width = 12, height = 10, dpi = 300)

# Density plot
beta_values <- as.data.frame(t(assay(data.met)))
beta_values$Group <- as.factor(data.met$definition)
ggplot(melt(beta_values), aes(x = value, color = Group, fill = Group)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  labs(title = paste("Density Plot of Beta Values:", project), x = "Beta Value", y = "Density") +
  theme_minimal(base_size = 14)
ggsave(filename = sprintf("%s_density_plot.png", project), bg = "white", width = 12, height = 10, dpi = 300)


# MDS plot
beta_values2 <- t(assay(data.met))
dist_matrix <- dist(beta_values2)
mds_result <- cmdscale(dist_matrix, k = 2)
mds_df <- as.data.frame(mds_result)
colnames(mds_df) <- c("MDS1", "MDS2")
mds_df$Group <- as.factor(data.met$definition)
ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = paste("MDS Plot of Methylation Data:", project), x = "MDS Dimension 1", y = "MDS Dimension 2") +
  theme_minimal(base_size = 14)
ggsave(filename = sprintf("%s_mds_plot.png", project), bg = "white", width = 12, height = 10, dpi = 300)

# --- 5. DIFFERENTIAL METHYLATION ANALYSIS (DMC & DMR) ---
# DMC Analysis
dmc <- TCGAanalyze_DMC(data = data.met,
                       groupCol = "definition",
                       group1 = "Primary_solid_Tumor",
                       group2 = "Solid_Tissue_Normal",
                       p.cut = 0.05,
                       diffmean.cut = 0.15,
                       save = FALSE,
                       cores = 1)

dmc_df <- as.data.frame(dmc)
meth_probes_df <- as.data.frame(data.met@rowRanges)[, c("seqnames", "start", "end", "strand")]
dmc_df <- merge(dmc_df, meth_probes_df, by = 0)
names(dmc_df)[1] <- "probeID"

# DMR Analysis
pheno <- data.met$definition
design.matrix <- model.matrix(~0 + pheno)
colnames(design.matrix) <- gsub("pheno", "", colnames(design.matrix))
cont.matrix <- makeContrasts(Cancer.vs.Normal = Primary_solid_Tumor - Solid_Tissue_Normal, levels = design.matrix)

GRset <- GenomicRatioSet(Beta = assay(data.met), gr = rowRanges(data.met), colData = colData(data.met))
annotation(GRset) <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")

myAnnotation <- cpg.annotate(object = GRset, datatype = "array", what = "Beta",
                             analysis.type = "differential", design = design.matrix,
                             contrasts = TRUE, cont.matrix = cont.matrix,
                             coef = "Cancer.vs.Normal")

dmrcoutput <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 0.05)
dmr_df <- as.data.frame(extractRanges(dmrcoutput))

print(paste("Analysis for", project, "is complete. Starting annotation..."))


# --- 6. ANNOTATION ---
# Download the liftOver chain file (only needs to be done once)
chain_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
chain_file <- "hg38ToHg19.over.chain.gz"
if (!file.exists(chain_file)) {
  download.file(chain_url, chain_file)
}
ch <- import.chain(chain_file)

# List of your hg38 BED files and their desired names
bed_files_hg38 <- c(
  # TF_Promoter = "TF.promoters.bed",
  miRNA_Gene_Body = "miRNA.gene_body.bed",
  # miRNA_Promoter = "miRNA.promoters.bed",
  # TF_Gene_Body = "TF.gene_body.bed"
)

# Initialize a list to hold the converted hg19 GRanges objects
annotations <- list()

# Loop through each file, convert it, and add to the list
for (name in names(bed_files_hg38)) {
  file_path <- bed_files_hg38[name]
  print(paste("Processing:", file_path))
  
  gr_hg38 <- import(file_path, format = "BED")
  gr_hg19_list <- liftOver(gr_hg38, ch)
  annotations[[name]] <- unlist(gr_hg19_list)
}

print("LiftOver complete. All annotation files are now in hg19.")

## 6.2 Perform Overlap Analysis ##
# Create GenomicRanges objects from the analysis results (already in memory)
dmc_gr <- makeGRangesFromDataFrame(dmc_df, keep.extra.columns = TRUE)
dmr_gr <- makeGRangesFromDataFrame(dmr_df, keep.extra.columns = TRUE)

# The 'annotations' list now contains the hg19 GRanges objects ready for use

# Function to add annotation based on overlap
find_overlaps_and_annotate <- function(target_gr, annotation_gr, annotation_name) {
  overlaps <- findOverlaps(target_gr, annotation_gr)
  
  if (!paste0("overlaps_", annotation_name) %in% names(mcols(target_gr))) {
    mcols(target_gr)[[paste0("overlaps_", annotation_name)]] <- "No"
  }
  
  mcols(target_gr)[[paste0("overlaps_", annotation_name)]][queryHits(overlaps)] <- "Yes"
  
  return(target_gr)
}

# Loop through each annotation type for DMCs and DMRs
for (name in names(annotations)) {
  dmc_gr <- find_overlaps_and_annotate(dmc_gr, annotations[[name]], name)
  dmr_gr <- find_overlaps_and_annotate(dmr_gr, annotations[[name]], name)
}


# --- 7. SAVE FINAL ANNOTATED RESULTS ---
annotated_dmc_df <- as.data.frame(dmc_gr)
annotated_dmr_df <- as.data.frame(dmr_gr)

# Save the final annotated files
write.csv(annotated_dmc_df, sprintf("%s_dmc_annotated.csv", project), row.names = FALSE)
write.csv(annotated_dmr_df, sprintf("%s_dmr_annotated.csv", project), row.names = FALSE)
