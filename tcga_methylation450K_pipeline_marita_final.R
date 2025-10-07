library("minfi")
library("limma")
library("ggpubr")
library("DMRcate")
library("reshape2")
library("missMethyl")
library("doParallel")
library("TCGAbiolinks")
library("RColorBrewer")
library("sesame")
library("sesameData")
library("DMRcatedata")
library("rtracklayer")
library("GenomicRanges")
library("dplyr")
library("tibble")
set.seed(1234)


dmc_cutoff <- 0.2
meandiff_cutoff <- 0.1

# Parallel backend is registered
cores_to_use <- 2
registerDoParallel(cores = cores_to_use)

cancer_type <- "TCGA-LIHC"
mirna_genes_bed <- "/Users/ross/Desktop/Projects/PasteurProjects/BSc theses/Methylation - Marita/miRNA_genes.sorted.bed"
main_dir <- "/Users/ross/Desktop/Projects/PasteurProjects/BSc theses/Methylation - Marita"

setwd(main_dir)
out_dir <- file.path(main_dir, "methylation_analysis_results", cancer_type)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)




########################### PART 1
# 1. Query & download
query_met <- GDCquery(project = cancer_type,
                      data.category = "DNA Methylation",
                      data.type = "Methylation Beta Value",
                      platform = "Illumina Human Methylation 450")

GDCdownload(query_met, directory = "GDCdata_beta_values", files.per.chunk = 20)
data.met <- GDCprepare(query_met, directory = "GDCdata_beta_values/")



# 2. Harmonize definitions
data.met$definition <- gsub("Primary solid Tumor", "Primary_solid_Tumor", data.met$definition)
data.met$definition <- gsub("Solid Tissue Normal", "Solid_Tissue_Normal", data.met$definition)

# 3. Keep only tumor + normal samples
keep_defs <- c("Primary_solid_Tumor", "Solid_Tissue_Normal")
data.sub  <- data.met[, data.met$definition %in% keep_defs]
coldat    <- as.data.frame(colData(data.sub))

# 4. Find patients with at least one of each
patients_all <- unique(coldat$patient)
tbl <- table(coldat$patient, coldat$definition)
patients_with_both <- rownames(tbl)[tbl[, "Primary_solid_Tumor"] > 0 & tbl[, "Solid_Tissue_Normal"] > 0]

# 5. Report & subset
if (length(patients_with_both) == 0) {
  cat(sprintf("No matched Primary_solid_Tumor & Solid_Tissue_Normal samples found in %s.\n", cancer_type))
  
} else {
  # Subset to only the paired patients
  data.paired <- data.sub[, colData(data.sub)$patient %in% patients_with_both]
  clinical.paired <- as.data.frame(colData(data.paired))
  
  # One row per patient
  clin.unique <- clinical.paired[!duplicated(clinical.paired$patient), ]
  
  cat(sprintf("Found %d paired patients in %s.\n", length(patients_with_both), cancer_type))
  cat("Number of paired patients: ", length(patients_with_both), "\n")
  avg_age_days  <- mean(as.numeric(clin.unique$age_at_diagnosis), na.rm = TRUE)
  avg_age_years <- avg_age_days / 365.25
  cat(sprintf("Average age of paired patients: %.1f days (%.1f years)\n", avg_age_days, avg_age_years))
  cat("Sex distribution of paired patients:\n")
  print(table(clin.unique$gender))
  rm(avg_age_days, avg_age_years)
  
  # overwrite for downstream
  data.met <- data.paired
  clinical.data <- clin.unique
  
}
rm(query_met, keep_defs, data.sub, coldat, patients_all, tbl, data.paired, clinical.paired, clin.unique)



########################### PART 2

# 1. Contingency test: is the distribution of Tumor vs. Normal samples different in M vs. F?
# Question: Are tumor vs. normal sample counts different in males vs. females?
tbl_gender_def <- table(Gender = colData(data.met)$gender, Definition = colData(data.met)$definition)

chisq_res <- chisq.test(tbl_gender_def)
print(chisq_res)

# 2. Global‐mean methylation by sex: t‐test & boxplot
# Question: Is the global mean methylation level itself different between sexes?
df_sex <- data.frame(SampleMean = colMeans(assay(data.met), na.rm = TRUE), Gender = factor(colData(data.met)$gender))

# boxplot + jitter + t-test p‐value
p <- ggboxplot(df_sex, x = "Gender", y = "SampleMean", color = "Gender", add = "jitter",
               ylab = "Mean DNA methylation (β-value)", xlab = "Sex", title = "Global mean methylation by sex") +
  stat_compare_means(method = "t.test")

print(p)
ggsave(file.path(out_dir, "1.mean_methylation_by_sex.png"), plot = p, bg = "white", width = 12, height = 10, dpi = 300)
rm(tbl_gender_def, chisq_res, df_sex, p)



########################### PART 3

# 1. Drop probes flagged by sesame
gr <- rowRanges(data.met)
keep_mask   <- !mcols(gr)$MASK_general  # Remove any probe failing mapping/SNP/repeat checks
keep_type   <- mcols(gr)$probeType == "cg"  # 2. Restrict to canonical CpG probes
keep_chr    <- !seqnames(gr) %in% c("chrX","chrY")  # 3. Drop probes on sex chromosomes
keep_complete <- rowSums(is.na(assay(data.met))) == 0  # 4. Drop any probes with missing values

# Combine
keep_probes <- keep_mask & keep_type & keep_chr & keep_complete
data.filt   <- data.met[keep_probes, ]

cat(sprintf("Filtered from %d → %d probes\n", nrow(data.met), nrow(data.filt)))

# Overwrite for plotting
data.met <- data.filt
rm(data.filt, gr, keep_mask, keep_type, keep_chr, keep_complete, keep_probes)



########################### PART 4

# Density plot
beta_df <- as.data.frame(t(assay(data.met)))
beta_df$Group <- factor(data.met$definition)
beta_melt <- melt(beta_df, id.vars = "Group", variable.name = "Probe", value.name = "BetaValue")

p1 <- ggplot(beta_melt, aes(x = BetaValue, color = Group, fill = Group)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  labs(title = "Density Plot of Beta Values", x = "Beta Value", y = "Density") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "top") 
ggsave(file.path(out_dir, "2.density_plot.png"), plot = p1, bg = "white", width = 12, height = 10, dpi = 300)


# MDS
beta_mat <- t(assay(data.met))
dist_mat <- dist(beta_mat)
mds_res <- cmdscale(dist_mat, k = 2)
mds_df  <- data.frame(MDS1 = mds_res[,1], MDS2 = mds_res[,2], Group = factor(data.met$definition))


p2 <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "MDS Plot of Methylation Data", x = "MDS Dimension 1", y = "MDS Dimension 2") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top") 
ggsave(file.path(out_dir, "3.mds_plot.png"), plot = p2, bg = "white", width = 12, height = 10, dpi = 300)
rm(beta_df, beta_melt, beta_mat, dist_mat, mds_res, mds_df, p1, p2)



########################### PART 5a

# 1. Run the DMC analysis
dmc <- TCGAanalyze_DMC(data          = data.met,
                       groupCol      = "definition",          # column in colData
                       group1        = "Primary_solid_Tumor", # tumor group
                       group2        = "Solid_Tissue_Normal", # normal group
                       p.cut         = 1,                     # FDR cutoff (no cutoff for now)
                       diffmean.cut  = 0,                     # allow all Δβ through
                       paired        = TRUE,                  # paired test
                       save          = FALSE,                 # no intermediate files
                       plot.filename = NULL,                  # no output
                       cores         = cores_to_use)


# 2. Results to data.frame and classify
dmc_df <- as.data.frame(dmc)
dmc_df$Status <- with(dmc_df, ifelse(p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal < 0.05 & mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal > dmc_cutoff, "Hypermethylated",
                                     ifelse(p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal < 0.05 & mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal < -dmc_cutoff, "Hypomethylated",
                                            "Not Significant")))

# Volcano plot
p_volc <- ggplot(dmc_df, aes(x = mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal, y = -log10(p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal), color = Status)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("Hypermethylated" = "#8ABFA3", "Hypomethylated" = "#F95454", "Not Significant" = "#B7B7B7")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-dmc_cutoff, dmc_cutoff),    linetype = "dashed") +
  labs(title = paste0("Volcano Plot of Differential Methylation: ", cancer_type),
       subtitle = paste0("With a Δβ cutoff of ", dmc_cutoff),
       x     = "Mean Δβ (Tumor vs. Normal)",
       y     = "-log10(Adjusted p-value)") +
  theme_minimal(base_size = 14) +
  theme(legend.position   = "top",
        axis.title        = element_text(size = 13),
        legend.title      = element_blank(),
        legend.text       = element_text(size = 11))
ggsave(file.path(out_dir, paste0(cancer_type, "_DMC_volcano.png")), plot = p_volc, bg = "white", width = 12, height = 10, dpi = 300)

# Quick general counts
cat("Hypermethylated CpGs: ", sum(dmc_df$Status == "Hypermethylated"), "\n")
cat("Hypomethylated  CpGs: ", sum(dmc_df$Status == "Hypomethylated"),  "\n")
cat("Total tested CpGs:    ", nrow(dmc_df), "\n")


####################### PART 5b

# 1. Grab the full GRanges of probes from your SummarizedExperiment
probe_gr <- rowRanges(data.met)

# 2. Subset to only those probes you tested (i.e. the ones in dmc_df)
probe_gr <- probe_gr[names(probe_gr) %in% rownames(dmc_df)]

# 3. Build probe annotation DF with Ensembl/RefSeq ID & HGNC symbol
probe_anno <- data.frame(ProbeID    = names(probe_gr),
                         Gene_ID    = probe_gr$gene,       # e.g. Ensembl/RefSeq IDs (semicolon-sep if multiple)
                         Gene_Name  = probe_gr$gene_HGNC,   # official HGNC symbols
                         Chromosome = as.character(seqnames(probe_gr)),
                         Start      = start(probe_gr),
                         End        = end(probe_gr),
                         Strand     = as.character(strand(probe_gr)),
                         stringsAsFactors = FALSE)

# 4. Merge coords + gene info into your DMC stats
dmc_annotated <- dmc_df %>% rownames_to_column("ProbeID") %>% left_join(probe_anno, by = "ProbeID")

# 5. Reorder & rename
dmc_tidy <- dmc_annotated %>%
  select(ProbeID,
         Gene_Name,
         Chromosome, Start, End, Strand,
         mean.Primary.solid.Tumor,
         mean.Solid.Tissue.Normal,
         mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal,
         p.value.Primary.solid.Tumor.Solid.Tissue.Normal,
         p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal,
         Status) %>%
  rename(Mean_Tumor          = mean.Primary.solid.Tumor,
         Mean_Normal         = mean.Solid.Tissue.Normal,
         Delta_Beta          = mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal,
         P_Value             = p.value.Primary.solid.Tumor.Solid.Tissue.Normal,
         Adj_P_Value         = p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal,
         Methylation_Status  = Status)


# 6. Output table
out_tsv <- file.path(out_dir, paste0(cancer_type, "_DMC_allProbes.tsv"))
write.table(dmc_tidy, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)


####################### PART 5c

# 1. Import your miRNA BED as a GRanges
mirna_gr <- import(mirna_genes_bed, format = "BED")
# miRNA ID is stored in a column named "miRNA_ID"
mcols(mirna_gr)$miRNA_ID <- mcols(mirna_gr)$name

# 2. Build GRanges of only those probes in dmc_df
probe_gr <- rowRanges(data.met)
dmc_probes <- rownames(dmc_df)
dmc_probe_gr <- probe_gr[dmc_probes]
# attach your statistics as metadata
mcols(dmc_probe_gr) <- cbind(mcols(dmc_probe_gr), dmc_df)

# 3. Find overlaps with miRNA genes (ignore strand)
hits <- findOverlaps(dmc_probe_gr, mirna_gr, ignore.strand = TRUE)
ov <- as.data.frame(hits)
ov$probeID <- names(dmc_probe_gr)[ov$queryHits]
ov$miRNA_ID <- mcols(mirna_gr)$miRNA_ID[ov$subjectHits]

# 4. Collapse multiple miRNA IDs per probe
miRNA_list <- ov %>% group_by(probeID) %>% summarise(miRNA = paste(unique(miRNA_ID), collapse = ";"))

# 5. Pull out only probes that hit any miRNA, and assemble a tidy data.frame
dmc_miRNA_tidy <- miRNA_list %>%
  left_join(as.data.frame(mcols(dmc_probe_gr)) %>% rownames_to_column("probeID"),by = "probeID") %>%
  # add back genomic coords
  left_join(data.frame(probeID    = names(dmc_probe_gr),
                       Chromosome = as.character(seqnames(dmc_probe_gr)),
                       Start      = start(dmc_probe_gr),
                       End        = end(dmc_probe_gr),
                       Strand     = as.character(strand(dmc_probe_gr)),
                       row.names  = NULL),
            by         = "probeID") %>%
  # reorder & rename columns
  select(ProbeID   = probeID,
         miRNA,
         Chromosome, Start, End, Strand,
         Mean_Tumor           = mean.Primary.solid.Tumor,
         Mean_Normal          = mean.Solid.Tissue.Normal,
         Delta_Beta           = mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal,
         P_Value              = p.value.Primary.solid.Tumor.Solid.Tissue.Normal,
         Adj_P_Value          = p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal,
         Methylation_Status   = Status)

# 6. Write to TSV
out_tsv <- file.path(out_dir, paste0(cancer_type, "_DMCs_in_miRNA.tsv"))
write.table(dmc_miRNA_tidy, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
rm(hits, ov, miRNA_list, dmc_miRNA_tidy, out_tsv, probe_gr, dmc_probe_gr, dmc_probes, dmc_df, dmc_cutoff, 
   dmc, dmc_annotated, dmc_tidy, p_volc, probe_anno)



########################### PART 6a

# 1. Create phenotype factor and design/contrast matrices
pheno <- factor(data.met$definition, levels = c("Solid_Tissue_Normal", "Primary_solid_Tumor"))
design.matrix <- model.matrix(~0 + pheno)
colnames(design.matrix) <- sub("pheno", "", colnames(design.matrix))
cont.matrix <- makeContrasts(Cancer.vs.Normal = Primary_solid_Tumor - Solid_Tissue_Normal, levels = design.matrix)

# 2. Build a GenomicRatioSet for DMRcate
GRset <- GenomicRatioSet(Beta    = assay(data.met),
                         gr      = rowRanges(data.met),
                         colData = colData(data.met))

# specify annotation / use the hg19 manifest for probe annotatio
annotation(GRset) <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")

# 3. Annotate CpGs and run DMRcate
myAnnotation <- cpg.annotate(object        = GRset,
                             datatype      = "array",
                             what          = "Beta",
                             analysis.type = "differential",
                             design        = design.matrix,
                             contrasts     = TRUE,
                             cont.matrix   = cont.matrix,
                             coef          = "Cancer.vs.Normal",
                             arraytype     = "450K",
                             annotation    =c(array="IlluminaHumanMethylation450k",
                                              annotation="ilmn12.hg19"))

dmrcoutput <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 1)

# 4. Extract significant regions
significantRegions <- extractRanges(dmrcoutput, genome = "hg38")
significantRegions_df <- as.data.frame(significantRegions)

# 5. Classify DMRs for volcano
significantRegions_df$Status <- with(significantRegions_df, 
                                     ifelse(min_smoothed_fdr < 0.05 & meandiff >  meandiff_cutoff, "Hypermethylated", 
                                            ifelse(min_smoothed_fdr < 0.05 & meandiff <  -meandiff_cutoff, "Hypomethylated",
                                                   "Not Significant")))

# 6. Volcano plot: median Δβ vs. -log10(FDR)


# 7. Print summary counts
cat("Hypermethylated DMRs:", sum(significantRegions_df$Status == "Hypermethylated"), "\n")
cat("Hypomethylated  DMRs:", sum(significantRegions_df$Status == "Hypomethylated"),  "\n")
cat("Total DMRs tested:   ", nrow(significantRegions_df), "\n")

# 9. Save full DMR table
dmr_csv <- file.path(out_dir, paste0(cancer_type, "_DMR_allRegions.csv"))
write.csv(significantRegions_df, dmr_csv, row.names = FALSE)


####################### PART 6b

# 1. Build a GRanges for DMRs
dmr_gr2 <- GRanges(seqnames = significantRegions_df$seqnames,
                   ranges   = IRanges(start = significantRegions_df$start,
                                      end   = significantRegions_df$end),
                   strand   = significantRegions_df$strand)

# 2. Overlap DMRs with miRNA genes (ignore strand)
hits <- findOverlaps(dmr_gr2, mirna_gr, ignore.strand = TRUE)
ov   <- as.data.frame(hits) %>% transmute(DMR_idx   = as.character(queryHits), miRNA_ID  = mcols(mirna_gr)$miRNA_ID[subjectHits])

# 3. Collapse multiple miRNA IDs per DMR
miRNA_list <- ov %>% group_by(DMR_idx) %>% summarize(miRNA = paste(unique(miRNA_ID), collapse = ";"))

# 4. Subset & annotate your DMR data.frame
dmr_miRNA <- significantRegions_df %>%
  rownames_to_column("DMR_idx") %>%
  inner_join(miRNA_list, by = "DMR_idx") %>%
  select(DMR_idx,
         miRNA,
         seqnames, start, end, strand,
         width,
         no.cpgs,
         meandiff,
         Stouffer,
         Fisher,
         min_smoothed_fdr,
         HMFDR,
         Status) %>%
  rename(Chromosome        = seqnames,
         Start             = start,
         End               = end,
         Strand            = strand,
         Region_Width      = width,
         Num_Probes        = no.cpgs,
         Median_DeltaBeta  = meandiff,
         Stouffer_stat     = Stouffer,
         Fisher_stat       = Fisher,
         Region_FDR        = min_smoothed_fdr,
         Region_HybridFDR  = HMFDR,
         DMR_Status        = Status)

# 5. Write out the miRNA‐overlapping DMRs
out_tsv <- file.path(out_dir, paste0(cancer_type, "_DMRs_in_miRNA.tsv"))
write.table(dmr_miRNA, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
