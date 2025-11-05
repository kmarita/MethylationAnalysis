library(GenomicRanges)

main_dir <- "~/methylation_project"
results_dir <- file.path(main_dir, "methylation_analysis_results")

cancer_types <- c("TCGA-BRCA","TCGA-BLCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", 
                  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
                  "TCGA-LIHC","TCGA-LUSC", "TCGA-LUAD", "TCGA-PAAD", "TCGA-PCPG", 
                  "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", 
                  "TCGA-THYM", "TCGA-UCEC")

read_miRNA_dmcs <- function(cancer_type, results_dir) {
  file_path <- file.path(results_dir, cancer_type, 
                         paste0(cancer_type, "_DMCs_in_miRNA.tsv"))
  
  if (!file.exists(file_path)) {
    cat(sprintf("File not found for %s\n", cancer_type))
    return(NULL)
  }
  
  dmcs <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dmcs_sig <- dmcs[dmcs$Methylation_Status != "Not Significant", ]
  
  if (nrow(dmcs_sig) == 0) {
    cat(sprintf("No significant DMCs for %s\n", cancer_type))
    return(NULL)
  }
  
  cat(sprintf("Loaded %d significant miRNA DMCs from %s\n", nrow(dmcs_sig), cancer_type))
  return(dmcs_sig)
}

all_dmcs_list <- list()
for (cancer in cancer_types) {
  dmcs <- read_miRNA_dmcs(cancer, results_dir)
  if (!is.null(dmcs)) {
    all_dmcs_list[[cancer]] <- dmcs
  }
}

if (length(all_dmcs_list) == 0) {
  stop("No miRNA DMC files found with significant results.")
}

cat(sprintf("\nTotal cancer types with significant miRNA DMCs: %d\n\n", length(all_dmcs_list)))

all_cpgs <- data.frame()

cancer_names <- names(all_dmcs_list)
for (cancer in cancer_names) {
  dmcs <- all_dmcs_list[[cancer]]
  
  if (nrow(dmcs) > 0) {
    cpg_info <- dmcs[, c("ProbeID", "Chromosome", "Start", "End", "miRNA", "Methylation_Status")]
    cpg_info$Cancer <- cancer
    cpg_info$Coord <- paste(cpg_info$Chromosome, cpg_info$Start, cpg_info$End, sep = "_")
    all_cpgs <- rbind(all_cpgs, cpg_info)
  }
}

cpg_coords <- unique(all_cpgs$Coord)
cpg_concordance <- data.frame()

for (coord in cpg_coords) {
  coord_data <- all_cpgs[all_cpgs$Coord == coord, ]
  num_cancers <- nrow(coord_data)
  
  if (num_cancers >= 2) {
    all_hyper <- all(coord_data$Methylation_Status == "Hypermethylated")
    all_hypo <- all(coord_data$Methylation_Status == "Hypomethylated")
    
    if (all_hyper) {
      concordance <- "All_Hypermethylated"
    } else if (all_hypo) {
      concordance <- "All_Hypomethylated"
    } else {
      concordance <- "Discordant"
    }
    
    cpg_concordance <- rbind(cpg_concordance, data.frame(
      Coord = coord,
      ProbeID = coord_data$ProbeID[1],
      Chromosome = coord_data$Chromosome[1],
      Start = coord_data$Start[1],
      End = coord_data$End[1],
      miRNA = coord_data$miRNA[1],
      Num_Cancers = num_cancers,
      Cancers = paste(sort(unique(coord_data$Cancer)), collapse = ";"),
      Concordance = concordance,
      stringsAsFactors = FALSE
    ))
  }
}

concordant_cpgs <- cpg_concordance[cpg_concordance$Concordance != "Discordant", ]

cat(sprintf("\nTotal common CpGs (in 2+ cancers): %d\n", nrow(cpg_concordance)))
cat(sprintf("Concordant CpGs (all same direction): %d\n", nrow(concordant_cpgs)))
cat(sprintf("  All Hypermethylated: %d\n", sum(concordant_cpgs$Concordance == "All_Hypermethylated")))
cat(sprintf("  All Hypomethylated: %d\n", sum(concordant_cpgs$Concordance == "All_Hypomethylated")))
cat(sprintf("Discordant CpGs (mixed directions): %d\n\n", sum(cpg_concordance$Concordance == "Discordant")))

if (nrow(concordant_cpgs) > 0) {
  
  miRNA_concordant_counts <- data.frame()
  
  for (i in 1:nrow(concordant_cpgs)) {
    miRNAs <- unlist(strsplit(concordant_cpgs$miRNA[i], ";"))
    miRNAs <- miRNAs[miRNAs != "" & !is.na(miRNAs)]
    
    for (mirna in miRNAs) {
      concordance_type <- concordant_cpgs$Concordance[i]
      num_cancers <- concordant_cpgs$Num_Cancers[i]
      
      miRNA_concordant_counts <- rbind(miRNA_concordant_counts,
                                       data.frame(miRNA = mirna,
                                                  Concordance = concordance_type,
                                                  Num_Cancers = num_cancers,
                                                  stringsAsFactors = FALSE))
    }
  }
  
  miRNA_summary <- data.frame()
  for (mirna in unique(miRNA_concordant_counts$miRNA)) {
    mirna_data <- miRNA_concordant_counts[miRNA_concordant_counts$miRNA == mirna, ]
    
    total_concordant <- nrow(mirna_data)
    all_hyper <- sum(mirna_data$Concordance == "All_Hypermethylated")
    all_hypo <- sum(mirna_data$Concordance == "All_Hypomethylated")
    num_cancers <- length(unique(unlist(strsplit(
      paste(concordant_cpgs$Cancers[concordant_cpgs$miRNA == mirna | 
                                      grepl(mirna, concordant_cpgs$miRNA)], collapse = ";"), 
      ";"))))
    
    miRNA_summary <- rbind(miRNA_summary,
                           data.frame(miRNA = mirna,
                                      Total_Concordant_CpGs = total_concordant,
                                      All_Hypermethylated_CpGs = all_hyper,
                                      All_Hypomethylated_CpGs = all_hypo,
                                      Num_Cancers = num_cancers,
                                      stringsAsFactors = FALSE))
  }
  
  miRNA_summary <- miRNA_summary[order(-miRNA_summary$Total_Concordant_CpGs), ]
  
  output_file1 <- file.path(results_dir, "concordant_miRNA_DMCs_summary.tsv")
  write.table(miRNA_summary, file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  output_file2 <- file.path(results_dir, "concordant_CpGs_detailed.tsv")
  write.table(concordant_cpgs, file = output_file2, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Analysis complete.\n")
  cat(sprintf("Concordant miRNA summary: %s\n", output_file1))
  cat(sprintf("Concordant CpGs detailed: %s\n", output_file2))
  
  cat("\nTop 30 miRNAs with most concordant CpGs:\n")
  print(head(miRNA_summary, 30))
  
} else {
  cat("No concordant CpGs found.\n")
}