library(GenomicRanges)
library(dplyr)

main_dir <- "~/methylation_project"
results_dir <- file.path(main_dir, "methylation_analysis_results")

cancer_types <- c("TCGA-BRCA","TCGA-BLCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", 
                  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP",
                  "TCGA-LIHC","TCGA-LUSC", "TCGA-LUAD", "TCGA-PAAD", "TCGA-PCPG", 
                  "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", 
                  "TCGA-THYM", "TCGA-UCEC")

read_miRNA_dmrs <- function(cancer_type, results_dir) {
  file_path <- file.path(results_dir, cancer_type, 
                         paste0(cancer_type, "_DMRs_in_miRNA.tsv"))
  
  if (!file.exists(file_path)) {
    cat(sprintf("File not found for %s\n", cancer_type))
    return(NULL)
  }
  
  dmrs <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dmrs_sig <- dmrs[dmrs$DMR_Status !="Not Significant", ]
  
  if (nrow(dmrs_sig) == 0) {
    cat(sprintf("No significant DMRs for %s\n", cancer_type))
    return(NULL)
  }
  
  cat(sprintf("Loaded %d significant miRNA DMRs from %s\n", nrow(dmrs_sig), cancer_type))
  return(dmrs_sig)
}

all_dmrs_list <- list()
for (cancer in cancer_types) {
  dmrs_sig <- read_miRNA_dmrs(cancer, results_dir)
  if (!is.null(dmrs_sig)) {
    all_dmrs_list[[cancer]] <- dmrs_sig
  }
}

if (length(all_dmrs_list) == 0) {
  stop("No miRNA DMR files found with significant results.")
}

cat(sprintf("\nTotal cancer types with significant miRNA DMRs: %d\n\n", length(all_dmrs_list)))

all_dmrs <- data.frame()

cancer_names <- names(all_dmrs_list)
for (cancer in cancer_names) {
  dmrs_sig <- all_dmrs_list[[cancer]]
  
  if (nrow(dmrs_sig) > 0) {
    dmr_info <- dmrs_sig[, c("DMR_idx", "Chromosome", "Start", "End", "miRNA", "DMR_Status")]
    dmr_info$Cancer <- cancer
    dmr_info$Coord <- paste(dmr_info$Chromosome, dmr_info$Start, dmr_info$End, sep = "_")
    all_dmrs <- rbind(all_dmrs, dmr_info)
  }
}

dmr_coords <- unique(all_dmrs$Coord)
dmr_concordance <- data.frame()

for (coord in dmr_coords) {
  coord_data <- all_dmrs[all_dmrs$Coord == coord, ]
  num_cancers <- nrow(coord_data)
  
  if (num_cancers >= 10) {
    all_hyper <- all(coord_data$DMR_Status == "Hypermethylated")
    all_hypo <- all(coord_data$DMR_Status == "Hypomethylated")
    
    if (all_hyper) {
      concordance <- "All_Hypermethylated"
    } else if (all_hypo) {
      concordance <- "All_Hypomethylated"
    } else {
      concordance <- "Discordant"
    }
    
    dmr_concordance <- rbind(dmr_concordance, data.frame(
      Coord = coord,
      DMR_idx = coord_data$DMR_idx[1],
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

concordant_dmrs <- dmr_concordance[dmr_concordance$Concordance != "Discordant", ]

cat(sprintf("\nTotal common DMRs (in 2+ cancers): %d\n", nrow(dmr_concordance)))
cat(sprintf("Concordant DMRs (all same direction): %d\n", nrow(concordant_dmrs)))
cat(sprintf("  All Hypermethylated: %d\n", sum(concordant_dmrs$Concordance == "All_Hypermethylated")))
cat(sprintf("  All Hypomethylated: %d\n", sum(concordant_dmrs$Concordance == "All_Hypomethylated")))
cat(sprintf("Discordant DMRs (mixed directions): %d\n\n", sum(dmr_concordance$Concordance == "Discordant")))

# Create miRNA x Cancer matrix for DMRs
if (nrow(concordant_dmrs) > 0) {
  
  # Create a list to store miRNA-cancer-DMR relationships
  mirna_cancer_dmr <- list()
  
  for (i in 1:nrow(concordant_dmrs)) {
    # Parse miRNAs (handle multiple miRNAs per DMR)
    miRNAs <- unlist(strsplit(concordant_dmrs$miRNA[i], ";"))
    miRNAs <- miRNAs[miRNAs != "" & !is.na(miRNAs)]
    miRNAs <- trimws(miRNAs) #removes spaces
    
    # Parse cancers
    cancers <- unlist(strsplit(concordant_dmrs$Cancers[i], ";"))
    cancers <- trimws(cancers)
    
    # Get DMR coordinate
    dmr_coord <- concordant_dmrs$Coord[i]
    
    # Store each miRNA-cancer-DMR combination
    for (mirna in miRNAs) {
      for (cancer in cancers) {
        key <- paste(mirna, cancer, sep = "||")
        if (is.null(mirna_cancer_dmr[[key]])) {
          mirna_cancer_dmr[[key]] <- character()
        }
        mirna_cancer_dmr[[key]] <- c(mirna_cancer_dmr[[key]], dmr_coord)
      }
    }
  }
  
  # Get unique miRNAs and cancers
  all_mirnas <- unique(unlist(lapply(names(mirna_cancer_dmr), function(x) strsplit(x, "\\|\\|")[[1]][1])))
  all_mirnas <- sort(all_mirnas)
  
  all_cancers <- sort(cancer_names)
  
  # Create matrix
  dmr_matrix <- matrix(0, nrow = length(all_mirnas), ncol = length(all_cancers),
                       dimnames = list(all_mirnas, all_cancers))
  
  # Fill matrix with counts of unique DMRs
  for (mirna in all_mirnas) {
    for (cancer in all_cancers) {
      key <- paste(mirna, cancer, sep = "||")
      if (!is.null(mirna_cancer_dmr[[key]])) {
        # Count unique DMRs for this miRNA-cancer combination
        dmr_matrix[mirna, cancer] <- length(unique(mirna_cancer_dmr[[key]]))
      }
    }
  }
  
  # Convert matrix to data frame for easier viewing/export
  dmr_matrix_df <- as.data.frame(dmr_matrix)
  dmr_matrix_df <- cbind(miRNA = rownames(dmr_matrix_df), dmr_matrix_df)
  rownames(dmr_matrix_df) <- NULL
  
  # Save matrix to file
  output_matrix_file <- file.path(results_dir, "miRNA_cancer_common_dmr_matrix.tsv")
  write.table(dmr_matrix_df, file = output_matrix_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("\nmiRNA x Cancer DMR matrix created with %d miRNAs and %d cancer types\n", 
              length(all_mirnas), length(all_cancers)))
  cat(sprintf("Matrix saved to: %s\n\n", output_matrix_file))
  
  # Print summary statistics
  cat("Matrix summary:\n")
  cat(sprintf("  Total miRNAs: %d\n", nrow(dmr_matrix)))
  cat(sprintf("  Total cancer types: %d\n", ncol(dmr_matrix)))
  cat(sprintf("  Total non-zero entries: %d\n", sum(dmr_matrix > 0)))
  cat(sprintf("  Max DMRs for any miRNA-cancer pair: %d\n", max(dmr_matrix)))
  
  # Print top miRNAs (by total DMRs across all cancers)
  mirna_totals <- rowSums(dmr_matrix)
  top_mirnas <- head(sort(mirna_totals, decreasing = TRUE), 20)
  cat("\nTop 20 miRNAs by total common DMRs across all cancers:\n")
  print(top_mirnas)
  
  # Also save the original summary outputs
  miRNA_concordant_counts <- data.frame()
  
  for (i in 1:nrow(concordant_dmrs)) {
    miRNAs <- unlist(strsplit(concordant_dmrs$miRNA[i], ";"))
    miRNAs <- miRNAs[miRNAs != "" & !is.na(miRNAs)]
    
    for (mirna in miRNAs) {
      concordance_type <- concordant_dmrs$Concordance[i]
      num_cancers <- concordant_dmrs$Num_Cancers[i]
      
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
      paste(concordant_dmrs$Cancers[concordant_dmrs$miRNA == mirna | 
                                      grepl(mirna, concordant_dmrs$miRNA)], collapse = ";"), 
      ";"))))
    
    miRNA_summary <- rbind(miRNA_summary,
                           data.frame(miRNA = mirna,
                                      Total_Concordant_DMRs = total_concordant,
                                      All_Hypermethylated_DMRs = all_hyper,
                                      All_Hypomethylated_DMRs = all_hypo,
                                      Num_Cancers = num_cancers,
                                      stringsAsFactors = FALSE))
  }
  
  miRNA_summary <- miRNA_summary[order(-miRNA_summary$Total_Concordant_DMRs), ]
  
  output_file1 <- file.path(results_dir, "concordant_miRNA_DMRs_summary.tsv")
  write.table(miRNA_summary, file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  output_file2 <- file.path(results_dir, "concordant_DMRs_detailed.tsv")
  write.table(concordant_dmrs, file = output_file2, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\nAnalysis complete.\n")
  cat(sprintf("miRNA-Cancer DMR matrix: %s\n", output_matrix_file))
  cat(sprintf("Concordant miRNA DMRs summary: %s\n", output_file1))
  cat(sprintf("Concordant DMRs detailed: %s\n", output_file2))
  
} else {
  cat("No concordant DMRs found.\n")
}
