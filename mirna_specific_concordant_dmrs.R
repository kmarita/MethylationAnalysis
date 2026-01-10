library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

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
  
  if (!  file.exists(file_path)) {
    cat(sprintf("File not found for %s\n", cancer_type))
    return(NULL)
  }
  
  dmrs <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dmrs_sig <- dmrs[dmrs$DMR_Status != "Not Significant", ]
  
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
  if (!  is.null(dmrs_sig)) {
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
  num_cancers <- length(unique(coord_data$Cancer))
  
  if (num_cancers >= 2) {
    all_hyper <- all(coord_data$DMR_Status == "Hypermethylated")
    all_hypo <- all(coord_data$DMR_Status == "Hypomethylated")
    
    if (all_hyper) {
      concordance <- "All_Hypermethylated"
    } else if (all_hypo) {
      concordance <- "All_Hypomethylated"
    } else {
      concordance <- "Discordant"
    }
    
    unique_mirnas <- unique(coord_data$miRNA)
    unique_mirnas <- unique_mirnas[unique_mirnas != "" & ! is.na(unique_mirnas)]
    
    dmr_concordance <- rbind(dmr_concordance, data.frame(
      Coord = coord,
      DMR_idx = paste(unique(coord_data$DMR_idx), collapse = ";"),
      Chromosome = coord_data$Chromosome[1],
      Start = coord_data$Start[1],
      End = coord_data$End[1],
      miRNA = paste(sort(unique_mirnas), collapse = ";"),
      Num_Cancers = num_cancers,
      Cancers = paste(sort(unique(coord_data$Cancer)), collapse = ";"),
      Concordance = concordance,
      stringsAsFactors = FALSE
    ))
  }
}

concordant_dmrs <- dmr_concordance[dmr_concordance$Concordance != "Discordant", ]

cat(sprintf("\nTotal DMR regions in 2+ cancers (exact coordinates): %d\n", nrow(dmr_concordance)))
cat(sprintf("Concordant DMRs (all same direction): %d\n", nrow(concordant_dmrs)))
cat(sprintf("  All Hypermethylated:   %d\n", sum(concordant_dmrs$Concordance == "All_Hypermethylated")))
cat(sprintf("  All Hypomethylated:  %d\n", sum(concordant_dmrs$Concordance == "All_Hypomethylated")))
cat(sprintf("Discordant DMRs (mixed directions): %d\n\n", sum(dmr_concordance$Concordance == "Discordant")))

cat("Distribution of DMRs by number of cancers:\n")
cancer_count_table <- table(dmr_concordance$Num_Cancers)
print(cancer_count_table)
cat("\n")

if (nrow(concordant_dmrs) > 0) {
  
  mirna_cancer_dmr <- list()
  
  for (i in 1:nrow(concordant_dmrs)) {
    miRNAs <- unlist(strsplit(concordant_dmrs$miRNA[i], ";"))
    miRNAs <- miRNAs[miRNAs != "" & !is.na(miRNAs)]
    miRNAs <- trimws(miRNAs)
    
    cancers <- unlist(strsplit(concordant_dmrs$Cancers[i], ";"))
    cancers <- trimws(cancers)
    
    dmr_coord <- concordant_dmrs$Coord[i]
    
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
  
  all_mirnas <- unique(unlist(strsplit(concordant_dmrs$miRNA, ";")))
  all_mirnas <- all_mirnas[all_mirnas != "" & !is.na(all_mirnas)]
  all_mirnas <- trimws(all_mirnas)
  all_mirnas <- sort(all_mirnas)
  
  all_cancers_with_dmrs <- unique(unlist(strsplit(concordant_dmrs$Cancers, ";")))
  all_cancers_with_dmrs <- trimws(all_cancers_with_dmrs)
  all_cancers_with_dmrs <- sort(all_cancers_with_dmrs)
  
  dmr_matrix <- matrix(0, nrow = length(all_mirnas), ncol = length(all_cancers_with_dmrs),
                       dimnames = list(all_mirnas, all_cancers_with_dmrs))
  
  for (mirna in all_mirnas) {
    for (cancer in all_cancers_with_dmrs) {
      key <- paste(mirna, cancer, sep = "||")
      if (!is.null(mirna_cancer_dmr[[key]])) {
        dmr_matrix[mirna, cancer] <- length(unique(mirna_cancer_dmr[[key]]))
      }
    }
  }
  
  dmr_matrix_df <- as.data.frame(dmr_matrix)
  dmr_matrix_df <- cbind(miRNA = rownames(dmr_matrix_df), dmr_matrix_df)
  rownames(dmr_matrix_df) <- NULL
  
  output_matrix_file <- file.path(results_dir, "miRNA_cancer_exact_dmr_matrixFINAL.tsv")
  write.table(dmr_matrix_df, file = output_matrix_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("\nmiRNA x Cancer DMR matrix created with %d miRNAs and %d cancer types\n", 
              length(all_mirnas), length(all_cancers_with_dmrs)))
  cat(sprintf("Matrix saved to: %s\n\n", output_matrix_file))
  
  cat("Matrix summary:\n")
  cat(sprintf("  Total miRNAs: %d\n", nrow(dmr_matrix)))
  cat(sprintf("  Total cancer types: %d\n", ncol(dmr_matrix)))
  cat(sprintf("  Total non-zero entries: %d\n", sum(dmr_matrix > 0)))
  cat(sprintf("  Max DMRs for any miRNA-cancer pair: %d\n", max(dmr_matrix)))
  
  mirna_totals <- rowSums(dmr_matrix)
  top_mirnas <- head(sort(mirna_totals, decreasing = TRUE), 20)
  cat("\nTop 20 miRNAs by total exact coordinate DMRs across all cancers:\n")
  print(top_mirnas)
  
  miRNA_concordant_counts <- data.frame()
  
  for (i in 1:nrow(concordant_dmrs)) {
    miRNAs <- unlist(strsplit(concordant_dmrs$miRNA[i], ";"))
    miRNAs <- miRNAs[miRNAs != "" & !is.na(miRNAs)]
    miRNAs <- trimws(miRNAs)
    
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
    
    mirna_rows <- grepl(paste0("(^|;)", mirna, "(;|$)"), concordant_dmrs$miRNA)
    num_cancers <- length(unique(unlist(strsplit(
      paste(concordant_dmrs$Cancers[mirna_rows], collapse = ";"), ";"))))
    
    miRNA_summary <- rbind(miRNA_summary,
                           data.frame(miRNA = mirna,
                                      Total_Concordant_DMRs = total_concordant,
                                      All_Hypermethylated_DMRs = all_hyper,
                                      All_Hypomethylated_DMRs = all_hypo,
                                      Num_Cancers = num_cancers,
                                      stringsAsFactors = FALSE))
  }
  
  miRNA_summary <- miRNA_summary[order(-miRNA_summary$Total_Concordant_DMRs), ]
  
  output_file1 <- file.path(results_dir, "concordant_miRNA_exact_DMRs_summaryFINAL.tsv")
  write.table(miRNA_summary, file = output_file1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  output_file2 <- file.path(results_dir, "concordant_exact_DMRs_detailedFINAL.tsv")
  write.table(concordant_dmrs, file = output_file2, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\nAnalysis complete.\n")
  cat(sprintf("miRNA-Cancer exact DMR matrix:   %s\n", output_matrix_file))
  cat(sprintf("Concordant miRNA exact DMRs summary:  %s\n", output_file1))
  cat(sprintf("Concordant exact DMRs detailed:  %s\n", output_file2))
  
  
  
  cancer_count_df <- as.data.frame(table(dmr_concordance$Num_Cancers))
  colnames(cancer_count_df) <- c("Num_Cancers", "Count")
  cancer_count_df$Num_Cancers <- as.numeric(as.character(cancer_count_df$Num_Cancers))
  
  p1 <- ggplot(cancer_count_df, aes(x = Num_Cancers, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +
    theme_minimal(base_size = 12) +
    labs(title = "Distribution of Exact Coordinate DMRs by Number of Cancers (2+)",
         x = "Number of Cancers",
         y = "Number of DMRs") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plot1_file <- file.path(results_dir, "exact_dmr_cancer_distributionFINAL.pdf")
  ggsave(plot1_file, p1, width = 10, height = 6)
  cat(sprintf("✓ Saved:  %s\n", plot1_file))
  
  if (length(top_mirnas) > 0) {
    top20_df <- data.frame(
      miRNA = names(top_mirnas),
      Total_DMRs = as.numeric(top_mirnas)
    )
    top20_df$miRNA <- factor(top20_df$miRNA, levels = top20_df$miRNA)
    
    p2 <- ggplot(top20_df, aes(x = miRNA, y = Total_DMRs)) +
      geom_bar(stat = "identity", fill = "coral", color = "black") +
      theme_minimal(base_size = 11) +
      labs(title = "Top 20 miRNAs by Exact Coordinate Concordant DMRs",
           x = "miRNA",
           y = "Total DMRs") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    plot2_file <- file.path(results_dir, "exact_top20_mirnas_dmrsFINAL.pdf")
    ggsave(plot2_file, p2, width = 12, height = 6)
    cat(sprintf("✓ Saved: %s\n", plot2_file))
  }
  
  chr_counts <- table(concordant_dmrs$Chromosome)
  chr_df <- data.frame(
    Chromosome = names(chr_counts),
    Count = as.numeric(chr_counts)
  )
  
  chr_df <- chr_df[!  chr_df$Chromosome %in% c("chrX", "chrY", "chrM"), ]
  
  if (nrow(chr_df) > 0) {
    chr_order <- paste0("chr", 1:22)
    chr_df$Chromosome <- factor(chr_df$Chromosome, 
                                levels = chr_order[chr_order %in% chr_df$Chromosome])
    chr_df <- chr_df[order(chr_df$Chromosome), ]
    
    p3 <- ggplot(chr_df, aes(x = Chromosome, y = Count)) +
      geom_bar(stat = "identity", fill = "mediumpurple", color = "black") +
      theme_minimal(base_size = 11) +
      labs(title = "Exact Coordinate Concordant DMRs Distribution Across Chromosomes",
           x = "Chromosome",
           y = "Number of DMRs") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    plot3_file <- file.path(results_dir, "exact_dmrs_by_chromosomeFINAL.pdf")
    ggsave(plot3_file, p3, width = 12, height = 6)
    cat(sprintf("✓ Saved: %s\n", plot3_file))
  }
  
  cancer_dmr_counts <- data.frame(
    Cancer = colnames(dmr_matrix),
    DMR_Count = colSums(dmr_matrix)
  )
  cancer_dmr_counts <- cancer_dmr_counts[order(-cancer_dmr_counts$DMR_Count), ]
  cancer_dmr_counts$Cancer <- factor(cancer_dmr_counts$Cancer, levels = cancer_dmr_counts$Cancer)
  
  p4 <- ggplot(cancer_dmr_counts, aes(x = Cancer, y = DMR_Count)) +
    geom_bar(stat = "identity", fill = "salmon", color = "black") +
    theme_minimal(base_size = 11) +
    labs(title = "Number of Exact Coordinate Concordant DMRs per Cancer Type",
         x = "Cancer Type",
         y = "Number of DMRs") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plot4_file <- file.path(results_dir, "exact_dmrs_per_cancerFINAL.pdf")
  ggsave(plot4_file, p4, width = 12, height = 6)
  cat(sprintf("✓ Saved: %s\n", plot4_file))
  
  if (nrow(dmr_matrix) > 0 && ncol(dmr_matrix) > 0) {
    top_n_mirnas <- 50
    if (nrow(dmr_matrix) > top_n_mirnas) {
      mirna_sums <- rowSums(dmr_matrix)
      top_mirnas_idx <- order(mirna_sums, decreasing = TRUE)[1:top_n_mirnas]
      dmr_matrix_top <- dmr_matrix[top_mirnas_idx, ]
    } else {
      dmr_matrix_top <- dmr_matrix
    }
    
    dmr_matrix_top <- dmr_matrix_top[rowSums(dmr_matrix_top) > 0, , drop = FALSE]
    dmr_matrix_top <- dmr_matrix_top[, colSums(dmr_matrix_top) > 0, drop = FALSE]
    
    if (nrow(dmr_matrix_top) > 1 && ncol(dmr_matrix_top) > 1) {
      plot5_file <- file.path(results_dir, "exact_miRNA_cancer_heatmapFINAL.pdf")
      pdf(plot5_file, width = 14, height = 12)
      
      pheatmap(dmr_matrix_top,
               color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize_row = 8,
               fontsize_col = 10,
               main = paste0("Exact Coordinate Concordant DMRs:   Top ", nrow(dmr_matrix_top), " miRNAs Across Cancer Types"),
               border_color = "grey60")
      
      dev.off()
      cat(sprintf("✓ Saved: %s\n", plot5_file))
    }
  }
  
  
} else {
  cat("No concordant DMRs found.\n")
}