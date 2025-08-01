library(minfi)                                                                  #Methylation data analysis
library(limma)                                                                  #linear models for diff analysis
library(ggpubr)                                                                  
library(DMRcate)                                                                #Differentially methylated regions detection
library(reshape2)                                 
library(missMethyl)                                                             #statistical analysis for methylation data
library(doParallel)                                                             #parallel processing
library(TCGAbiolinks)                                                           #tcga Data downloading
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(sesame)
library(sesameData)
library(DMRcatedata)

setwd("C:/Users/konst/OneDrive/Documents/MethylationAnalyis")

query_met <- GDCquery(project = "TCGA-UVM",                                    #the dataset
                      data.category = "DNA Methylation",                        #data category
                      data.type = "Methylation Beta Value",                     #data type
                      platform = "Illumina Human Methylation 450")
head(query_met$results[[1]]) 

GDCdownload(query_met, directory = "GDCdata_beta_values", files.per.chunk = 20)                       #download data 
data.met <- GDCprepare(query_met, directory = "GDCdata_beta_values/")           #load them into R, data.met is an object containing methylation data and metadata



paired.sample.ids <- data.met[, duplicated(data.met$patient)]$patient           #identify patients tha have both tumor and normal samples
data.met <- data.met[, !((data.met$patient %in% paired.sample.ids) & (data.met$definition == "Primary solid Tumor"))] #remove from data.met the primary tumors of the paired patients
data.met <- data.met[, !data.met$definition == "Recurrent Solid Tumor"]         #remore recurrent tumors
                                                                                # !!! Now the data set has: normal samples from patients with both tumor and normal samples and primary tumors from patients that only have tumors!!!!

data.met <- data.met[!data.met@rowRanges$probeType == "rs", ]                   # remove probes that associated with SNPs because we may confuse gentic variantes with epigenetic changes (methylation)
data.met <- data.met[!data.met@rowRanges@seqnames %in% c("chrX", "chrY")]       #remove probes on sex chroms -> better quality of data


data.met <- data.met[rowSums(is.na(assay(data.met))) == 0, ]                    #also remove probes with missing values(NAs)


clinical_data <- data.met@colData                                               #extract metadata so in Clinical_data i have clean methylation data 
cat("Number of paired samples:", nrow(clinical_data), "\n")
cat("Average age of paired samples:", mean(as.numeric(clinical_data$age_at_diagnosis), na.rm = TRUE), "\n")
cat("Sex distribution of paired samples:\n")
print(table(clinical_data$gender))


chisq.test(table(clinical_data[c("gender", "definition")]))

data.met$definition <- gsub("Primary solid Tumor", "Primary_solid_Tumor", data.met$definition) #replace the gaps in labels with _ 
data.met$definition <- gsub("Solid Tissue Normal", "Solid_Tissue_Normal", data.met$definition)

# box plots
df <- data.frame("Sample.mean" = colMeans(assay(data.met), na.rm = TRUE), "Groups" = data.met$definition)

ggboxplot(data = df, y = "Sample.mean", x = "Groups", color = "Groups", add = "jitter",
          ylab = "Mean DNA methylation (beta-values)", xlab = "",
          title = "Mean DNA Methylation Levels Across Groups") + 
          stat_compare_means()
ggsave(filename = "mean_methylation.png", bg = "white", width = 12, height = 10, dpi = 300)

#density plot (distribution of beta values across all probes/samples)
beta_values <- as.data.frame(t(assay(data.met)))
beta_values$Group <- as.factor(data.met$definition)


ggplot(melt(beta_values), aes(x = value, color = Group, fill = Group)) +
      geom_density(alpha = 0.3, size = 1.2) +
      labs(title = "Density Plot of Beta Values", x = "Beta Value", y = "Density") +
      theme_minimal(base_size = 14) +
      scale_color_brewer(palette = "Set2") +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "top") 
ggsave(filename = "density_plot.png", bg = "white", width = 12, height = 10, dpi = 300)

#MDS Plot (check batch effects or global methylation differences)
beta_values2 <- t(assay(data.met))
dist_matrix <- dist(beta_values2)
mds_result <- cmdscale(dist_matrix, k = 2)
mds_df <- as.data.frame(mds_result)
colnames(mds_df) <- c("MDS1", "MDS2")
mds_df$Group <- as.factor(data.met$definition)


ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Group)) +
      geom_point(size = 2, alpha = 0.7) +
      labs(title = "MDS Plot of Methylation Data", x = "MDS Dimension 1", y = "MDS Dimension 2") +
      theme_minimal(base_size = 14) +
      scale_color_brewer(palette = "Set1") +
      theme(legend.position = "top") 
ggsave(filename = "mds_plot.png", bg = "white", width = 12, height = 10, dpi = 300)


############################################################ DIFFERENTIAL METHYLATION ANALYSIS################################################
dmc <- TCGAanalyze_DMC(data = data.met,                                         #dmc has all the differentially methylated CpG island regions stored in it
                      groupCol = "definition",
                      group1 = "Primary_solid_Tumor",
                      group2 = "Solid_Tissue_Normal",
                      p.cut = 0.05,                                             #p-value=0.05
                      diffmean.cut = 0.15,                                      #It ensures that only CpG sites with a mean methylation difference ≥ 0.15 (15%)(biologically significant) between the two groups are retained
                      save = FALSE,
                      legend = "State",
                      plot.filename = "UVM_vs_Normal_metvolcano.png",
                      cores = 1)


dmc_df <- as.data.frame(dmc)
dmc_df$MethylationStatus <- ifelse(dmc_df$p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal < 0.05 &                                       #classification of DM CpGs 
                                     dmc_df$mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal > 0.15, "Hypermethylated",               # if p_val<0.05 & mean>0.15 = hypermethylated                           
                                   ifelse(dmc_df$p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal < 0.05 &                                # if p_val<0.05 & mean<-0.15 = hypomethylated
                                            dmc_df$mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal < -0.15, "Hypomethylated", "Not Significant")) #else insignificant 

#volcano plot of the most significant CpGs
ggplot(dmc_df, aes(x = mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal, y = -log10(p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal), color = MethylationStatus)) +
       geom_point(alpha = 0.8, size = 1.5) +
       scale_color_manual(values = c("Hypermethylated" = "#8ABFA3", "Hypomethylated" = "#F95454", "Not Significant" = "#B7B7B7")) +
       labs(title = "Volcano Plot of Differential Methylation", 
            x = "Mean Difference (Tumor vs. Normal)", 
            y = "-log10(adjusted p-value)") +
       theme_minimal(base_size = 14) +
       theme(legend.position = "top") +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
       geom_vline(xintercept = c(-0.15, 0.15), linetype = "dashed")
ggsave(filename = "volcano_plot.png", bg = "white", width = 12, height = 10, dpi = 300)


pheno <- data.met$definition                                                                                            #extract the phenotype (sample group labels) from the data.met object.
design.matrix <- model.matrix(~0 + pheno)                                                                               #Create a design matrix for linear modeling without an intercept (~0 +).
colnames(design.matrix) <- gsub("pheno", "", colnames(design.matrix))                                                   # Cleans up column names by removing the "pheno" prefix.                                                                                                                                     
cont.matrix <- makeContrasts(Cancer.vs.Normal = Primary_solid_Tumor - Solid_Tissue_Normal, levels = design.matrix)      #Defines the contrast for differential analysis: tumor vs. normal.  

GRset <- GenomicRatioSet(Beta = assay(data.met), gr = rowRanges(data.met), colData = colData(data.met))                 #Creates a GenomicRatioSet object from the methylation data.
annotation(GRset) <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")                              #Adds platform and genome annotation to the GRset object.


myAnnotation <- DMRcate::cpg.annotate(object = GRset, arraytype = "450K", datatype = "array", what = "Beta",           #Annotates CpG sites for differential methylation analysis.
                                      analysis.type = "differential", design = design.matrix, contrasts = TRUE,      
                                      cont.matrix = cont.matrix, coef = "Cancer.vs.Normal")

dmrcoutput <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 0.05)                                             #Identifies differentially methylated regions (DMRs).
significantRegions <- extractRanges(dmrcoutput)                                                                       #lambda = 1000: smoothing window size (in bp).
significantRegions <- as.data.frame(significantRegions)                                                                                                                     #C = 2: scaling factor for bandwidth.
                                                                                                                      #pcutoff = 0.05: significance threshold.

write.csv(dmc_df, "UVM_vs_Normal_dmc.results.csv", row.names = FALSE)
write.csv(significantRegions, "UVM_vs_Normal_dmr.results.csv", row.names = FALSE)

rm(list = setdiff(ls(), c("data.met", "clinical_data", "significantRegions", "dmc_df")))

