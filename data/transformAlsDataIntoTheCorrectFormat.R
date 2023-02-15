library(biomaRt)
library(ggplot2)
library(knitr)
library(rmarkdown)
library(ggthemes)
library(plotly)
library(DT)
library(data.table)
library(dplyr)
library(shinythemes)
library(RColorBrewer)
library(markdown)
library(tibble)
library(stringr)


# final data  ################################################################
# exprs_data <- readRDS("./data/gxprs_base_pheno_alt_names.rds")
# exprsData <- as.data.frame(exprs_data)
# write.csv(head(exprsData), "exprsData.csv")

# raw data  ################################################################
# raw_data <- readRDS("./data/gxprs_raw_base_pheno_alt_names_md5.rds")
# rawData <- as.data.frame(raw_data)
# write.csv(head(rawData), "rawData.csv")

# raw data niave process: bk, log2, rsn  ################################################################
# raw_data_bklogrsn <- readRDS("./data/gxprs_bklog2rsnCombat_347_31431.rds")
# rawDataBklogrsn <- as.data.frame(raw_data_bklogrsn)
# head(rawDataBklogrsn)
# write.csv(head(rawDataBklogrsn), "rawDataBklogrsn.csv")

# limma results  ################################################################
# limma_df <- readRDS("./data/MRC_LBB_Limma.rds")
# limmaDf <- as.data.frame(limma_df)
# write.csv(head(limmaDf), "limmaDf.csv")

# limma_df_raw <- readRDS("./data/MRC_LBB_Raw_Limma.rds")
# limmaDfRaw <- as.data.frame(limma_df_raw)
# write.csv(head(limmaDfRaw), "limmaDfRaw.csv")

########## ALS data
bbExperimentDesign <- read.table("data/Matrices/BBsamples.design.updated.sv1.171subjects.txt", sep = '\t', header = TRUE)
bbExpressionMatrix <- read.table("data/Matrices/BBNormMatrix.txt", sep = '\t', header = TRUE)
bbDEResults <- read.table("data/Matrices/BBgenesDEresOriginal171.txt", sep = '\t', header = TRUE)

taExperimentDesign <- read.table("data/Matrices/TargetAlssamples.design.updated.site.sv1.234subjects.txt", sep = '\t', header = TRUE)
taExpressionMatrix <- read.table("data/Matrices/TargetALSNormMatrix.txt", sep = '\t', header = TRUE)
taDEResults <- read.table("data/Matrices/TargetAlsgenesDEres234subjectsTargetALS_site copy.txt", sep = '\t', header = TRUE)

# Update the expression datasets column names
colnames(taExpressionMatrix) <- c("gene", head(colnames(taExpressionMatrix), -1))

# order the expression data and results by ensemble gene ID
taExpressionMatrix <- taExpressionMatrix[order(taExpressionMatrix[,"gene"]),]
taDEResults[,"gene"] <- rownames(taDEResults)
taDEResults <- taDEResults[taDEResults[,"gene"] %in% taExpressionMatrix[,"gene"],]
taDEResults <- taDEResults[order(taDEResults[,"gene"]),]

bbExpressionMatrix[,"gene"] <- rownames(bbExpressionMatrix)
bbExpressionMatrix <- bbExpressionMatrix[order(bbExpressionMatrix[,"gene"]),]
bbDEResults[,"gene"] <- rownames(bbDEResults)
bbDEResults <- bbDEResults[order(bbDEResults[,"gene"]),]

# Check that the expression data and DE results are in the same order
all(taDEResults[,"gene"] == taExpressionMatrix[,"gene"])
all(bbExpressionMatrix[,"gene"] == bbDEResults[,"gene"])

### Identify gene symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- taExpressionMatrix[,1]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "external_gene_name", "description", "entrezgene_id"),values=genes,mart= mart)

G_list <- G_list[ !duplicated(G_list[, c("ensembl_gene_id", "external_gene_name")], fromLast=T),]

# Merge datasets
taExpressionMatrix <- merge(taExpressionMatrix, G_list[, c("ensembl_gene_id", "external_gene_name")], by.x = "gene" , by.y = "ensembl_gene_id", all.x = TRUE)
taDEResults <- merge(taDEResults, G_list[, c("ensembl_gene_id", "external_gene_name", "entrezgene_id")], by.x = "gene" , by.y = "ensembl_gene_id", all.x = TRUE)

# Rows without gene symbols
taExpressionMatrix[is.na(taExpressionMatrix[,"external_gene_name"]), 1]

### Identify gene symbols
genes <- bbExpressionMatrix[,"gene"]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "external_gene_name", "description", "entrezgene_id"),values=genes,mart= mart)
G_list <- G_list[ !duplicated(G_list[, c("ensembl_gene_id", "external_gene_name")], fromLast=T),]

# Merge datasets
bbExpressionMatrix <- merge(bbExpressionMatrix, G_list[, c("ensembl_gene_id", "external_gene_name")], by.x = "gene" , by.y = "ensembl_gene_id", all.x = TRUE)
bbDEResults <- merge(bbDEResults, G_list[, c("ensembl_gene_id", "external_gene_name", "entrezgene_id")], by.x = "gene" , by.y = "ensembl_gene_id", all.x = TRUE)

# Identify Cases and Controls
bbExperimentDesign[bbExperimentDesign[, "Status"] == 1, "Phenotype"] <- " Control"
bbExperimentDesign[bbExperimentDesign[, "Status"] == 2, "Phenotype"] <- "ALS"

taExperimentDesign[taExperimentDesign[, "Status"] == 1, "Phenotype"] <- " Control"
taExperimentDesign[taExperimentDesign[, "Status"] == 2, "Phenotype"] <- "ALS"

# Fill in missing genes with ensemble IDs
bbExpressionMatrix[is.na(bbExpressionMatrix[,"external_gene_name"]), "external_gene_name"] <- bbExpressionMatrix[is.na(bbExpressionMatrix[,"external_gene_name"]), "gene"]
taExpressionMatrix[is.na(taExpressionMatrix[,"external_gene_name"]), "external_gene_name"] <- taExpressionMatrix[is.na(taExpressionMatrix[,"external_gene_name"]), "gene"]
bbExpressionMatrix[bbExpressionMatrix[,"external_gene_name"] == "", "external_gene_name"] <- bbExpressionMatrix[bbExpressionMatrix[,"external_gene_name"] == "", "gene"]
taExpressionMatrix[taExpressionMatrix[,"external_gene_name"] == "", "external_gene_name"] <- taExpressionMatrix[taExpressionMatrix[,"external_gene_name"] == "", "gene"]

bbDEResults[is.na(bbDEResults[,"external_gene_name"]), "external_gene_name"] <- bbDEResults[is.na(bbDEResults[,"external_gene_name"]), "gene"]
taDEResults[is.na(taDEResults[,"external_gene_name"]), "external_gene_name"] <- taDEResults[is.na(taDEResults[,"external_gene_name"]), "gene"]
bbDEResults[bbDEResults[,"external_gene_name"] == "", "external_gene_name"] <- bbDEResults[bbDEResults[,"external_gene_name"] == "", "gene"]
taDEResults[taDEResults[,"external_gene_name"] == "", "external_gene_name"] <- taDEResults[taDEResults[,"external_gene_name"] == "", "gene"]


# Update Duplicate Genes to include Ensemble ID
taExpressionMatrix[duplicated(taExpressionMatrix[, "external_gene_name"]), "external_gene_name"] <- paste0(taExpressionMatrix[duplicated(taExpressionMatrix[, "external_gene_name"]), "external_gene_name"], " (", taExpressionMatrix[duplicated(taExpressionMatrix[, "external_gene_name"]), "gene"], ")")
bbExpressionMatrix[duplicated(bbExpressionMatrix[, "external_gene_name"]), "external_gene_name"] <- paste0(bbExpressionMatrix[duplicated(bbExpressionMatrix[, "external_gene_name"]), "external_gene_name"], " (", bbExpressionMatrix[duplicated(bbExpressionMatrix[, "external_gene_name"]), "gene"], ")")

taDEResults[duplicated(taDEResults[, "external_gene_name"]), "external_gene_name"] <- paste0(taDEResults[duplicated(taDEResults[, "external_gene_name"]), "external_gene_name"], " (", taDEResults[duplicated(taDEResults[, "external_gene_name"]), "gene"], ")")
bbDEResults[duplicated(bbDEResults[, "external_gene_name"]), "external_gene_name"] <- paste0(bbDEResults[duplicated(bbDEResults[, "external_gene_name"]), "external_gene_name"], " (", bbDEResults[duplicated(bbDEResults[, "external_gene_name"]), "gene"], ")")

# Update row names
rownames(bbExpressionMatrix) <- bbExpressionMatrix[,"external_gene_name"]
rownames(taExpressionMatrix) <- taExpressionMatrix[,"external_gene_name"]
rownames(taDEResults) <- taDEResults[,"external_gene_name"]
rownames(bbDEResults) <- bbDEResults[,"external_gene_name"]

bbExpressionMatrix <- bbExpressionMatrix[, tail(head(colnames(bbExpressionMatrix), -1), -1)]
taExpressionMatrix <- taExpressionMatrix[, tail(head(colnames(taExpressionMatrix), -1), -1)]

# Transpose expression data
transposedBbExpressionMatrix <- t(bbExpressionMatrix)
transposedTaExpressionMatrix <- t(taExpressionMatrix)

transposedBbExpressionMatrix <- as.data.frame(transposedBbExpressionMatrix)
transposedTaExpressionMatrix <- as.data.frame(transposedTaExpressionMatrix)


# Add phenotypic data
transposedBbExpressionMatrix <- transposedBbExpressionMatrix[order(rownames(transposedBbExpressionMatrix)), ]
bbExperimentDesign <- bbExperimentDesign[order(rownames(bbExperimentDesign)), ]

transposedTaExpressionMatrix <- transposedTaExpressionMatrix[order(rownames(transposedTaExpressionMatrix)), ]
taExperimentDesign <- taExperimentDesign[order(rownames(taExperimentDesign)), ]

colNamestransposedBbExpressionMatrix <- colnames(transposedBbExpressionMatrix)
colNamestransposedTaExpressionMatrix <- colnames(transposedTaExpressionMatrix)

nrow(transposedBbExpressionMatrix)
nrow(bbExperimentDesign)

transposedBbExpressionMatrix <- cbind(transposedBbExpressionMatrix, bbExperimentDesign$Phenotype)
transposedTaExpressionMatrix <- cbind(transposedTaExpressionMatrix, taExperimentDesign$Phenotype)

colnames(transposedBbExpressionMatrix)[ncol(transposedBbExpressionMatrix)] <- "Phenotype"
colnames(transposedTaExpressionMatrix)[ncol(transposedTaExpressionMatrix)] <- "Phenotype"

####

transposedBbExpressionMatrix["Dataset"] <- rep("KCL Brain Bank", nrow(transposedBbExpressionMatrix))
transposedTaExpressionMatrix["Dataset"] <- rep("Target ALS", nrow(transposedTaExpressionMatrix))

transposedBbExpressionMatrix["Tissue"] <- rep("Motor Cortex", nrow(transposedBbExpressionMatrix))
transposedTaExpressionMatrix["Tissue"] <- rep("Motor Cortex", nrow(transposedTaExpressionMatrix))



transposedTaExpressionMatrix <- transposedTaExpressionMatrix[, c("Phenotype", "Dataset", "Tissue", colNamestransposedTaExpressionMatrix)]
transposedBbExpressionMatrix <- transposedBbExpressionMatrix[, c("Phenotype", "Dataset", "Tissue", colNamestransposedBbExpressionMatrix)]

# Save datasets
#write.csv(transposedBbExpressionMatrix, "data/brainBainExpressionData.csv")
#write.csv(transposedTaExpressionMatrix, "data/targetAlsExpressionData.csv")

# Load datasets
#transposedTaExpressionMatrix <- read.csv("data/brainBainExpressionData.csv")
#transposedBbExpressionMatrix <- read.csv("data/targetAlsExpressionData.csv")

# Update exression datasets to have the same columns
colnamesToAddToBb <- setdiff(colnames(transposedTaExpressionMatrix), colnames(transposedBbExpressionMatrix))
colnamesToAddToTa <- setdiff(colnames(transposedBbExpressionMatrix), colnames(transposedTaExpressionMatrix))

transposedTaExpressionMatrix <- as.data.frame(transposedTaExpressionMatrix)
transposedBbExpressionMatrix <- as.data.frame(transposedBbExpressionMatrix)

transposedTaExpressionMatrix[colnamesToAddToTa] <- rep(NA, nrow(transposedTaExpressionMatrix))
transposedBbExpressionMatrix[colnamesToAddToBb] <- rep(NA, nrow(transposedBbExpressionMatrix))

# Reorder the columns
transposedTaExpressionMatrix <- transposedTaExpressionMatrix[,colnames(transposedBbExpressionMatrix)]

# Merge datasets
combinedExpressionData <- rbind(transposedBbExpressionMatrix, transposedTaExpressionMatrix)

# Save datast
#write.csv(combinedExpressionData, "data/combinedExpressionData.csv")

# Convert to tibble and save
combinedExpressionData <- as_tibble(combinedExpressionData)
saveRDS(combinedExpressionData, "data/combinedExpressionData.rds")

##### Format Differential Gene Expression ResultS
head(bbDEResults)
head(taDEResults)

bbDEResults["Dataset"] <- "KCL Brain Bank"
taDEResults["Dataset"] <- "Target ALS"

bbDEResults["Tissue"] <- "Motor Cortex"
taDEResults["Tissue"] <- "Motor Cortex"

colnames(bbDEResults)[which(names(bbDEResults) == "external_gene_name")] <- "Gene"
colnames(taDEResults)[which(names(taDEResults) == "external_gene_name")] <- "Gene"

colnames(bbDEResults)[which(names(bbDEResults) == "entrezgene_id")] <- "entrez_id"
colnames(taDEResults)[which(names(taDEResults) == "entrezgene_id")] <- "entrez_id"

colnames(bbDEResults)[which(names(bbDEResults) == "padj")] <- "adj.P.Val"
colnames(taDEResults)[which(names(taDEResults) == "padj")] <- "adj.P.Val"

colnames(bbDEResults)[which(names(bbDEResults) == "log2FoldChange")] <- "logFC"
colnames(taDEResults)[which(names(taDEResults) == "log2FoldChange")] <- "logFC"

columns <- c("Dataset", "Gene", "Tissue", "gene", "entrez_id", "logFC", "pvalue", "adj.P.Val")
bbDEResults <- bbDEResults[, columns]
taDEResults <- taDEResults[, columns]

combinedResults <- rbind(bbDEResults, taDEResults)

# Save datast
#write.csv(combinedResults, "data/combinedResults.csv")
# combinedResults <- read.csv("data/combinedResults.csv", row.names = 1)
head(combinedResults)

# Convert to tibble and save
combinedResults <- as.data.frame(combinedResults)

combinedResults <- combinedResults %>%
  mutate(Log_Fold_Change = round(logFC, 2)) %>%
  mutate(P_Value = signif(pvalue, 3)) %>%
  mutate(FDR_P_Value = signif(adj.P.Val, 3)) %>%
  mutate(Gene_Symbol = gsub(" .*$", "", Gene)) %>%
  mutate(
    Gene_Symbol_URL = paste(
      '<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=',
      Gene_Symbol,
      ' "target="_blank"',
      '">',
      Gene_Symbol,
      '</a>',
      sep = ""
    )
  ) %>%
  mutate(Ensembl_ID = gene) %>%
  mutate(
    Ensembl_ID_URL = paste(
      '<a href="http://www.ensembl.org/id/',
      Ensembl_ID,
      ' "target="_blank"',
      '">',
      Ensembl_ID,
      '</a>',
      sep = ""
    )
  ) %>%
  mutate(Entrez_ID = entrez_id) %>%
  mutate(
    Entrez_ID_URL = paste(
      '<a href="https://www.ncbi.nlm.nih.gov/gene/',
      entrez_id,
      ' "target="_blank"',
      '">',
      entrez_id,
      '</a>',
      sep = ""
    )
  ) %>%
  dplyr::select(
    Dataset,
    Gene,
    Tissue,
    Gene_Symbol,
    Gene_Symbol_URL,
    Entrez_ID,
    Entrez_ID_URL,
    Ensembl_ID,
    Ensembl_ID_URL,
    Log_Fold_Change,
    P_Value,
    FDR_P_Value
  )

rownames(combinedResults) <- seq(nrow(combinedResults))

head(combinedResults)

saveRDS(combinedResults, "data/combinedResults.rds")

combinedResults <- readRDS("data/combinedResults.rds")
head(combinedResults)

combinedExpressionData <- readRDS("data/combinedExpressionData.rds")
combinedExpressionData[1:5,1:5]

combinedExpressionData[combinedExpressionData$Phenotype == "Case", "Phenotype"] <- "ALS"
combinedExpressionData[combinedExpressionData$Phenotype == "Control", "Phenotype"] <- " Control"

saveRDS(combinedExpressionData, "data/combinedExpressionData1.rds")

