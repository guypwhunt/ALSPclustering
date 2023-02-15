# libraries
library(shiny)
library(ggplot2)
library(shinydashboard)
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
library(shinycssloaders)
library(MASS)
library(stringr)
library(ggpubr)
library(rlang)
library(shinybusy)
library(DESeq2)
library(RNAAgeCalc)

# Functions
columnToRowNames <- function(df, index) {
  rownames(df) <- df[, index]
  df <- df[,-index]
  return(df)
}

# Set Seed
set.seed(1)

objects <- NULL
input <- NULL

# Data Objects
input$tissue_type <- c("brain", "blood", "muscle", "nerve")
input$tissue_type <- input$tissue_type[1]
input$population <- c("All Races", "North European/European")
input$population <- input$population[1]

# 1. preprocessing of uploaded dataset
# 1a. Load data (FRONT END - upload files for both of these)
dataset_raw_counts <- fread("data/zucca_raw_counts.txt") %>%
  as.data.frame() %>% columnToRowNames(1)


dataset_phenotype <- fread("data/zucca_phenotype.txt") %>%
  as.data.frame() %>% columnToRowNames(1)

input$phenotypicColumn <- colnames(dataset_phenotype)[1]
input$age_column <- colnames(dataset_phenotype)[2]


# 1b. DESeq2 processing (FRONT END - should be raw counts only so emphasise and also should be in ensembl format)
# load data
objects$geneExpressionDf <- dataset_raw_counts
objects$phenotypicDf <- dataset_phenotype


###### Workflow #####

# create DEseq2 and normalise with size factors
dds_genes <-
  DESeqDataSetFromMatrix(
    countData = objects$geneExpressionDf,
    colData = objects$phenotypicDf,
    design = ~ 1
  ) %>%
  estimateSizeFactors()

# remove counts less than 10 in more than 5 individuals
idx <- rowSums(counts(dds_genes, normalized = TRUE) >= 5) >= 10
dds_genes <- dds_genes[idx, ]


# perform VSD stabilisation on whole matrix
dataset_vsd <- vst(dds_genes, blind = TRUE) %>%
  assay()

# 2. Cluster probability on uploaded dataset (dataset_vsd)
# 2a. load in the informative KCL brainbank matrix
brainbank_informative <-
  fread("data/brainbankinform_VSDremXY794genes.txt") %>%
  as.data.frame()  %>% columnToRowNames(1) %>% t()

# 2b. Find the intersection between brainbank_informative and the dataset and create subsets of the two datasets with these
dataset_informative_shared <-
  dataset_vsd[rownames(dataset_vsd) %in% rownames(brainbank_informative), ]
brainbank_informative_shared <-
  brainbank_informative[rownames(brainbank_informative) %in% rownames(dataset_informative_shared), ]
brainbank_t <- data.frame(t(brainbank_informative))


# 2c. Run LDA on the subsets - which format do they have to be in? transpose back?!
# assign the shared datasets to variables
brainbank_LDA <- data.frame(t(brainbank_informative_shared))
brainbank_LDA$Cluster <- brainbank_t$Cluster
dataset_LDA <- data.frame(t(dataset_informative_shared))

# extract number of informative genes as a variable
informative_gene_no <- length(rownames(dataset_informative_shared))

# train an LDA model using BrainBank data and cluster assignments
brainbank_LDA[1:informative_gene_no] <-
  scale(brainbank_LDA[1:informative_gene_no])
train <- brainbank_LDA
model <- lda(Cluster ~ ., data = train)

# predict dataset assignments using brainbank training set
# str(dataset_LDA)
dataset_LDA[1:informative_gene_no] <-
  scale(dataset_LDA[1:informative_gene_no])
test <- dataset_LDA
predicted <- predict(model, test)

lda_class <- predicted$x
plot_dat <- cbind(test, lda_class)
lda_plot <-
  ggplot(plot_dat, aes(LD1, LD2)) + geom_point(aes(color = predicted$class), size =
                                                 0.7) + theme_bw() + theme(legend.position = "right") + xlab("Linear Discriminant 1 (LD1)") + ylab("Linear Discriminant 2 (LD2)") + stat_ellipse(aes(group = predicted$class, color = predicted$class))
lda_plot


colnames(predicted$posterior) <-
  c(
    "Probability of being in Cluster 1",
    "Probability of being in Cluster 2",
    "Probability of being in Cluster 3"
  )
class(predicted$class)

objects$ldaResults <- cbind(predicted$x, predicted$posterior) %>%
  as.data.frame()

objects$ldaResults[, "Cluster Assignment"] <-
  as.data.frame(predicted$class)

objects$ldaResults <-
  cbind(objects$phenotypicDf, objects$ldaResults)

objects$ldaResultsMachineFriendly <- objects$ldaResults

colnames(objects$ldaResultsMachineFriendly) <-
  str_replace_all(
    string = colnames(objects$ldaResultsMachineFriendly),
    pattern = " ",
    repl = ""
  )

phenotypicColumMachineFriendly <-
  str_replace_all(
    string = input$phenotypicColumn,
    pattern = " ",
    repl = ""
  )

plot_ly(data = objects$ldaResultsMachineFriendly) %>%
  add_trace(
    x = ~ LD1,
    y = ~ LD2,
    color = ~ ClusterAssignment,
    text = paste0(
      phenotypicColumMachineFriendly,
      " : ",
      objects$ldaResultsMachineFriendly[, phenotypicColumMachineFriendly]
    ),
    type = 'scatter',
    mode = 'markers',
    legendgroup = "Cluster",
    showlegend = T
  ) %>%
  add_annotations(
    text = "Clusters",
    xref = "paper",
    yref = "paper",
    x = 1.02,
    xanchor = "left",
    y = 0.9,
    yanchor = "bottom",
    # Same y as legend below
    legendtitle = TRUE,
    showarrow = FALSE
  ) %>%
  #Increase distance between groups in Legend
  layout(legend = list(
    tracegroupgap = 50,
    y = 0.9,
    yanchor = "top"
  ))



######### Transcriptimal analysis

# run age analysis and add to the phenotype and/or results table
if(input$population == "All Races") {
  population <- "all"
} else {
  population <- "caucasian"
}

res <- predict_age(exprdata = objects$geneExpressionDf, tissue = input$tissue_type, exprtype = "counts", idtype = "ENSEMBL", stype = population)

RNAAge <- res$RNAAge
RNAAge_accel <- res$RNAAge - objects$ldaResultsMachineFriendly[ , input$age_column]

objects$ldaResultsMachineFriendly$TranscriptionalAge <- RNAAge
objects$ldaResultsMachineFriendly$TranscriptionalAgeAcceleration <- RNAAge_accel

objects$ldaResults[, "Transcriptional Age"] <- RNAAge
objects$ldaResults[, "Transcriptional Age Acceleration"] <- RNAAge_accel


######
# load in data
input$phenotypicComparisonColumn <- "Age_Onset"
input$phenotypicComparisonCovariants <- "Sex"

phenotypicComparisonColumn <- str_replace_all(string = input$phenotypicComparisonColumn,
                                              pattern = " ",
                                              repl = "")

phenotypicComparisonCovariants <- str_replace_all(string = input$phenotypicComparisonCovariants,
                                              pattern = " ",
                                              repl = "")

phenotypic_variable <- objects$ldaResultsMachineFriendly[, phenotypicComparisonColumn]
 # whatever phenotype is chosen

columnsToExclude <- c("LD1", "LD2", "ProbabilityofbeinginCluster1",
                      "ProbabilityofbeinginCluster2", "ProbabilityofbeinginCluster3"#,
                      # "TranscriptionalAge", "TranscriptionalAgeAcceleration")
                      )


for(column in colnames(objects$ldaResultsMachineFriendly)) {
  if(!column %in% c("LD1", "LD2", "ProbabilityofbeinginCluster1",
                    "ProbabilityofbeinginCluster2", "ProbabilityofbeinginCluster3",
                    "TranscriptionalAge", "TranscriptionalAgeAcceleration")){

    typeOfColumn <- typeof(objects$ldaResultsMachineFriendly[, column])

    isColumnFactor <- !is.factor(objects$ldaResultsMachineFriendly[, column])

    if(typeOfColumn == "integer" & isColumnFactor){
      if(max(objects$ldaResultsMachineFriendly[, column]) < 6){
        if(max(objects$ldaResultsMachineFriendly[, column]) > -1){
          objects$ldaResultsMachineFriendly[, column] <- as.factor(objects$ldaResultsMachineFriendly[, column])
        }
      }
    } else if(typeOfColumn != "double" & isColumnFactor){
      objects$ldaResultsMachineFriendly[, column] <- as.factor(objects$ldaResultsMachineFriendly[, column])
    }
  }
}

# descriptive statistics for the variable
objects$phenotypicComparisonSummaryResults <-
  group_by(
    objects$ldaResultsMachineFriendly,
    objects$ldaResultsMachineFriendly$ClusterAssignment
  ) %>%
  summarise(
    count = n(),
    mean = mean((!!sym(phenotypicComparisonColumn)), na.rm = TRUE),
    sd = sd((!!sym(phenotypicComparisonColumn)), na.rm = TRUE),
    median = median((!!sym(phenotypicComparisonColumn)), na.rm = TRUE),
    IQR = IQR((!!sym(phenotypicComparisonColumn)), na.rm = TRUE)
  ) %>% as.data.frame()

colnames(objects$phenotypicComparisonSummaryResults)[1] <- "cluster"

# normality of each variable assessed with shapiro wilk
normality <- shapiro.test(objects$ldaResultsMachineFriendly[, phenotypicComparisonColumn])

# variable is log-transformed before test if not normal
if (normality$p.value < 0.05){
  phenotypic_variable_transformed <- log(phenotypic_variable)
} else {
  phenotypic_variable_transformed <- phenotypic_variable
}

equation <- "phenotypic_variable_transformed ~ objects$ldaResultsMachineFriendly$ClusterAssignment"

covariantsFlag <- length(phenotypicComparisonCovariants) > 1
if(isFALSE(covariantsFlag)) {
  covariantsFlag <- phenotypicComparisonCovariants != "NA"
}

if(covariantsFlag) {
  for (covariant in phenotypicComparisonCovariants) {
    equation <- paste0(equation, " + ", covariant)
  }
}

res <- aov(
  as.formula(equation),
  data = objects$ldaResultsMachineFriendly
)

# one-way ANOVA
anova_p_value <- summary(res)[[1]][["Pr(>F)"]][1]

# tukey's test
pairwise_res <- TukeyHSD(res)
objects$phenotypicComparisonPairwiseResults <- as.data.frame(pairwise_res$`objects$ldaResultsMachineFriendly$ClusterAssignment`)

tukey_p_values <- objects$phenotypicComparisonPairwiseResults[,4]

columnNamesOrder <- colnames(objects$phenotypicComparisonPairwiseResults)
objects$phenotypicComparisonPairwiseResults$comparison <- rownames(objects$phenotypicComparisonPairwiseResults)
objects$phenotypicComparisonPairwiseResults <- objects$phenotypicComparisonPairwiseResults[,
                                                                                           c("comparison", columnNamesOrder)]

colnames(objects$phenotypicComparisonPairwiseResults) <- c("comparison", "difference", "lower confidence interval limit", "upper confidence interval limit", "adjusted p-value")

tukey_first_group <- c("2", "3", "3")
tukey_second_group <- c("1", "1", "2")

# put the p values in a handy table for the plotting of the p values
p_value_results <- data.frame(x1 = tukey_first_group, x2= tukey_second_group, x3 = tukey_p_values)
colnames(p_value_results) <- c("group1", "group2", "p.adj")

# plot (unsure how to do the "name of phenotypic variable" thing as it wont let you do phenotypic_data$Cluster for example)
pheno_plot <- ggboxplot(objects$ldaResultsMachineFriendly, x = "ClusterAssignment", y = phenotypicComparisonColumn, color = "ClusterAssignment", title = sprintf("ANOVA p-value: %s", anova_p_value))


plot_ly(data = objects$ldaResultsMachineFriendly) %>%
  add_trace(
    y = objects$ldaResultsMachineFriendly[, phenotypicComparisonColumn],
    x = ~ ClusterAssignment,
    color = ~ ClusterAssignment,
    type = 'box',
    legendgroup = "Cluster",
    showlegend = T
  ) %>%
  add_annotations(
    text = "Clusters",
    xref = "paper",
    yref = "paper",
    x = 1.02,
    xanchor = "left",
    y = 0.9,
    yanchor = "bottom",
    # Same y as legend below
    legendtitle = TRUE,
    showarrow = FALSE
  ) %>%
  #Increase distance between groups in Legend
  layout(
    legend = list(
      tracegroupgap = 50,
      y = 0.9,
      yanchor = "top"
    ),
    xaxis = list(title = 'Cluster Assignment'),
    yaxis = list(title = input$phenotypicComparisonColumn)
  )
