## app.R #####
try(library(shiny))
try(library(ggplot2))
try(library(shinydashboard))
try(library(knitr))
try(library(rmarkdown))
try(library(ggthemes))
try(library(plotly))
try(library(DT))
try(library(data.table))
try(library(dplyr))
try(library(shinythemes))
try(library(RColorBrewer))
try(library(markdown))
try(library(shinycssloaders))
try(library(MASS))
try(library(stringr))
try(library(ggpubr))
try(library(rlang))
try(library(shinybusy))
try(library(DESeq2))
try(library(RNAAgeCalc))

columnToRowNames <- function(df, index) {
  rownames(df) <- df[, index]
  df <- df[,-index]
  return(df)
}

# Set Seed
set.seed(1)

# Data Objects


# UI ####
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "ALS Phenotypic Clustering (Spargo et al., 2023)",
                  titleWidth = 600),
  ## dashboardSidebar  ####
  dashboardSidebar(
    ## siderbar menu ####
    sidebarMenu(
      menuItem(
        "Tutorial & Example Datasets",
        tabName = "exampleDataPage",
        icon = icon("fas fa-file")
      ),
      menuItem(
        "Phenotypic Clustering",
        tabName = "clusteringPage",
        icon = icon("fas fa-users")
      ),
      menuItem(
        "Transcriptional Age Analysis",
        tabName = "transcriptionalAgePage",
        icon = icon("fas fa-hourglass")
      ),
      menuItem(
        "Phenotypic Comparison",
        tabName = "phenotypicComparisonPage",
        icon = icon("fas fa-balance-scale")
      ),
      menuItem(
        "Full Results Table",
        tabName = "resultsPage",
        icon = icon("fas fa-table")
      ),
      menuItem("README", tabName = "readme", icon = icon("fas fa-info"))
    )
  ),
  ## dashboardBody ####
  dashboardBody(tabItems(
    tabItem(
      tabName = "clusteringPage",
      add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        br(),
        # Sidebar panel for inputs #####################################################
        sidebarPanel(
          # File upload for gene expression data
          fileInput(inputId = "geneExpressionFile",
                    label = "Upload Raw Gene Counts File with Ensembl IDs"),
          br(),
          actionButton("clusteringButton", "Perform Clustering")
        ),
        sidebarPanel(
          fileInput(inputId = "phenotypicFile",
                    label = "Upload Phenotypic File")
        ),

        sidebarPanel(
          selectInput(
            "phenotypicColumn",
            "Select a Phenotypic Column",
            NULL,
            selected = NULL,
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          )
        )
      ),
      box(
        width = "105%",
        height = "105%",
        status = "primary",
        solidHeader = TRUE,
        title = "Linear Discriminant Analysis Plot",
        #withSpinner(
        plotlyOutput(
          "ldaPlot",
          width = "100%",
          height = "600px",
          inline = TRUE
        )
        #)
      )
    ),
    tabItem(
      tabName = "transcriptionalAgePage",
      add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        br(),
        # Sidebar panel for inputs #####################################################
        sidebarPanel(
          # File upload for gene expression data
          selectInput(
            "tissue_type",
            "Select the Tissue Type",
            c("brain", "blood", "muscle", "nerve"),
            selected = "brain",
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          ),
          br(),
          actionButton("transcriptionalAgeButton", "Perform Analysis")

        ),
        sidebarPanel(
          selectInput(
            "population",
            "Select the Population Type",
            c("All Races", "North European/European"),
            selected = "All Races",
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          )
        ),
        sidebarPanel(
          selectInput(
            "age_column",
            "Select the Age Column",
            NULL,
            selected = NULL,
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          )
        )
      ),
      box(
        width = "105%",
        height = "105%",
        status = "primary",
        solidHeader = TRUE,
        title = "Transcriptional Analysis Plot",
        #withSpinner(
        plotlyOutput(
          "taPlot",
          width = "100%",
          height = "600px",
          inline = TRUE
        )
        #)
      )
    ),
    tabItem(
      tabName = "phenotypicComparisonPage",
      add_busy_spinner(spin = "fading-circle"),
      fluidRow(
        br(),
        sidebarPanel(
          selectInput(
            "phenotypicComparisonColumn",
            "Select a Continous Column to Compare Between Clusters",
            NULL,
            selected = NULL,
            multiple = FALSE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          ),
          br(),
          actionButton("phenotypicComparisonButton", "Perform Analysis")
        ),
        sidebarPanel(
          selectInput(
            "phenotypicComparisonCovariants",
            "Select Categorical Covariant(s)",
            NULL,
            selected = NULL,
            multiple = TRUE,
            selectize = TRUE,
            width = NULL,
            size = NULL
          )
        )
      ),
      box(
        width = "105%",
        height = "105%",
        status = "primary",
        solidHeader = TRUE,
        title = "Phenotypic Comparison Plot",
        #withSpinner(
        plotlyOutput(
          "pcPlot",
          width = "100%",
          height = "600px",
          inline = TRUE
        )
      ),
      br(),
      h2("Cluster Summary Statistic Table"),
      br(),
      fluidRow(column(
        DT::dataTableOutput("clusterSummaryStatisticTable"),
        width = 12
      )),
      br(),
      h2("Pairwise Test Results Table"),
      br(),
      fluidRow(column(
        DT::dataTableOutput("pairwiseTestResultsTable"),
        width = 12
      )),
      br(),
      fluidRow(
        downloadButton(
          "downloadClusterSummaryStatisticTable",
          "Download Summary Stats"
        ),
        downloadButton(
          "downloadPairwiseTestResultsTable",
          "Download Pairwise Results"
        )
      )
    ),
    tabItem(
      tabName = "resultsPage",
      br(),
      h2("Full Table of Clustering Results"),
      br(),
      fluidRow(
        column(
          DT::dataTableOutput("fulltable"),
          width = 12,
          style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
        )
      ),
      br(),
      downloadButton("downloadLdaResults", "Download")
    ),
    tabItem(tabName = "readme",
            includeMarkdown("./data/README.md")),
    tabItem(
      tabName = "exampleDataPage",
      br(),
      h2("Tutorial Video"),
      br(),
      HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/BM6fjephXiU" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'),
      h2("An Example Phenotypic File"),
      br(),
      downloadButton("downloadPFile", "Download Phenotypic File")
    )
  ))
)

# server ####
server <- function(input, output, session) {
  objects <-
    reactiveValues(
      geneExpressionDf = NULL,
      phenotypicDf = NULL,
      ldaResults = NULL
    )

  updatePhenotypicSelectInputs <- function(columnNames) {
    updateSelectInput(session,
                      "phenotypicColumn",
                      choices = columnNames,
                      selected = columnNames[1])

    updateSelectInput(
      session,
      "age_column",
      choices = c("NA", columnNames),
      selected = "NA"
    )

    updateSelectInput(session,
                      "phenotypicComparisonColumn",
                      choices = c(columnNames))

    updateSelectInput(
      session,
      "phenotypicComparisonCovariants",
      choices = c("NA", columnNames),
      selected = "NA"
    )
  }

  updatePhenotypicComparisonsSelectInputs <- function(columnNames) {
    updateSelectInput(session,
                      "phenotypicComparisonColumn",
                      choices = c(columnNames))

    updateSelectInput(
      session,
      "phenotypicComparisonCovariants",
      choices = c("NA", columnNames),
      selected = "NA"
    )
  }

  observeEvent(input$geneExpressionFile, {
    tryCatch({
      geneExpressionFile <- input$geneExpressionFile
      geneExpressionFileExt <-
        tools::file_ext(geneExpressionFile$datapath)

      req(geneExpressionFile)

      geneExpressionFlag <- is.null(geneExpressionFile)

      if (!geneExpressionFlag) {
        objects$geneExpressionDf <-
          fread(geneExpressionFile$datapath) %>%
          as.data.frame() %>% columnToRowNames(1)

        objects$geneExpressionDf <- objects$geneExpressionDf[, order(colnames(objects$geneExpressionDf))]
      }
    }, error = function(e) {
      showNotification(
        'An error occurred importing the phenotypic file. Please see the subsequent message for details.',
        type = "error",
        duration = 15
      )
      try(showNotification(toString(e), type = "error", duration = 30))
    })
  })

  observeEvent(input$phenotypicFile, {
    tryCatch({
      phenotypicFile = input$phenotypicFile
      phenotypicFileExt <- tools::file_ext(phenotypicFile$datapath)

      req(phenotypicFile)

      phenotypicFlag <- is.null(phenotypicFile)

      if (!phenotypicFlag) {
        #validate(need(phenotypicFileExt == "csv", "Please upload a csv file"))

        objects$phenotypicDf <- fread(phenotypicFile$datapath) %>%
          as.data.frame() %>% columnToRowNames(1)

        objects$phenotypicDf <- objects$phenotypicDf[order(rownames(objects$phenotypicDf)),]

        columnNames <- colnames(objects$phenotypicDf)

        updatePhenotypicSelectInputs(columnNames)
      }
    }, error = function(e) {
      showNotification(
        'An error occurred importing the phenotypic file. Please see the subsequent message for details.',
        type = "error",
        duration = 15
      )
      try(showNotification(toString(e), type = "error", duration = 30))
    })
  })

  observeEvent(input$clusteringButton, {
    tryCatch({
      if (is.null(objects$phenotypicDf)) {
        rowNamesList <- colnames(objects$geneExpressionDf)

        objects$phenotypicDf <-
          data.frame(ID = rowNamesList, row.names = rowNamesList)

        updatePhenotypicSelectInputs(c("ID"))
      }

      # create DEseq2 and normalise with size factors
      dds_genes <-
        DESeqDataSetFromMatrix(
          countData = objects$geneExpressionDf,
          colData = objects$phenotypicDf,
          design = ~ 1
        ) %>%
        estimateSizeFactors()

      # remove counts less than 10 in more than 5 individuals
      idx <-
        rowSums(counts(dds_genes, normalized = TRUE) >= 5) >= 10
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
      informative_gene_no <-
        length(rownames(dataset_informative_shared))

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

      colnames(predicted$posterior) <-
        c(
          "Probability of being in Cluster 1",
          "Probability of being in Cluster 2",
          "Probability of being in Cluster 3"
        )
      class(predicted$class)

      objects$ldaResults <-
        cbind(predicted$x, predicted$posterior) %>%
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

      output$fulltable <- DT::renderDataTable({
        as.data.frame(objects$ldaResults)
      })

      observeEvent(input$phenotypicColumn, {
        output$ldaPlot <- renderPlotly({
          req(input$phenotypicColumn)
          objects$phenotypicColumMachineFriendly <-
            reactive(str_replace_all(
              string = input$phenotypicColumn,
              pattern = " ",
              repl = ""
            ))

          plot_ly(data = objects$ldaResultsMachineFriendly) %>%
            add_trace(
              x = ~ LD1,
              y = ~ LD2,
              color = ~ ClusterAssignment,
              marker = list(size = 12),
              text = paste0(
                objects$phenotypicColumMachineFriendly(),
                " : ",
                objects$ldaResultsMachineFriendly[, objects$phenotypicColumMachineFriendly()]
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
            layout(
              legend = list(
                tracegroupgap = 50,
                y = 0.9,
                yanchor = "top"
              ),
              xaxis = list(title = 'Linear Discriminant 1'),
              yaxis = list(title = 'Linear Discriminant 2')
            )
        })
      })
    }, error = function(e) {
      showNotification(
        'An error occurred during the clustering analysis. Please see the subsequent message for details.',
        type = "error",
        duration = 15
      )
      try(showNotification(toString(e), type = "error", duration = 30))
    })
  })

  observeEvent(input$transcriptionalAgeButton, {
    tryCatch({
      if (input$population == "All Races") {
        population <- "all"
      } else {
        population <- "caucasian"
      }

      res <-
        predict_age(
          exprdata = objects$geneExpressionDf,
          tissue = input$tissue_type,
          exprtype = "counts",
          idtype = "ENSEMBL",
          stype = population
        )

      RNAAge <- res$RNAAge
      objects$ldaResultsMachineFriendly$TranscriptionalAge <- RNAAge
      objects$ldaResults[, "Transcriptional Age"] <- RNAAge

      if (!input$age_column == "NA") {
        RNAAge_accel <-
          res$RNAAge - objects$ldaResultsMachineFriendly[, input$age_column]
        objects$ldaResultsMachineFriendly$TranscriptionalAgeAcceleration <-
          RNAAge_accel
        objects$ldaResults[, "Transcriptional Age Acceleration"] <-
          RNAAge_accel

        output$taPlot <- renderPlotly({
          plot_ly(data = objects$ldaResultsMachineFriendly) %>%
            add_trace(
              y = ~ TranscriptionalAge,
              x = ~ objects$ldaResultsMachineFriendly[, str_replace_all(
                string = input$age_column,
                pattern = " ",
                repl = ""
              )],
              marker = list(size = 12),
              text = paste0(
                "Transcriptional Age Acceleration : ",
                objects$ldaResultsMachineFriendly$TranscriptionalAgeAcceleration
              ),
              color = ~ ClusterAssignment,
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
            layout(
              legend = list(
                tracegroupgap = 50,
                y = 0.9,
                yanchor = "top"
              ),
              xaxis = list(title = 'Age'),
              yaxis = list(title = 'Transcriptional Age')
            )
        })
      } else {
        output$taPlot <- renderPlotly({
          plot_ly(data = objects$ldaResultsMachineFriendly) %>%
            add_trace(
              y = ~ TranscriptionalAge,
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
              yaxis = list(title = 'Transcriptional Age')
            )
        })

      }

      ldaResultsColumnNames <- colnames(objects$ldaResults)

      updatePhenotypicComparisonsSelectInputs(ldaResultsColumnNames[!ldaResultsColumnNames %in%
                                                                      c(
                                                                        "LD1",
                                                                        "LD2",
                                                                        "Probability of being in Cluster 1",
                                                                        "Probability of being in Cluster 2",
                                                                        "Probability of being in Cluster 3",
                                                                        "Cluster Assignment"
                                                                      )])

    }, error = function(e) {
      showNotification(
        'An error occurred during the clustering. Please see the subsequent message for details.',
        type = "error",
        duration = 15
      )
      message(toString(e), type = "error")
      try(showNotification(toString(e), type = "error", duration = 30))
    })
  })

  observeEvent(input$phenotypicComparisonButton, {
    tryCatch({
      phenotypicComparisonColumn <-
        str_replace_all(
          string = input$phenotypicComparisonColumn,
          pattern = " ",
          repl = ""
        )

      phenotypicComparisonCovariants <-
        str_replace_all(
          string = input$phenotypicComparisonCovariants,
          pattern = " ",
          repl = ""
        )

      phenotypic_variable <-
        objects$ldaResultsMachineFriendly[, phenotypicComparisonColumn]
      # whatever phenotype is chosen

      for (column in phenotypicComparisonCovariants) {
        if (!column %in% c(
          "NA",
          "LD1",
          "LD2",
          "ProbabilityofbeinginCluster1",
          "ProbabilityofbeinginCluster2",
          "ProbabilityofbeinginCluster3",
          "TranscriptionalAge",
          "TranscriptionalAgeAcceleration"
        )) {
          typeOfColumn <- typeof(objects$ldaResultsMachineFriendly[, column])

          isColumnFactor <-
            !is.factor(objects$ldaResultsMachineFriendly[, column])

          if (typeOfColumn == "integer" & isColumnFactor) {
            if (max(objects$ldaResultsMachineFriendly[, column]) < 6) {
              if (max(objects$ldaResultsMachineFriendly[, column]) > -1) {
                objects$ldaResultsMachineFriendly[, column] <-
                  as.factor(objects$ldaResultsMachineFriendly[, column])
              }
            }
          } else if (typeOfColumn != "double" & isColumnFactor) {
            objects$ldaResultsMachineFriendly[, column] <-
              as.factor(objects$ldaResultsMachineFriendly[, column])
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
          mean = mean((
            !!sym(phenotypicComparisonColumn)
          ), na.rm = TRUE),
          sd = sd((
            !!sym(phenotypicComparisonColumn)
          ), na.rm = TRUE),
          median = median((
            !!sym(phenotypicComparisonColumn)
          ), na.rm = TRUE),
          IQR = IQR((
            !!sym(phenotypicComparisonColumn)
          ), na.rm = TRUE)
        ) %>% as.data.frame()

      colnames(objects$phenotypicComparisonSummaryResults)[1] <-
        "cluster"

      # normality of each variable assessed with shapiro wilk
      normality <-
        shapiro.test(objects$ldaResultsMachineFriendly[, phenotypicComparisonColumn])

      # variable is log-transformed before test if not normal
      if (normality$p.value < 0.05) {
        phenotypic_variable_transformed <- log(phenotypic_variable)
      } else {
        phenotypic_variable_transformed <- phenotypic_variable
      }

      # one-way ANOVA
      equation <-
        "phenotypic_variable_transformed ~ objects$ldaResultsMachineFriendly$ClusterAssignment"

      covariantsFlag <- length(phenotypicComparisonCovariants) > 1
      if (isFALSE(covariantsFlag)) {
        covariantsFlag <- phenotypicComparisonCovariants != "NA"
      }

      if (covariantsFlag) {
        for (covariant in phenotypicComparisonCovariants) {
          equation <- paste0(equation, " + ", covariant)
        }
      }

      res <- aov(as.formula(equation),
                 data = objects$ldaResultsMachineFriendly)

      anova_p_value <- summary(res)[[1]][["Pr(>F)"]][1]

      # tukey's test
      pairwise_res <- TukeyHSD(res)

      objects$phenotypicComparisonPairwiseResults <-
        as.data.frame(pairwise_res$`objects$ldaResultsMachineFriendly$ClusterAssignment`)

      tukey_p_values <-
        objects$phenotypicComparisonPairwiseResults[, 4]

      columnNamesOrder <-
        colnames(objects$phenotypicComparisonPairwiseResults)
      objects$phenotypicComparisonPairwiseResults$comparison <-
        rownames(objects$phenotypicComparisonPairwiseResults)
      objects$phenotypicComparisonPairwiseResults <-
        objects$phenotypicComparisonPairwiseResults[,
                                                    c("comparison", columnNamesOrder)]
      colnames(objects$phenotypicComparisonPairwiseResults) <-
        c(
          "comparison",
          "difference",
          "lower confidence interval limit",
          "upper confidence interval limit",
          "adjusted p-value"
        )

      tukey_first_group <- c("2", "3", "3")
      tukey_second_group <- c("1", "1", "2")

      # put the p values in a handy table for the plotting of the p values
      p_value_results <-
        data.frame(x1 = tukey_first_group, x2 = tukey_second_group, x3 = tukey_p_values)
      colnames(p_value_results) <- c("group1", "group2", "p.adj")


      output$pcPlot <- renderPlotly({
        req(input$phenotypicComparisonColumn)
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
      })

      output$clusterSummaryStatisticTable <- DT::renderDataTable({
        objects$phenotypicComparisonSummaryResults
      })

      output$pairwiseTestResultsTable <- DT::renderDataTable({
        objects$phenotypicComparisonPairwiseResults
      })

    }, error = function(e) {
      showNotification(
        'An error occurred during the phenotypic comparison. Please see the subsequent message for details.',
        type = "error",
        duration = 15
      )
      message(toString(e), type = "error")
      try(showNotification(toString(e), type = "error", duration = 30))
    })
  })


  observeEvent(objects$ldaResults, {
    try({
      output$fulltable <- DT::renderDataTable({
        datatable(objects$ldaResults, options = list(paging = FALSE))
      })
    })
  })

  output$downloadLdaResults <- try({
    downloadHandler(
      filename = function() {
        'ldaResults.csv'
      },
      content = function(con) {
        req(objects$ldaResults)
        fwrite(objects$ldaResults, con, row.names = TRUE)
      }
    )
  })

  output$downloadClusterSummaryStatisticTable <- try({
    downloadHandler(
      filename = function() {
        'summaryStatistics.csv'
      },
      content = function(con) {
        req(objects$ldaResults)
        fwrite(objects$phenotypicComparisonSummaryResults,
               con,
               row.names = TRUE)
      }
    )
  })

  output$downloadPairwiseTestResultsTable <- try({
    downloadHandler(
      filename = function() {
        'pairwiseRests.csv'
      },
      content = function(con) {
        req(objects$ldaResults)
        fwrite(objects$phenotypicComparisonPairwiseResults,
               con,
               row.names = TRUE)
      }
    )
  })


  output$downloadGEFile <- try({
    downloadHandler(
      filename = function() {
        'geneExpressionFile.csv'
      },
      content = function(con) {
        fwrite(fread("data/expressionData.csv"), con)
      }
    )
  })

  output$downloadPFile <- try({
    downloadHandler(
      filename = function() {
        'phenotypicFile.csv'
      },
      content = function(con) {
        fwrite(fread("data/phenotypicData.csv"), con)
      }
    )
  })
}

# run app ####
shinyApp(ui, server)
