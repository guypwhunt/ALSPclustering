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
library(caret)
library(xgboost)

# Load Machine Learning Models
firstVisitMLModel <- readRDS("data/data_pheno/XGBoost_AUC_maybeFinal/results/XGBfinalTune.Rds")
fullMLModel <- readRDS("data/data_LCAvars/XGBoost_AUC_final/results/XGBfinalTune.Rds")
nonCensoredMLModel <- readRDS("data/data_LCAvars_noncens/XGBoost_AUC_final/results/XGBfinalTune.Rds")

delay_panel<- read.delim("data/diagDelay_panel.tsv")

columnToRowNames <- function(df, index) {
  rownames(df) <- df[, index]
  df <- df[,-index]
  return(df)
}

scaleToPanel <- function(x,ID,panel){
  rowID<- panel$ID==ID        #define the panel row to scale against
  x <- x-panel$mean[rowID]    #center on mean
  x <- x/panel$sd[rowID]      #Scale to reference variance

  x
}

# Set Seed
set.seed(1)

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
        "Phenotypic Comparison",
        tabName = "phenotypicComparisonPage",
        icon = icon("fas fa-balance-scale")
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
          fileInput(inputId = "phenotypicFile",
                    label = "Upload Phenotypic File"),
          br(),
          br(),
          br(),
          br(),
          actionButton("clusteringButton", "Perform Clustering")
        ),
        sidebarPanel(
          radioButtons("mlModelChoice",
                       "Select a Machine Learning that is based on:",
                       choices = c("Data Available at First Vist",
                                   "Uncensored Patient's Data",
                                   "Censored and Uncensored Patient's Data"),
                       selected = "Data Available at First Vist"),
          br(),
          selectInput(
            "countryOfOrigin",
            "Select the Country of Diganosis for all the Patients",
            c("BE", "CH", "ES", "FR", "GB", "IE", "IL", "IT", "NL", "PT",
              "SE", "TR", "US", "Other"),
            selected = "Other")
          # selectInput(
          #   "phenotypicColumn",
          #   "Select a Phenotypic Column",
          #   NULL,
          #   selected = NULL,
          #   multiple = FALSE,
          #   selectize = TRUE,
          #   width = NULL,
          #   size = NULL
          # )
        )),
      br(),
      h2("Full Table of Clustering Results"),

      fluidRow(
        column(
          DT::dataTableOutput("fulltable"),
          width = 12,
          style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
        )
      ),
      br(),
      downloadButton("downloadClusteringResults", "Download"),
      br(),
      br()#,
      # box(
      #   width = "105%",
      #   height = "105%",
      #   status = "primary",
      #   solidHeader = TRUE,
      #   title = "Linear Discriminant Analysis Plot",
      #   plotlyOutput(
      #     "ldaPlot",
      #     width = "100%",
      #     height = "600px",
      #     inline = TRUE
      #   )
      # )
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
      phenotypicDf = NULL
    )

  observeEvent(input$phenotypicFile, {
    tryCatch({
      phenotypicFile = input$phenotypicFile

      req(phenotypicFile)

      phenotypicFlag <- is.null(phenotypicFile)

      if (!phenotypicFlag) {
        objects$phenotypicDf <- fread(phenotypicFile$datapath) %>%
          as.data.frame() %>% columnToRowNames(1)

        objects$phenotypicDf <- objects$phenotypicDf[order(rownames(objects$phenotypicDf)),]

        columnNames <- colnames(objects$phenotypicDf)

        updateSelectInput(session,
                          "phenotypicComparisonColumn",
                          choices = columnNames)

        updateSelectInput(
          session,
          "phenotypicComparisonCovariants",
          choices = c("NA", columnNames),
          selected = "NA"
        )

        # updateSelectInput(session,
        #                   "phenotypicColumn",
        #                   choices = columnNames,
        #                   selected = columnNames[1])
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

      if(input$mlModelChoice == "Data Available at First Vist"){
        mlModel <- firstVisitMLModel
      } else if (input$mlModelChoice == "Uncensored Patient's Data") {
        mlModel <- nonCensoredMLModel
      } else {
        mlModel <- fullMLModel
      }

      if(!input$countryOfOrigin == "Other") {

        objects$phenotypicDf <- objects$phenotypicDf %>%
          mutate(
            CTYDEL = scaleToPanel(CTYDEL,ID = input$countryOfOrigin,
                                  panel=delay_panel)
          )

      } else {
        showNotification("Diagnositic delay has not been standised.
                         Please ensure you have manually standised it.",
                         type = "warning")
      }


      objects$phenotypicDf$Cluster <- predict(mlModel, newdata = objects$phenotypicDf) %>%
        as.character() %>% str_replace_all("X", "")

      subselectedDataset <- objects$phenotypicDf[, c(
        "CTYDEL", "Sex_at_birth", "Site_of_Onset", "Age_at_onset_years",
        "Time_to_death_or_last_followup_years", 	"Phenotype_1",
        "Phenotype_2", "Phenotype_3", "Cluster")]

      # model <- lda(Cluster ~ ., data = subselectedDataset)
      #
      # predicted <- predict(model, subselectedDataset)
      #
      # objects$ldaResults <-
      #   cbind(predicted$x, predicted$posterior) %>%
      #   as.data.frame()
      #
      # objects$ldaResultsMachineFriendly <- objects$ldaResults
      #
      # colnames(objects$ldaResultsMachineFriendly) <-
      #   str_replace_all(
      #     string = colnames(objects$ldaResultsMachineFriendly),
      #     pattern = " ",
      #     repl = ""
      #   )
      #
      # observeEvent(input$phenotypicColumn, {
      #   output$ldaPlot <- renderPlotly({
      #     req(input$phenotypicColumn)
      #     objects$phenotypicColumMachineFriendly <-
      #       reactive(str_replace_all(
      #         string = input$phenotypicColumn,
      #         pattern = " ",
      #         repl = ""
      #       ))
      #
      #     plot_ly(data = objects$ldaResultsMachineFriendly) %>%
      #       add_trace(
      #         x = ~ LD1,
      #         y = ~ LD2,
      #         color = ~ ClusterAssignment,
      #         marker = list(size = 12),
      #         text = paste0(
      #           objects$phenotypicColumMachineFriendly(),
      #           " : ",
      #           objects$ldaResultsMachineFriendly[, objects$phenotypicColumMachineFriendly()]
      #         ),
      #         type = 'scatter',
      #         mode = 'markers',
      #         legendgroup = "Cluster",
      #         showlegend = T
      #       ) %>%
      #       add_annotations(
      #         text = "Clusters",
      #         xref = "paper",
      #         yref = "paper",
      #         x = 1.02,
      #         xanchor = "left",
      #         y = 0.9,
      #         yanchor = "bottom",
      #         # Same y as legend below
      #         legendtitle = TRUE,
      #         showarrow = FALSE
      #       ) %>%
      #       #Increase distance between groups in Legend
      #       layout(
      #         legend = list(
      #           tracegroupgap = 50,
      #           y = 0.9,
      #           yanchor = "top"
      #         ),
      #         xaxis = list(title = 'Linear Discriminant 1'),
      #         yaxis = list(title = 'Linear Discriminant 2')
      #       )
      #   })
      # })

    }, error = function(e) {
      showNotification(
        'An error occurred during the clustering analysis. Please see the subsequent message for details.',
        type = "error",
        duration = 15
      )
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
        objects$phenotypicDf[, phenotypicComparisonColumn]
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
          typeOfColumn <- typeof(objects$phenotypicDf[, column])

          isColumnFactor <-
            !is.factor(objects$phenotypicDf[, column])

          if (typeOfColumn == "integer" & isColumnFactor) {
            if (max(objects$phenotypicDf[, column]) < 6) {
              if (max(objects$phenotypicDf[, column]) > -1) {
                objects$phenotypicDf[, column] <-
                  as.factor(objects$phenotypicDf[, column])
              }
            }
          } else if (typeOfColumn != "double" & isColumnFactor) {
            objects$phenotypicDf[, column] <-
              as.factor(objects$phenotypicDf[, column])
          }
        }
      }

      # descriptive statistics for the variable
      objects$phenotypicComparisonSummaryResults <-
        group_by(
          objects$phenotypicDf,
          objects$phenotypicDf$Cluster
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
        shapiro.test(objects$phenotypicDf[, phenotypicComparisonColumn])

      # variable is log-transformed before test if not normal
      if (normality$p.value < 0.05) {
        phenotypic_variable_transformed <- log(phenotypic_variable)
      } else {
        phenotypic_variable_transformed <- phenotypic_variable
      }

      # one-way ANOVA
      equation <-
        "phenotypic_variable_transformed ~ objects$phenotypicDf$Cluster"

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
                 data = objects$phenotypicDf)

      anova_p_value <- summary(res)[[1]][["Pr(>F)"]][1]

      # tukey's test
      pairwise_res <- TukeyHSD(res)

      objects$phenotypicComparisonPairwiseResults <-
        as.data.frame(pairwise_res$`objects$phenotypicDf$Cluster`)

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

      output$pcPlot <- renderPlotly({
        req(input$phenotypicComparisonColumn)
        plot_ly(data = objects$phenotypicDf) %>%
          add_trace(
            y = objects$phenotypicDf[, phenotypicComparisonColumn],
            x = ~ Cluster,
            color = ~ Cluster,
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



  observeEvent(objects$phenotypicDf, {
    try({
      output$fulltable <- DT::renderDataTable({
        datatable(objects$phenotypicDf, options = list(paging = FALSE))
      })
    })
  })

  output$downloadClusteringResults <- try({
    downloadHandler(
      filename = function() {
        'clusteringResults.csv'
      },
      content = function(con) {
        req(objects$phenotypicDf)
        fwrite(objects$phenotypicDf, con, row.names = TRUE)
      }
    )
  })

  output$downloadClusterSummaryStatisticTable <- try({
    downloadHandler(
      filename = function() {
        'summaryStatistics.csv'
      },
      content = function(con) {
        req(objects$phenotypicDf)
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
        req(objects$phenotypicDf)
        fwrite(objects$phenotypicComparisonPairwiseResults,
               con,
               row.names = TRUE)
      }
    )
  })

  output$downloadPFile <- try({
    downloadHandler(
      filename = function() {
        'phenotypicFile.csv'
      },
      content = function(con) {
        fwrite(fread("data/exampleData.csv"), con)
      }
    )
  })
}

# run app ####
shinyApp(ui, server)
