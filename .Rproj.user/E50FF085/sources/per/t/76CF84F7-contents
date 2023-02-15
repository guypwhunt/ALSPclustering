try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18BCells"
columnNames <- append("fileName", gpr18BCellsColumnNames)

markersOrCells <- markersOrCellsClassification[3]

clusterNames <-  clusterColumns[3:4]

n.cores <- 1

# clusterName <- clusterNames[1]
# clusterNames <- c("clusters_phenograph")
# markersOrCell <- markersOrCells[1]
my.cluster <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

df <-
  fread(file=paste0(
    "./data/",
    directoryName,
    '/clusteringOutput/clusteringOutputs.csv'
  ))
df <- as.data.frame(df)

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    try(source("R/01_functions.R"))
    try(source("R/00_datasets.R"))

    loadlibraries()

    performAllDifferentialAbundanceTests(df,
                                         directoryName,
                                         columnNames,
                                         clusterName,
                                         markersOrCell)
  }
}
