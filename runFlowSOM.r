# Command line arguments:
# 1. clean_data.fcs
# 2. number of meta clusters
# 3. flag to include method name as column
# 4. output directory
# 5. output file name for cell/cluster assignment
# 6. output file name for cluster mean feature values

# install packages
if(!require('flowCore')) {install.packages('flowCore')}
if(!require('FlowSOM')) {install.packages('FlowSOM')}

# load libraries
library("flowCore")
library("FlowSOM")


# get arguments
args <- commandArgs(trailingOnly=TRUE)

# read in data
data <- read.FCS(args[1])

# get the number of columns
num_cols <- length(colnames(data))

# check if logicle transformation is necessary
maxs <- vector() # initialize vec
for(i in 2:num_cols) { # loop through column indices (excluding the first one which is cell ID)
    maxs <- append(maxs,max(exprs(data)[,i])) # add max of column to vec
}
max = max(maxs) # get the max of all the maxs

# if the highest expression value is greater than 1000, logicle transform the data
if (max > 1000) {
    cols <- colnames(data)
    cols <- cols[cols != 'CellID'] # do not logicle transform cell IDs
    logicleTrans <- estimateLogicle(data, cols) # automatically estimate the logicle transformation based on the data
    data <- transform(data, logicleTrans) # apply logicle transformation
}

# run FlowSOM, cluster using all columns besides first (assuming it is the cell ID column)
fSOM <- FlowSOM(data, colsToUse=c(2:num_cols), nClus=as.integer(args[2]), compensate=FALSE, spillover=NULL)

# get cluster assignments
Cluster <- GetMetaclusters(fSOM)

# get raw input data and add clusters
data_raw <- cbind(Cluster, exprs(data))

# make cells.csv
cells <- data_raw[,c('CellID','Cluster')]
if (as.logical(args[3])) { # inlcude method column
    Method <- rep(c('FlowSOM'),nrow(cells))
    cells <- cbind(cells, Method)
}
write.table(cells,file=paste(args[4], args[5], sep='/'), row.names=FALSE, quote=FALSE, sep=',') # write data to csv

# make clusters.csv
# the averages are not log transformed, but the original values
clusterData <- aggregate(subset(data_raw, select=-c(CellID)), list(data_raw[,'Cluster']), mean) # group feature/expression data by cluster and find mean expression for each cluster, remove CellID column
clusterData <- subset(clusterData, select=-c(Group.1)) # remove group number column because is identical to community assignation number
if (as.logical(args[3])) { # inlcude method column
    Method <- rep(c('FlowSOM'),nrow(clusterData))
    clusterData <- cbind(clusterData, Method)
}
write.table(clusterData,file=paste(args[4], args[6], sep='/'),row.names=FALSE,quote=FALSE,sep=',') # write data to csv