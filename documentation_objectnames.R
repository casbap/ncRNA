DOCUMENTATION FOR MAJOR CHANGES IN PROJECT (thesis progression)

#Object renaming conventions
orgData <- unzipped data from DEE2
qcData <- QC data from DEE2
qcSummary <- quality data
totalSummary <- not quality data
qcPass <- pass data from QC
totalPass <- not pass data from QC
qcWarn <- warn data from QC
totalWarn <- not warn data from QC
qcFail <- fail data from QC
totalFail <- not fail data from QC
qualitySummary <- data frame combo of pass, warn, fail, summary
databasePass <- PASS samples from QC data
databasePassList <- list form of database_pass
orgPass <- filtered data using generated list from database_pass_list
orgPassWide <- converted matrix into wide form [save]
orgMetadata <- metadata from DEE2
orgMetadataPass <- filtered data to only include PASS samples
orgCountMetadata <- filtered genecount data and filtered metadata into a list

agg <- aggregate [save]
srx_agg(x, counts = "GeneCounts")
srx_agg(hsapiensCountMetadata)

dgeObj <- DGEList(agg)
DGElist object from to store read counts & info
DGEList(counts = matrix(0, 0, 0), lib.size = colSums(counts), norm.factors = rep(1,ncol(counts)), samples = NULL, group = NULL, genes = NULL, remove.zeros = FALSE)

**normAgg <- normalisation aggregate**  [save]
normAgg <- cpm(agg) 
**normAggScale** **<- normalised aggreated, scaled**  [save]

correlate <- correlation matrix (spearman)
correlate_melt <- molten data frame (stacks wide data columns into single column)

distClust <- distance matrix of normalised aggregate (scaled) data
hClust <- heirarhial clustering of clust
correlate_distClustVcoph <- correlation of distance matrix and cophenetic distance

dend <- dendogram
cutClust <- cutree of hClust
cutClustlength <- length of unique cutClusters (size)
cutClustvalues <- length clusters sorted  [save]
cutClustvalues_dynamic <- dynamic method (cuttreeHybrid) from library [save]

orgGO <- organism gene ontology data
orgGOmatrixwide <- matrix of GO data converted into wide form [save]
diffCountGO <- Gene IDs in the gene count data, but not in GO matrix
diffGOCount <- Gene IDs in the GO matrix, but not in gene count data
GOtable <- removal of diffGOCount [save]
GOtrain <- blinded index GOtable [save]