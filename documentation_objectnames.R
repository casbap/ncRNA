DOCUMENTATION FOR MAJOR CHANGES IN PROJECT (thesis progression)
#Object renaming conventions

  #DATAPREP

# orgData <- unzipped data from DEE2
# qcData <- QC data from DEE2
#   qcSummary <- quality data
#   totalSummary <- not quality data
#   qcPass <- pass data from QC
#   totalPass <- not pass data from QC
#   qcWarn <- warn data from QC
#   totalWarn <- not warn data from QC
#   qcFail <- fail data from QC
#   totalFail <- not fail data from QC
#   qualitySummary <- data frame combo of pass, warn, fail, summary
#   databasePass <- PASS samples from QC data
#    databasePassList <- list form of database_pass
#       orgPass <- filtered data using generated list from database_pass_list
#         orgPassWide <- converted matrix into wide form [save]
# orgMetadata <- metadata from DEE2
#   orgMetadataPass <- filtered data to only include PASS samples
#   orgCountMetadata <- filtered genecount data and filtered metadata into a list
# 
# agg <- aggregate [save]
#   srx_agg(x, counts = "GeneCounts")
#   srx_agg(hsapiensCountMetadata)
# 
# dgeObj <- DGEList(agg)
#   DGElist object from to store read counts & info
#   DGEList(counts = matrix(0, 0, 0), lib.size = colSums(counts), norm.factors = rep(1,ncol(counts)), samples = NULL, group = NULL, genes = NULL, remove.zeros = FALSE)
# 
# normAgg <- normalisation aggregate  [save]
#   normAgg <- cpm(agg) 
# normAggScale <- normalised aggreated, scaled  [save]
# 
# correlate <- correlation matrix (spearman)
#   correlate_melt <- molten data frame (stacks wide data columns into single column)
# 
# distClust <- distance matrix of normalised aggregate (scaled) data
# hClust <- heirarhial clustering of clust
# correlate_distClustVcoph <- correlation of distance matrix and cophenetic distance
# 
# dend <- dendogram
# cutClust <- cutree of hClust
#   cutClustlength <- length of unique cutClusters (size)
#   cutClustvalues <- length clusters sorted  [save]
#   cutClustvalues_dynamic <- dynamic method (cuttreeHybrid) from library [save]


  #GOPREP

# ensembl <- use Ensembl genomic data & biotypes

# orgGO <- organism gene ontology data
# orGOlong <- long-form genome ids
#   orgGOmatrixwide <- matrix of GO data converted into wide form [save] --> converts to dataframe
# diffCountGO <- Gene IDs in the gene count data, but not in GO matrix
# diffGOCount <- Gene IDs in the GO matrix, but not in gene count data
# GO_table <- removal of diffGOCount [save]
# GO_train <- blinded index GOtable [save]
# GO_test <- blinded data from GO_table
# normAggScale_test <- rows in normAggScale also in GO_test (from ensembl data)


  #IMPUTATIONBLINDED

# cl_list <- cluster total and minimum cluster size used in dynamic tree cut
# thresh_list <- thresholds 0.1-0.9
# sample_kfold <- optimise_impute(cl_list, thresh_list, cutClustvalues_dynamic, normAggScale, GO_train)
# kfold_dynamic_all [save]
# 
# data_list <- lapply(sample_kfold, function(x) ...
#                     kfold_dynamic_sum <- summary of mean prediction scores for each parameter [save]


  #CLUSTERS

# logcounts <- agg.rds from DataPrep
# clusterPass <- optimised cluster size [save]
# clusterPasslength <- length of clusters and number of clusters


  #CLUSTERSWITHGO

# orgMetadata <- metadata from DEE2
#   databasePassList <- column names of orgMetadataPass
#   orgMetadataPass <- filtered data to only include PASS samples
# srxAgg <- sequencing run aggregate function
# orgCountMetadata <- filtered genecount data and filtered metadata into a list
# agg <- srxAgg and orgCountMetadata dataframes in one function
#   aggFilterCols <- samples > 1 million reads
#   aggFilterRows <- genes > 0 counts
# 
# cpm <- cpm(agg)    aggregate normalisation of library size bias
# keepthresh <- genes that have at least two TRUEs in each rows of the threshold function
# keepCount <- highly expressed genes
# 
# dgeObj <- convert counts to DGEList object
# logcounts <- log2 counts per million from dgeObj
# 
# distClust <- distance matrix of normalised aggregate (scaled) data
# hClust <- heirarhial clustering of clust
# cluster <- optimising cluster size
#   clusterLength <- length of unique cluster
# 
# orgGO <- organism gene ontology data
#   orgGOmatrix
#   orgGOmatrixwide <- matrix of GO data [save]
# 
# fracBlind <- fraction of number of genes to be blinded
# numGenes <- not aggregate * fracBlind
# blind <- total genes to be blinded
# dictBlind <- dataframe of blind & rand [save]
# aggRand <- agg
# indBlind <- in dict blind
# 
# corrClust <- genes per cluster and correlation values
# corrAllClusters <- corrClust for all clusters (normalised gene counts and matrix of genes with the cluster number)
# 
# GOperClust <- list of dataframe containing GO terms belonging to a cluster
#   GOClust <- initial form of GOperClust
#   GOallClusters <- GOperClust for all clusters (matrix of fenes belonging to a GO term and matrix of genes with a cluster number)
# 
# weightcorrClust <- weight of genes to be correlated 
#   cluster1GO <- weightcorrClust function output
#   weightcorr2Clust <- sum of corr values used as denominator in normalized_weighted_go function
#     cluster1GO2 <- weightcorr2Clust function output
# 
# cluster <- data(cluster)      from chunk 8
# cluster1
#   cluster1_list
#   cluster1GO_nzero <- non-zero columns in GO dictionary
#   corrCluster1 <- correlation in cluster 1 (same as corrClust)
#   cluster1GO <- GO terms for cluster 1 (same as GOperClust)
#   corrblind1 <- weighted correlation for cluster 1
#   geneblind1 <- blinded gene form for cluster 1 
#   corrblind1df <- dataframe of corrblind1 with subset of row containing geneblind1
#   weightGOClust1 <- merged GO terms dataframe (cluster 1) and correlation values for geneblind1
#     normweightGOClust1 <- normalised weightGOClust1
#   cluster1wGO <- GO terms for cluster 1 but weighted
# 
# inpGO <- threshold from cluster1wGO for normweightGOClust1 > 0.1
#   cluster1wGOdf <- cluster1wGO as dataframe
# Clust1thresh <- threshold (0.005) of cluster 1
#   inpClust1df <- dataframe 1
#   inpClust1df2 <- dataframe 2
# 
# cluster1inp <- filtered GO matrix to only include genes of interest and GO terms from inputed dataframe
#   cluster1inp1 <- ordered by row and column names
#   inpClust1trim <- filtered inpClust1df to include cluster1inp
# Clust1inpVSinp <- input dataframe (cluster1inp) VS resultant dataframe (inpClust1trim)
# Clust1neg <- inpClust1trim - cluster1inp     negative values for threshold > 0.005
# 
# cluster1GO2df <- dataframe using weightcorr2Clust
# Clust1thresh <- 0.01
# Clust1inpVSinp2 <-  input dataframe (cluster1inp) VS resultant dataframe (inpClust1trim)
# 
# clust3GO <- GO terms for cluster 3
#   clust3blind <- blind genes (list form) in cluster
#   cluster3wGO2df <-  cluster3wGO as dataframe
# Clust3thresh <- 0.01
#   impClust3df <- cluster3wGO2df to Clust3thresh as dataframe
#   cluster3inp <- filtered GO matrix to only include genes of interest and GO terms from inputed dataframe 
# 
# orgCPGO <- enriched data
# orgCPGOdf <- orgCPGO as dataframe
# Clust3inpVSinp3 <- input dataframe (cluster3inp) VS resultant dataframe (inpClust3trim)
# Clust3neg <- inpClust3trim - cluster3inp     negative values for threshold > 0.005


  #IMPUTATION

# cluster1 (taken from ClustersWithGO)
#   cluster1_list <- list of cluster1
#   corrClust1 <- corrClust for cluster 1
#   cluster1GO <- GOperClust for cluster 1
#     cluster1GO_list <- genes of interest for cluster 1 (row names)
#     diffClust1GO <- GeneIDs in cluster not present in GO list
#   cluster1wGO <- weighted values for blinded genes in cluster 1
#     cluster1wGOdf <- dataframe of cluster1wGO
#       noGOs <- no GO term values in cluster1wGOdf
# 
# Clust1thresh <- 0.02
# impClust1df <- threshold for cluster1wGOdf
# 
# newRownames <- noGOs * 0
# input_mat <- input matrix       Add GeneIDs with no GO terms to the cluster 1 GO terms and set them to zero
# Clustnegneg <- impClust1df - input_mat
# 
# weightGOallClusters <- nested gene clusters contianing data for cluster analysis [save]
# GOtermsperClust <- list of GO terms per cluster from original database before blinding [save]
# 
# clusterVal <- min-max cut value of clusters [save]
# nullOnt <- number of gene IDs with no GO IDs 
