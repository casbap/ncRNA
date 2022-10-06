### ------------- Data Preparation

# Sequencing Run Aggregate function (srx_agg)
# corresponds srr forms to its experiment, srx, for filtered gene counts and metadata 
# 
# Inputs:
#      Genomic Metadata (From DEE2)
# 
# used in DataPrep (chunk 5)

srx_agg <- function(x,counts="GeneCounts") {
  IDX=which(names(x) %in% "GeneCounts")
  mds<-x$MetadataSummary
  n=nrow(x[[IDX]])
  SRX_dat <- vapply(X=unique(mds$SRX_accession) ,function(srx) {
    srrs<-rownames(mds)[which(mds$SRX_accession %in% srx)]
    if (length(srrs)>1) {
      rowSums(x[[IDX]][,srrs])
    } else {
      x[[IDX]][,srrs]
    }
  } , numeric(n))
  rownames(SRX_dat) <- rownames(x[[IDX]])
  colnames(SRX_dat) <- unique(mds$SRX_accession)
  SRX_dat
}


# Cluster lengths (cl_cut_dynamic)
# Uses cutreeHybrid() function from the 'dynamicTreeCut' library to iteratively cut a 
# dendrogram based on a minimum cluster size.
#
# Outputs:
#      list (cutClustvalues_dynamic) grouped by total clusters which gives:
#          (1) a table of GeneIDs assigned to a cluster number,
#          (2) the cutree denominator to achieve that total cluster, and
#          (3) the tally of the total Gene IDs belonging to a cluster number.
# Inputs:
#      hr (hClust) = hierarchical clustering object
#      cl (distClust) = distance matrix used to build 'hClust'
#      min and max = minimum and maximum value of the sequence of values used
#                    in the 'minClusterSize' attribute of the cutreeHybrid function 
#      interval = intervals for the sequence of values between min and max
# 
# Used in DataPrep (Chunk 8)

cl_cut_dynamic <- function(hr, cl, min, max, interval){
  
  cutClust_cuts <- list()
  
  cut <- seq(from = min, to = max, by = interval)
  
  for (i in cut){
    
    dyn_tree <- cutreeHybrid(
      # Input data: basic tree cutting
      dendro=hr, distM=as.matrix(cl),
      # Branch cut criteria and options
      minClusterSize = i, deepSplit = 4,
      # PAM stage options
      pamStage = TRUE, pamRespectsDendro = FALSE,
      respectSmallClusters = TRUE,
      # Various options
      verbose = 1, indent = 0)
    
    dyn_clust <- as.data.frame(dyn_tree$labels)
    colnames(dyn_clust) <- "ClusterNumber"
    dyn_clust$ensembl_gene_id <- hr$labels
    
    # Tally the total GeneIDs assigned per cluster
    cutClust_length <- length(unique(dyn_tree$labels))
    tally <- list()
    for (j in 1:cutClust_length) {
      tot <- length(which(dyn_clust$ClusterNumber == j))
      tally[[j]] <- tot
    }
    
    tallyDF <- t(as.data.frame(tally))
    rownames(tallyDF) <- c(1:cutClust_length)
    
    
    cutClust_cuts[[paste0(cutClust_length)]][[paste0("GeneID_assignments")]] <- dyn_clust
    cutClust_cuts[[paste0(cutClust_length)]][[paste0("min_cl_size")]] <- i
    cutClust_cuts[[paste0(cutClust_length)]][["Tally"]] <- tallyDF
  }
  return(cutClust_cuts)
}


### ------------- Cluster Analysis


# Correlation of gene counts within a single cluster (corr_per_clust)
#
# Outputs:
#      nested list of the genes grouped per cluster and their corresponding correlation values (corrClust)
# Inputs:
#      x = normAggLog (normalized log genes from RNA seq)
#      y = cutClust (matrix of genes with associated cluster number)
#      clust_total = total number of clusters (derived from cutClustvalues_dynamic)
#
# Used in Imputation (chunk 2)

corr_per_clust <- function(x, y, clust_total){
  corr_cl <- list()
  for (i in 1:clust_total){
    clusterX <- y[y$ClusterNumber == i,]
    cluster_list <- as.list(clusterX$GeneID)
    cluster <- x[rownames(x) %in% cluster_list,]
    corr_result <- cor(t(cluster))
    corr_cl[[paste0("Cluster", i)]] <- corr_result
  }
  return(corr_cl)
}


# Association between GO terms and genes within a cluster (GO_per_cl)
#
# Output:
#      list of data frames containing all GO terms associated with the genes 
#      belonging to a particular cluster (GOperClust)
# Inputs:
#      x = GO_table (matrix of genes belonging to a GO term)
#      y = cutClust (matrix of genes w/ cluster number)
#      clust_total = total number of clusters
#
# Used in Imputation (chunk 2)

GO_per_cl <- function(x,y,clust_total){
  GO_cl <- list()
  for (i in 1:clust_total){
    cluster_list <- names(y[y == i])
    cluster_GOterms <- x[x$ensembl_gene_id %in% cluster_list,]
    rownames(cluster_GOterms)<- cluster_GOterms[,1] 
    cluster_GOterms[,1] <- c()
    cluster_GOterms <- cluster_GOterms[,which(colSums(cluster_GOterms) > 0)]
    GO_cl[[paste0("Cluster", i)]] <- cluster_GOterms
  }
  return(GO_cl)
}


# Association between GO terms and genes within a cluster prior to blinding (GO_per_cutClustlength)
# GO_per_cl function modified to accommodate different input data structure of the blinded genes.
#
# Outputs:
#       nested list of GO terms (column names) grouped per cluster from the 
#       original input data frame *before* blinding (GOperClustall)
# Inputs:
#       x = list object containing the input data frame before blinding (GO_table)
#       clust_total = total number of clusters (derived from cutClustvalues_dynamic)
#
# Used in Imputation (chunk 2) to make final list of GO terms per cluster

GO_per_cl_list <- function(x,clust_total){
  GO_names <- list()
  
  for (i in 1:clust_total){
    input <- x[[i]][[1]]
    GO_list <- colnames(input)
    
    GO_names[[paste0("Cluster", i)]] <- GO_list
  }
  return(GO_names)
}


# Association between GO terms and genes within a cluster after blinding (GO_per_cl_blinded)
# GO_per_cl function modified to accommodate the different input data structure of the blinded genes.
#
# Outputs:
#        list of data frames containing all GO terms (from blinded data) 
#        associated with the genes belonging to a cluster (GO_blindedCl)
# Inputs:
#        x = Annotation matrix after blinding (GO_test)
#        y = matrix of genes with cluster number (cutClust)
#        GO_per_cutClustlength = a nested list object containing the GO terms (column names)
#                         of the org genes (before blinding) in each cluster 
#        clust_total = total number of clusters (derived from cutClustvalues_dynamic)
#
# Used in to modify the 'impute' function and for cross validation (cross_val) process

GO_per_cl_blinded <- function(x,y,GO_per_cutClustlength,clust_total){
  GO_cl <- list()
  for (i in 1:length(GO_per_cutClustlength)){
    clusterX <- y[y$ClusterNumber == i,]
    cluster_list <- as.list(clusterX$ensembl_gene_id)
    cluster_GOterms <- x[x$ensembl_gene_id %in% cluster_list,]
    rownames(cluster_GOterms)<- cluster_GOterms[,1] 
    cluster_GOterms[,1] <- c()
    cluster_GOterms <- cluster_GOterms[, colnames(cluster_GOterms) %in% GO_per_cl_list[[i]] ]
    GO_cl[[paste0("Cluster", i)]] <- cluster_GOterms
  }
  return(GO_cl)
}


### ------------- Imputation functions

# Gene pairs correlation weights (edgelist)
#
# Outputs:
#     list that contains a dataframe with the pair of genes that makes up an edge (listed from and to) 
#     with a corresponding correlation value treated as the edge weight (edgelist)
#     This list is created using the igraph library.
# Inputs:
#     corrClust = a list of correlation matrices per cluster (from 'corr_per_clust' function)
#     clust_total = total number of clusters (derived from cutClustvalues_dynamic)
#
# Used in Imputation (chunk 2) 

edgelist <- function(corrClust, clust_total){
  e_list <- list()
  
  for(i in 1:clust_total){
    # Create igraph object
    cor_g <- graph_from_adjacency_matrix(corrClust[[i]], mode='directed', weighted = 'correlation', 
                                         diag = TRUE)
    # Extract edge list
    cor_edge_list <- get.data.frame(cor_g, 'edges')
    
    e_list[[paste0("Cluster", i)]] <- cor_edge_list
  }
  return(e_list)
}


# Gene counts within clusters (gene_list)
#
# Outputs:
#      list of genes for use in cluster analysis imputation stages
# Inputs:
#      corrClust ('output of corr_per_clust' function)
#      clust_total (length of cutClustvalues_dynamic object)
#
# Used in Imputation (chunk 2)

gene_list <- function(corrClust)
for (i in 1:clust_total){
  print(i)
  corr_cl <- corrClust[[i]]
corr_cl <- corr_cl[order(rownames(corr_cl)),]
corr_cl <- corr_cl[,order(colnames(corr_cl))]
gene_list <- rownames(corr_cl)

}

# Correlations values within a certain threshold connecting genes and GO terms (weightcorrClust)
# Weighted values for all blinded genes in a single cluster.
# 
# Outputs: (clusterXwGO), where X = cluster number
#      list of genes belonging to a cluster with a tally of how many correlation values passed the threshold
#      that connects the gene to a GO term
# Inputs:
#      gene_list = list of genes to be correlated
#      edgelist = an edge list (from igraph Library) of the correlation values in a cluster
#      GO_cl = a list of GO terms for a cluster (from GO_per_clust function)
#      thresh = threshold from 0 to 1
#
# Used in Imputation (chunk 2) 


weightcorrClust <- function(gene_list, edgelist, GO_cl, thresh){
  wcorr_result <- list()
  
  for (gene in gene_list){
    corr_value_df <- edgelist[edgelist$from %in% gene, c(2,3)]
    rownames(corr_value_df) <- corr_value_df$to
    weighted_go <- merge(x = corr_value_df, y = GO_cl, by = "row.names", all.x = TRUE)
    rownames(weighted_go) <- weighted_go[,1]
    weighted_go[,1:2] <- c()
  # turn NA to zero for GeneIDs with no GO terms
    weighted_go[is.na(weighted_go)] <- 0
    weighted_go <- weighted_go[,1]*weighted_go[,2:ncol(weighted_go)]
    wGO_thresh <- (as.matrix(weighted_go) > thresh)*1 
    imputed_go <- colSums(wGO_thresh)
    wcorr_result[[gene]] <- imputed_go
  }
  setNames(wcorr_result, paste0(gene))
  return(wcorr_result)
}


# Impute function
# 
# Outputs:
#      nested list of gene clusters containing the data generated for cluster analysis (weight)
# Inputs:
#      GOperClust = GO terms per cluster from GO_per_cl
#      GO_clall = GO terms per cluster from 'GO_per_cl_blinded' function
#      corrClust = correlation values grouped by cluster
#      clust_total = total number of clusters (derived from cutClustvalues_dynamic)
#      cor_edge_list = an edge list (from igraph Library) from the correlation 
#                      values per cluster
#      thresh = optimal threshold value from 0 to 1 ()
#
# Used in Imputation (chunk 2) 

impute <- function (GOperClust, GO_clall, corrClust, clust_total, cor_edge_list, thresh){
  wGO_list <- list()
  
  for (i in 1:clust_total){
    print(i)
    corr_cl <- corrClust[[i]]
    GO_cl <- GO_clall[[i]]
    GO_cl_orig <- GOperClust[[i]]
    
    #If all the genes in the cluster have no GOs and the matrix is empty
    if (is_empty(corr_cl) == TRUE || isTRUE(corr_cl == 0) || 
        is_empty(GO_cl_orig) == TRUE || isTRUE(GO_cl_orig == 0) ){
      wGO_list[[paste0("Cluster", i)]][[paste0("Comment")]] <- "No Gene Ontologies found"
      wGO_list[[paste0("Cluster", i)]][[paste0("Input")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Output")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Blinded_GO")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Diff")]] <- 0
      
      next
    }
    
    # should add gene that are not included in the GO data otherwise merge will yield weird results
    diff_go_genes <- setdiff(gene_list,rownames(GO_cl))
    if(length(diff_go_genes) > 0) {
      add <- data.frame(matrix(0,nrow=length(diff_go_genes),ncol=ncol(GO_cl)))
      rownames(add) <- diff_go_genes
      colnames(add) <- colnames(GO_cl)
      GO_cl <- rbind(GO_cl,add)
    }
    GO_cl <- GO_cl[order(rownames(GO_cl)),]
    GO_cl <- GO_cl[,order(colnames(GO_cl))]
    
    # Do the same process above for the original annotation matrix to use for model validation
    diff_go_genes2 <- setdiff(gene_list,rownames(GO_cl_orig))
    # Check and add non-annotated genes
    if(length(diff_go_genes2) > 0) {
      add <- data.frame(matrix(0,nrow=length(diff_go_genes2),ncol=ncol(GO_cl_orig)))
      rownames(add) <- diff_go_genes2
      colnames(add) <- colnames(GO_cl_orig)
      GO_cl_orig <- rbind(GO_cl_orig,add)
    }
    # Check and filter GO terms from blinded GO matrix (GO_cl)
    diff_GO_table <- setdiff(colnames(GO_cl_orig),colnames(GO_cl))
    if(length(diff_GO_table) > 0){
      GO_cl_orig <- GO_cl_orig[, colnames(GO_cl_orig) %in% colnames(GO_cl)]
    }
    
    GO_cl_orig <- GO_cl_orig[order(rownames(GO_cl_orig)),]
    GO_cl_orig <- GO_cl_orig[,order(colnames(GO_cl_orig))]
    
    edgelist <- cor_edge_list[[i]]
    wGO_cl <- weightcorrClust(gene_list, edgelist, GO_cl, thresh)
    wGO_df <- as.data.frame(do.call(rbind, wGO_cl))
    wGO_thresh <- (as.matrix(wGO_df) > 0)*1 
    wGO_thresh <- wGO_thresh[order(rownames(wGO_thresh)),]
    wGO_thresh <- wGO_thresh[,order(colnames(wGO_thresh))]
    
    cl_subtract <- wGO_thresh - GO_cl
    
    wGO_list[[paste0("Cluster", i)]][[paste0("Input")]] <- GO_cl_orig
    wGO_list[[paste0("Cluster", i)]][[paste0("Output")]] <- wGO_thresh
    wGO_list[[paste0("Cluster", i)]][[paste0("Blinded_GO")]] <- GO_cl
    wGO_list[[paste0("Cluster", i)]][[paste0("Diff")]] <- cl_subtract
    
  }
  return(wGO_list)
}


### ------------- Measures of performance 

# Initial performance measures
#
# Outputs:
#      nested list of vectors that measures the performance of the imputation using the 
#      blinded and original data frames grouped per cluster.
# Inputs:
#      imputed_df = data frame, the imputed binary matrix from `impute()`
#      clust_total = total number of clusters (derived from cutClustvalues_dynamic)


stats_cl <- function(imputed_df, clust_total){
  stats_list <- list()
  
  for (i in 1:clust_total){
    
    input <- imputed_df[[i]][["Input"]]
    blind <- imputed_df[[i]][["Output"]]
    
    diff <- input - blind 
    sum <- input + blind
    
    TP <- length(which(sum == 2)) #True positive
    TN <- length(which(sum == 0)) #True negative
    FP <- length(which(diff == -1)) #False positive
    FN <- length(which(diff == 1)) #False negative
    
    # False negative rate (FNR)
    FNR <- FN/(FN+TP)
    # Sensitivity/True Positive Rate (TPR)
    TPR <- 1 - FNR
    # False Positive Rate (FPR)
    FPR <- FP/(FP+TN)
    # Specificity/True Negative Rate (TNR)
    TNR <- 1 - FPR
    # Precision/Positive Predictive Value (PPV)
    PPV <- TP/(TP+FP)
    # Accuracy (ACC)
    ACC <- (TP+TN)/(TP+TN+FP+FN)
    # F1 score (is the harmonic mean of precision and sensitivity)
    F1 <- (2*TP)/((2*TP)+FP+FN)
    # Recall
    Recall <- TP/(TP+FN)
    
    stats_list[[paste0("Cluster", i)]][[paste0("TP")]] <- TP
    stats_list[[paste0("Cluster", i)]][[paste0("TN")]] <- TN
    stats_list[[paste0("Cluster", i)]][[paste0("FP")]] <- FP
    stats_list[[paste0("Cluster", i)]][[paste0("FN")]] <- FN
    
    stats_list[[paste0("Cluster", i)]][[paste0("Sensitivity(TPR)")]] <- TPR
    stats_list[[paste0("Cluster", i)]][[paste0("Specificity(TNR)")]] <- TNR
    stats_list[[paste0("Cluster", i)]][[paste0("FPR")]] <- FPR
    stats_list[[paste0("Cluster", i)]][[paste0("FNR")]] <- FNR
    stats_list[[paste0("Cluster", i)]][[paste0("Precision(PPV)")]] <- PPV
    stats_list[[paste0("Cluster", i)]][[paste0("Accuracy(ACC)")]] <- ACC
    stats_list[[paste0("Cluster", i)]][[paste0("F1_Score")]] <- F1
    stats_list[[paste0("Cluster", i)]][[paste0("Recall")]] <- Recall
  }
  return(stats_list)
}


# Dynamic performance measures 
#
# Outputs:
#     list of values measuring the performance of the whole imputation using the blinded 
#     and original data frame ()
# Inputs:
#     stats_cl = list, statistical measures of performance per cluster

stats_all <- function(stats_cl){
  
  stats_total <- list()
  
  # Transform nested list to data frame
  stats_df <- as.data.frame(t(as.data.frame(
    lapply(stats_cl, unlist))))
  
  TP <- sum(stats_df$TP) #True positive
  TN <- sum(stats_df$TN) #True negative
  FP <- sum(stats_df$FP) #False positive
  FN <- sum(stats_df$FN) #False negative
  
  # False negative rate (FNR)
  FNR <- FN/(FN+TP)
  # Sensitivity/True Positive Rate (TPR)
  TPR <- 1 - FNR
  # False Positive Rate (FPR)
  FPR <- FP/(FP+TN)
  # Specificity/True Negative Rate (TNR)
  TNR <- 1 - FPR
  # Precision/Positive Predictive Value (PPV)
  PPV <- TP/(TP+FP)
  # Negative Predictive Value (NPV)
  NPV <- TN/(TN+FN)
  # accuracy (ACC)
  ACC <- (TP+TN)/(TP+TN+FP+FN)
  # balanced accuracy (BA)
  BA <- (TPR+TNR)/2
  # F1 score (is the harmonic mean of precision and sensitivity)
  F1 <- (2*TP)/((2*TP)+FP+FN)
  Recall <- TP/(TP+FN)
  
  stats_total[["Stats_df"]] <- stats_df
  
  stats_total[["TP_all"]] <- TP
  stats_total[["TN_all"]] <- TN
  stats_total[["FP_all"]] <- FP
  stats_total[["FN_all"]] <- FN
  
  stats_total[["Sensitivity"]] <- TPR
  stats_total[["Specificity"]] <- TNR
  stats_total[["FPR"]] <- FPR
  stats_total[["FNR"]] <- FNR
  stats_total[["Precision"]] <- PPV
  stats_total[["Accuracy"]] <- ACC
  stats_total[["F1_Score"]] <- F1
  stats_total[["Recall"]] <- Recall
  
  return(stats_total)
}


# Statistical Analysis of cross validation (cross_val) 
# Outputs:
#       list of the measures of performance from the stats_all function applied to a k-fold cross validation process
# Inputs:
#       n = number of folds (for k_fold, use n = 10)
#       GO_table = original GO annotation matrix
#       cutClust = matrix of genes w/ cluster number
#       GOperClust = GO terms per cluster from GO_per_cl
#       GOperClustall = a nested list object containing the GO terms (column names) before blinding in each cluster
#       corrClust = correlation values grouped by cluster
#       clust_total = total number of clusters (derived from cutClustvalues_dynamic)
#       thresh = threshold value from 0 to 1 

cross_val <- function(n, GO_table, cutClust, GOperClust,
                      GOperClustall, corrClust,
                      clust_total, thresh){
  
  stats <- list()
  
  # partition the data
  nfolds <- n
  set.seed(42)
  folds <- sample(1:nfolds, nrow(GO_table), replace = TRUE)
  
  dict_folds <- data.frame(GO_table$ensembl_gene_id, folds)
  
  for (z in 1:n){
    
    which_fold <- z
    print(paste0("Fold ", z))
    ind_blinded <- which(dict_folds$folds == which_fold)
    GO_blinded <- GO_table
    GO_blinded[ind_blinded, 2:ncol(GO_blinded)] = 0
    
    # Get the GO terms by cluster using the blinded data and the GO list per cluster
    GO_blindedCl <- GO_per_cl_blinded(GO_blinded, cutClust, 
                                      GOperClustall, clust_total)
    
    # Get the edge list per cluster
    cor_edge_list <- edge_list(correlate=corrClust, clust_total)
    
    # Impute function
    wGO_blinded <- impute(GOperClust, GO_clall=GO_blindedCl,corrClust, clust_total, cor_edge_list, thresh)
    
    # stats
    stats_perCl <- stats_cl(wGO_blinded, clust_total)
    stats_final <- stats_all(stats_perCl)
    
    stats[[paste0("Fold", z)]] <- stats_final
  }
  stats[["Index_folds"]] <- dict_folds
  return(stats)
  
}

### ------------- Optimisation 

# Outputs:
#      data frame of the averaged value of a specific measure of performance from all cross validation folds
# Inputs:
#      kfold_all = a data frame containing all the measures of performance for each fold (from 'optimise_impute' function)
#      stat_type =  a number indicating the measure of performance to be averaged from 1 to 10 in the ff order: 
#                   Total Positive (TP), Total Negative (TN), False Positive (FP), False Negative (FN), 
#                   Sensitivity (TPR), Specificity (TNR), Precision (PPV), and F1 Score (F1)

mean_Csweep <- function(kfold_all, stat_type) {
  mean_df <- data.frame()
  
  mean_df <- do.call(rbind.data.frame, lapply(kfold_all[[1]], "[", stat_type))
  colnames(mean_df)[1] <- names(kfold_all[1])
  
  for (i in 2:length(kfold_all)){
    col_name <- names(kfold_all[i])
    mean_df[col_name] <- do.call(rbind.data.frame, lapply(kfold_all[[i]], "[", stat_type))
  }
  
  Mean <- colMeans(mean_df)
  mean_df <- rbind(mean_df, Mean)
  rownames(mean_df)[rownames(mean_df) == "11"] <- "Mean"
  
  return(mean_df)
}


# Summary of mean_Csweep function // averages within the k-fold optimisation parameters (summary_Csweep)
#
# Outputs:
#      data frame of the averaged value of all measures of performances for all folds from stats_all()
# Inputs:
#      kfold_all = a data frame containing all the measures of performance for each fold (from 'optimise_impute' function)
#
# Used in ImputationBlinded (Chunk 2) to solve for the mean prediction scores for each parameter 

summary_Csweep <- function(kfold_all){
  summary_df <- data.frame()
  summary <- mean_Csweep(kfold_all, stat_type=6)
  row_name <- names(kfold_all[[1]][[1]][6])
  summary_df <- summary[11,]
  rownames(summary_df)[rownames(summary_df) == "Mean"] <- row_name
  
  stat_type_list <- c(7:length(kfold_all[[1]][[1]]))
  
  for (i in stat_type_list){
    summary <- mean_Csweep(kfold_all, stat_type=i)
    row_name <- names(kfold_all[[1]][[1]][i])
    summary_df[nrow(summary_df)+1,] <- summary[11,]
    rownames(summary_df)[rownames(summary_df) == "Mean"] <- row_name
  }
  return(summary_df)
}


# Prediction scores for statistical analysis of best thresholds per cluster (optimse_impute)
# 
# Outputs: (optimse_impute)
#      list of prediction scores from a 10-fold validation process. Scores are grouped by a parameter pair
#      consisting of a total cluster and a threshold value. Each threshold value will be applied to each total cluster value
# Inputs:
#      cutClust = a list of the total clusters (extracted lengths from cutClust) that forms cl_list (a list of unique clusters)
#      thresh = a numerical list denoting the target threshold value from 0 to 1
#      cutClustvalues_dynamic = the output of the cl_cut_dynamic() function which gives a table of GeneIDs assigned to clusters
#      normAggLog = normalized gene normAggLog from RNA seq
#      GO_table = the Gene Onltology matrix for all GeneIDs
#      
# Used in Imputation Blinded (chunk 1) within the sample k_fold stage

optimise_impute <- function(cl_list, thresh, cutClustvalues_dynamic, normAggLog, GO_table)
  {
  scores <- list()
  
  for (i in cl_list){
    print(i)
    for (j in thresh){
      print(j)
      cl_tot <- as.character(i)
      cl <- cutClustvalues_dynamic[[cl_tot]][["GeneID_assignments"]]
      
      corrClust <- corr_per_clust(normAggLog, cl, i)
      # For faster runtime, set cluster 2 to 0
      corrClust[2] <- 0
      
      GOperClust <- GO_per_cl(GO_table, cl, i)
      GOperClustall <- GO_per_cl_list(GOperClust, i)
      
      sc <- cross_val(n=10, GO_table=GO_table, cutClust=cl, GOperClust=GOperClust, GOperClustall=GOperClustall, corrClust=corrClust, clust_total=i, thresh=j)
      
      scores[[paste0(i, "_", j)]] <- sc
    }
  }
  return(scores)
}


# A set of thresholds optimised for the cross validation process
#
# Outputs:
#      list of mean performance values for a series of thresholds determined by a set interval (thresh)
# Inputs: 
#      interval function = interval from 0 to 1 used to generate a sequence of values that will serve as a threshold [set manually for now]
#      GO_table = original GO annotation matrix
#      cutClust = matrix of genes w/ cluster number
#      GOperClust = GO terms per cluster from GO_per_cl
#      GOperClustall = a nested list object containing the GO terms (column names) before blinding in each cluster
#      corrClust = correlation values grouped by cluster
#      clust_total = total number of clusters (derived from cutClustvalues_dynamic)
#      thresh = threshold calculated from 0 to 1 [will be optimised]
#
# Used in Imputation (chunk 2) 

cross_val_thresh <- function(interval){
  
  kfold_global <- list()
  
  thresh <- seq(from = 0, to = 1, by = interval)
  
  for (thresh in thresh){
    
    kfold <- cross_val(n=10, GO_table, cutClust, GOperClust,
                       GOperClustall, corrClust,
                       clust_total, thresh)
    
    kfold_df <- as.data.frame(do.call(rbind, kfold))
    
    Sensitivity <- mean(unlist(kfold_df$Sensitivity))
    Specificity <- mean(unlist(kfold_df$Specificity))
    Precision <- mean(unlist(kfold_df$Precision))
    Accuracy <- mean(unlist(kfold_df$Accuracy))
    F1_Score <- mean(unlist(kfold_df$`F1 Score`))
    
    
    kfold_global[[paste0(thresh)]][["DF"]] <- kfold_df
    kfold_global[[paste0(thresh)]][["Sensitivity"]] <- Sensitivity
    kfold_global[[paste0(thresh)]][["Specificity"]] <- Specificity
    kfold_global[[paste0(thresh)]][["Precision"]] <- Precision
    kfold_global[[paste0(thresh)]][["Accuracy"]] <- Accuracy
    kfold_global[[paste0(thresh)]][["F1_Score"]] <- F1_Score
  }
  return(kfold_global)
}
