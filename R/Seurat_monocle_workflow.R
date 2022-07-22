# Namespaces
.AA_NAMESPACE = function() c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
.DNA_NAMESPACE = function() c('A', 'T', 'G', 'C')
.RNA_NAMESPACE = function() c('A', 'U', 'G', 'C')

#===============================================================================

#' Create a Seurat object from score matrix, then, dimensionality reduction and
#' classification of the score matrix.
#'
#' @param score_matrix Matrix. Score matrix of each site.
#' @param nfeatures Number of features to select as top variable features.
#' @param dims_SNN Dimensions of reduction to use as input with SNN.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param dims_UMAP Dimensions of reduction to use as input with UMAP.
#' @param pca.reduction.name dimensional reduction name with PCA, Default:"pca".
#' @param pca.reduction.key dimensional reduction key, specifies the string before the number for the dimension names with PCA. Default:"PC_".
#' @param umap.reduction.key dimensional reduction key, specifies the string before the number for the dimension names with UMAP. Default:"UMAP_".
#' @param ... Other parameters of CreateSeuratObject, FindVariableFeatures, ScaleData, RunPCA, FindNeighbors, FindClusters and RunUMAP functions in Seurat package.
#' @import Seurat
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readRNAStringSet
#' @importFrom Biostrings readAAStringSet
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings RNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings consensusMatrix
#' @importFrom methods as
#' @importFrom methods new
#' @importFrom methods slotNames
#' @importFrom stats aggregate
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats prcomp
#' @importFrom stats quantile
#' @importFrom stats wilcox.test
#' @importFrom utils combn
#' @importFrom utils data
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom VGAM negbinomial.size
#' @importFrom VGAM tobit
#' @importFrom VGAM uninormal
#' @importFrom xtable label
#' @importFrom SeuratObject VariableFeatures
#' @return score matrix of each site. Rownames of the repeated sequence name will be suffixed.
#' @seealso \code{\link{CreateSeuratObject},\link{FindVariableFeatures},\link{ScaleData},
#'                \link{RunPCA},\link{FindNeighbors},\link{FindClusters},\link{RunUMAP}}
#' @export
#' @examples
#' require(ggseqlogo)
#' data(ggseqlogo_sample)
#' test_aa <- get_alignment_score(alignment = seqs_aa[[1]],type = "AA")
#' my_Seurat <- Seurat_workflow(test_aa)
#-------------------------------------------------------------------------------

##Perform dimensionality reduction processing on a given data set,
##and the input object is score matrix
Seurat_workflow <- function(score_matrix,nfeatures = 2000,dims_SNN = 1:5,
                            resolution = 0.1,dims_UMAP = 1:5,
                            pca.reduction.name = "pca",
                            pca.reduction.key = "PC_",
                            umap.reduction.key = "UMAP_",...){
  if (is.null(rownames(score_matrix))){
    rownames(score_matrix) <- 1:nrow(score_matrix)
  }
  if (is.null(colnames(score_matrix))){
    colnames(score_matrix) <- 1:ncol(score_matrix)
  }
  score <- as(as.matrix(t(score_matrix)), "dgCMatrix")
  my_seurat <- Seurat::CreateSeuratObject(counts = score,...)
  my_seurat <- Seurat::FindVariableFeatures(my_seurat, nfeatures = nfeatures,...)
  all.genes <- rownames(my_seurat)
  my_seurat <- Seurat::ScaleData(my_seurat, features = all.genes,...)
  my_seurat <- Seurat::RunPCA(my_seurat, features = VariableFeatures(object = my_seurat),reduction.name = pca.reduction.name,
							  reduction.key = pca.reduction.key,...)
  my_seurat <- Seurat::FindNeighbors(my_seurat, dims = dims_SNN,...)
  my_seurat <- Seurat::FindClusters(my_seurat, resolution = resolution,...)
  my_seurat <- Seurat::RunUMAP(my_seurat, dims = dims_UMAP,reduction.key = umap.reduction.key,...)
  #my_seurat <- Seurat::RunTSNE(my_seurat, dims = 1:5)
  return(my_seurat)
}
#===============================================================================
#' import metadata from data frame
#'
#' @param seurat Seurat object.
#' @param metadata Number of features to select as top variable features.
#' @import Seurat
#' @return Seurat object with imported metadata
#' @export
#-------------------------------------------------------------------------------
import_metadata <- function(seurat,metadata){
	stopifnot("metadata must have the same rows" =
              nrow(seurat@meta.data) == nrow(metadata))
  if (is.null(colnames(seurat)) | is.null(metadata)){
    warning("no rownames in seurat object or new metadata")}
    else{
	    if (any(colnames(seurat) != rownames(metadata))){
	      warning("rownames of new metadata do not match rownames of old metadata")}
      if (any(colnames(seurat@meta.data) %in% colnames(metadata))){
        warning(paste0("colnames ",
                        paste(colnames(seurat@meta.data)[colnames(seurat@meta.data) %in% colnames(metadata)],collapse = ","),
                        " of new metadata appeared in colnames of old metadata"))
            }}
	seurat@meta.data <- cbind(seurat@meta.data,metadata)
	return(seurat)
}
#===============================================================================
#' Import a seurat or scatter/scran CellDataSet object and convert it to a monocle cds.
#' @description This function is modified by the \code{\link{importCDS}} function.
#' @param otherCDS The object you would like to convert into a monocle cds.
#' @param import_all Whether or not to import all the slots in seurat or scatter. Default is FALSE (or only keep minimal dataset).
#' @import Seurat monocle
#' @return a new monocle cell dataset object converted from other objects (Scatter or Seurat).
#' @export
#-------------------------------------------------------------------------------
##Import the objects generated by the Seurat package into the Monocle package
monocle_import <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts

    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }

    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    },
    #warning = function(w) { },
    error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)

      message("This Seurat object doesn't provide any meta data");
      pd
    })

    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0

    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)

    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL

        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS

      } else {
        # mist_list <- list(ident = ident)
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }
    var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)
    monocle_cds@auxClusteringData$seurat <- mist_list

  }
  return(monocle_cds)
}
#===============================================================================
#' Dimensionality reduction, classification and pseudotime analysis of the score matrix.
#'
#' @param seurat Seurat object.
#' @param method A character string specifying the algorithm to use for dimensionality reduction.
#' @param ... Other parameters of reduceDimension and orderCells functions in monocle package.
#' @import monocle
#' @return CellDataSet object with dimensionality reduction, classification and pseudotime analysis results.
#' @seealso \code{\link{reduceDimension},\link{orderCells}}
#' @export
#-------------------------------------------------------------------------------
##Pseudotime analysis applying monocle
monocle_workflow <- function(seurat,method = "DDRTree",ncenter = 5,...){
  monocle <- monocle_import(seurat)
  if (any(is.na(monocle$Size_Factor))){monocle$Size_Factor <- rep(1,length(monocle$Size_Factor))}
  monocle@expressionFamily@vfamily <- monocle@expressionFamily@vfamily[1]
  monocle  <- monocle::reduceDimension(monocle, max_components = 2,
                                       reduction_method = method,ncenter = ncenter,...)
  monocle  <- monocle::orderCells(monocle,...)
  return(monocle)
}
#===============================================================================
#' Find specific sites for each gene family.
#'
#' @param seurat Seurat object.
#' @param ident Character. One of colnames of Seurat object metadata. If NULL, use the active.ident of the Seurat object as the grouping information.
#' @param only.pos Only return positive marker site. Default:TRUE.
#' @param logfc Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default:1.
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default: 0.25.
#' @param ... Other parameters of Seurat::FindAllMarkers.
#' @import Seurat
#' @return Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
#' @seealso \code{\link{FindAllMarkers}}
#' @export
#-------------------------------------------------------------------------------
find_Specific_site <- function(seurat,ident = NULL,only.pos = TRUE,logfc = 0.1,min.pct = 0,...){
  if (is.null(ident)){
    all_markers <- Seurat::FindAllMarkers(seurat,only.pos = only.pos,logfc.threshold = logfc,min.pct = min.pct,...)
  } else{
    stopifnot("Please imput a Seurat object" =
              "Seurat" %in% class(seurat))
    stopifnot("ident must be length 1 vector" =
              length(ident) == 1)
    stopifnot("ident do not in metadata" =
              ident %in% colnames(seurat@meta.data))
    metadata <- seurat@meta.data
    my_ident <- as.factor(metadata[,ident])
    names(my_ident) <- rownames(metadata)
    seurat@active.ident <- my_ident
    all_markers <- Seurat::FindAllMarkers(seurat,only.pos = only.pos,logfc.threshold = logfc,min.pct = min.pct,...)
  }
  return(all_markers)
}
#===============================================================================
#' Get bits of each site.
#'
#' @param seq_prob Seurat object.
#' @param type Character. One of DNA, RNA and AA.
#' @param freq Logical. If TRUE, the final result will be multiplied by the ratio of valid characters. Default:TRUE.
#' @param base Logical. If TRUE, only keep the score of AA_STANDARD/DNA_BASES/RNA_BASES for protein/DNA/RNA sequences. Default:TRUE.
#' @param prob Logical. If TRUE, return seq prob of each site. Default:FALSE.
#' @return bits/seq prob of each site
#-------------------------------------------------------------------------------
##get_seq_bits
get_bits <- function(seq_prob,type = c("DNA","RNA","AA"),freq = T,base = T,prob = F){
  if (base){
    seq_prob <- switch(type,DNA = seq_prob[DNA_BASES,],
                       RNA = seq_prob[RNA_BASES,],
                       AA = seq_prob[AA_STANDARD,])
  }
  if (!prob){
    bits <- switch(type,DNA = apply(seq_prob,2,function(x)(2 + sum(x[x > 0]*(log2(x[x > 0]))))),
                   RNA = apply(seq_prob,2,function(x)(2 + sum(x[x > 0]*(log2(x[x > 0]))))),
                   AA = apply(seq_prob,2,function(x)(log2(20) + sum(x[x > 0]*(log2(x[x > 0]))))))
  } else {
    bits <- seq_prob
  }
  if (freq){
    bits <- bits * colSums(seq_prob)
  }
  bits
}
#===============================================================================
#' Get sequence conservation of each site.
#'
#' @param seqs input sequence alignment result, for DNA/RNA/protein sequences, DNAStringSet/RNAStringSet/AAStringSet or character class were required. All sequences must have the same length.
#' @param alignment_file The path of input sequence alignment result.
#' @param type Character. One of DNA, RNA and AA.
#' @param group Vector. Grouping information of each sequence. The length of group vector is consistent with the number of sequences.
#' @param ... Other parameters of \code{\link{consensusMatrix}}.
#' @return sequence conservation matrix of each site.
#' @export
#-------------------------------------------------------------------------------
get_seq_conservation <- function(seqs = NULL,alignment_file = NULL,type = c("DNA","RNA","AA"),group = NULL,...){
  if (is.null(seqs)){
    my_seq_alignment <- switch(type,DNA = Biostrings::readDNAStringSet(alignment_file),
                               RNA = Biostrings::readRNAStringSet(alignment_file),
                               AA = Biostrings::readAAStringSet(alignment_file))
  } else if (is.null(alignment_file)){
    if ('character' %in% class(seqs)){
      my_seq_alignment <- switch(type,DNA = Biostrings::DNAStringSet(seqs),
                                 RNA = Biostrings::RNAStringSet(seqs),
                                 AA = Biostrings::AAStringSet(seqs))
    } else if(any(c('AAStringSet','DNAStringSet','RNAStringSet') %in% class(seqs))){
      my_seq_alignment <- seqs
    } else {stop('Wrong imput sequence type')}
  } else {stop("No imput sequences")}
  if (is.null(group)){
    my_seq_prob <- Biostrings::consensusMatrix(my_seq_alignment,as.prob = TRUE,...)
    my_seq_prob <- switch(type,DNA = my_seq_prob[DNA_BASES,],
                          RNA = my_seq_prob[RNA_BASES,],
                          AA = my_seq_prob[AA_STANDARD,])
    my_seq_bits <- switch(type,DNA = get_bits(my_seq_prob,type = "DNA",freq = T),
                          RNA = get_bits(my_seq_prob,type = "RNA",freq = T),
                          AA = get_bits(my_seq_prob,type = "AA",freq = T))

  } else {
    stopifnot("grouped seqs is not equal to aligned seqs" =
                length(group) == length(my_seq_alignment))
    if (any(is.na(group))){warning("The grouping information of some objects is invalid")}
    alignment_list <- split(my_seq_alignment,group)
    my_seq_prob <- lapply(alignment_list,function(x)
      Biostrings::consensusMatrix(x,as.prob = TRUE,...))
    my_seq_bits <- switch(type,DNA = do.call(rbind,lapply(my_seq_prob,function(x)get_bits(x,type = "DNA",freq = T,...))),
                          RNA = do.call(rbind,lapply(my_seq_prob,function(x)get_bits(x,type = "RNA",freq = T,...))),
                          AA = do.call(rbind,lapply(my_seq_prob,function(x)get_bits(x,type = "AA",freq = T,...))))
  }
  return(my_seq_bits)
}
#===============================================================================
## scEpath
## edmonds algorithm
## This function is transplanted from the Matlab package scEpath.
#' Get max value incoming edge and its index from edge matrix.
#'
#' @param G Edge matrix.
#' @param VERTEX The node number in G.
#' @return A list contain EDGE_VALUE and EDGE_ADDED
#-------------------------------------------------------------------------------
index_of_max_value_incoming_edge <- function(G,VERTEX){
  if (length(G) == 0){
    EDGE_VALUE = -1
    EDGE_INDEX = NaN
    result <- list(EDGE_VALUE = EDGE_VALUE,EDGE_ADDED = EDGE_ADDED)
    return(result)
  }
  INDICES <- which(G[,2] == VERTEX)
  VALUES <- G[,3][INDICES]
  EDGE_VALUE <- min(VALUES,na.rm = T)
  EDGE_ADDED <- INDICES[which.min(VALUES)]
  result <- list(EDGE_VALUE = EDGE_VALUE,EDGE_ADDED = EDGE_ADDED)
  return(result)
}
#===============================================================================
## This function is transplanted from the Matlab package scEpath.
#' Get shortest paths from edge matrix.
#'
#' @param G Edge matrix.
#' @param S Numeric. The Start node number in G.
#' @param D Numeric. The End node number in G.
#' @importFrom igraph shortest.paths
#' @importFrom igraph add_edges
#' @importFrom igraph make_empty_graph
#' @return A list contain dist and path between Start node number and End node number.
#-------------------------------------------------------------------------------
iscycle <- function(G,S,D){
  if (!is.null(G)){
    MAXN <- max(c(S,D,max(G[,1:2])))
    G <- rbind(G,c(MAXN,MAXN,0))
    G2 <- add_edges(make_empty_graph(MAXN), t(G[,1:2]), weight=G[,3])
    dist <- c(shortest.paths(G2,S,D))
    path <- unlist(shortest_paths(G2,S,D)$vpath)
    if (is.infinite(dist) | length(path) < 2){
      dist <- c(shortest.paths(G2,D,S))
      path <- unlist(shortest_paths(G2,D,S)$vpath)
    }
  } else {
    dist <- Inf
    path <- numeric()
  }
  result <- list(dist = dist,path = path)
  return(result)
}
#===============================================================================
## This function is transplanted from the Matlab package scEpath.
#' obtaining relabelled graphs.
#'
#' @param ED List contain graphical structure.
#' @return returns the new graphical structure and the mapping between old and new graphs.
#-------------------------------------------------------------------------------
relabelgraph <- function(ED){
  NODES <- ED$V
  UNSHRUNK <- setdiff(NODES,ED$VERTICESINCKT)
  SHRUNK <- ED$VERTICESINCKT
  destcount2 <- 2
  MAPVERT <- matrix(0,length(NODES),2)
  for (i in 1:length(NODES)){
    if (NODES[i] %in% SHRUNK){
      MAPVERT[i,] <- c(NODES[i],1)
    } else {
      MAPVERT[i,] <- c(NODES[i],destcount2)
      destcount2 = destcount2 + 1
    }
  }
  NEW_GRAPH <- list()
  NEW_GRAPH$V <- sort(unique(MAPVERT[,2]))
  destcount <- 1
  NEW_GRAPH$E <- list()
  MAPEDGE <- list()
  for (i in 1:nrow(ED$E)){
    ARC <- unlist(ED$E[i,])
    SP <- ARC[1];EP <- ARC[2];WT <- ARC[3]
    if ((!SP %in% SHRUNK) & (!EP %in% SHRUNK)){
      ISP <- MAPVERT[which(MAPVERT[,1] == SP),2]
      IEP <- MAPVERT[which(MAPVERT[,1] == EP),2]
      NEW_GRAPH$E[destcount] <- list(c(ISP,IEP,WT))
      if (i %in% ED$BE){
        NEW_GRAPH$BE <- c(NEW_GRAPH$BE,destcount)
      }
      MAPEDGE[i] <- list(c(i,destcount))
      destcount=destcount+1

    }else if ((SP %in% SHRUNK) & (!EP %in% SHRUNK)){
      ISP <- MAPVERT[which(MAPVERT[,1] == SP),2]
      IEP <- MAPVERT[which(MAPVERT[,1] == EP),2]
      NEW_GRAPH$E[destcount] <- list(c(ISP,IEP,WT))
      if (i %in% ED$BE){
        NEW_GRAPH$BE <- c(NEW_GRAPH$BE,destcount)
      }
      MAPEDGE[i] <- list(c(i,destcount))
      destcount=destcount+1
    }else if ((!SP %in% SHRUNK) & (EP %in% SHRUNK)){
      ISP <- MAPVERT[which(MAPVERT[,1] == SP),2]
      IEP <- MAPVERT[which(MAPVERT[,1] == EP),2]
      INCIDENTCKTWT <- finde0e1(ED = ED,EP = EP)
      #print(INCIDENTCKTWT)
      NEWWT <- WT - INCIDENTCKTWT
      NEW_GRAPH$E[destcount]=list(c(ISP,IEP,NEWWT))
      if (i %in% ED$BE){
        NEW_GRAPH$BE <- c(NEW_GRAPH$BE,destcount)
      }
      MAPEDGE[i] <- list(c(i,destcount))
      destcount=destcount+1
    }else {
      MAPEDGE[i] <- list(c(i,-1))
    }
  }
  MAPEDGE <- do.call(rbind,MAPEDGE)
  count = 1;NEW_GRAPH$BV <- c()
  for (i in 1:length(ED$BV)){
    if (!ED$BV[i] %in% SHRUNK){
      IND <- which(MAPVERT[,1] == ED$BV[i])
      NEW_GRAPH$BV[count] <- MAPVERT[IND,2]
      count=count+1
    }
  }
  NEW_GRAPH$BV <- unlist(NEW_GRAPH$BV)
  if (count > 1){
    NEW_GRAPH$BV <- sort(unique(NEW_GRAPH$BV))
  }
  NEW_GRAPH$E <- as.matrix(do.call(rbind,NEW_GRAPH$E))
  result <- list(NEW_GRAPH = NEW_GRAPH,
                 MAPVERT = MAPVERT,
                 MAPEDGE = MAPEDGE)
}
#===============================================================================
## This function is transplanted from the Matlab package scEpath.
#' Extract the weight of the end node of edges which are in circuit.
#'
#' @param ED List contain graphical structure.
#' @param EP Numeric. The end node number in graphs.
#' @return the weight of the end node of edges which are in circuit.
#-------------------------------------------------------------------------------
finde0e1 <- function(ED,EP){
  INBE <- ED$E[ED$BE,,drop = F]
  count = 1
  ISCKTEDGE <- list()
  for (i in 1:length(ED$BE)){
    if (all(INBE[i,1:2] %in% ED$VERTICESINCKT)){
      ISCKTEDGE[count] <- list(INBE[i,])
      count <- count + 1
    }
  }
  ISCKTEDGE <- do.call(rbind,ISCKTEDGE)
  min_wt_in_ckt <- min(ISCKTEDGE[,3])
  INDEXofINCIDENT <- which(ISCKTEDGE[,2]==EP)
  incident_wt_in_ckt <- ISCKTEDGE[INDEXofINCIDENT,3]
  return(incident_wt_in_ckt)
}
#===============================================================================
#' Getting circuit with min dist.
#'
#' @param ED List contain graphical structure.
#' @return circuit with min dist.
#-------------------------------------------------------------------------------
reconstruct_2 <- function(ED){
  MAXI <- length(ED)
  for (i in MAXI:2){
    RECONBUCKET <- c()
    EXPANDED <- 0
    for (j in 1:length(ED[[i]]$BE)){
      INDEX <- which(ED[[i-1]]$MAPPINGEDGE[,2] == ED[[i]]$BE[j])
      TOBEADDED1 <- ED[[i-1]]$MAPPINGEDGE[INDEX,1]
      RECONBUCKET <- c(RECONBUCKET,TOBEADDED1)
    }
    INDICES_OF_CKT <- find_ckt_edges_to_add(ED,i-1,RECONBUCKET);
    RECONBUCKET <- c(RECONBUCKET,INDICES_OF_CKT)
    ED[[i-1]]$BE <- RECONBUCKET
  }
  TREEMAX <- ED[[1]]$BE
  return(TREEMAX)
}
#===============================================================================
#' Getting circuit with min dist.
#'
#' @param ED List contain graphical structure.
#' @param level Number of iterations
#' @param RECONBUCKET edge in circuit
#-------------------------------------------------------------------------------
find_ckt_edges_to_add <- function(ED,level,RECONBUCKET){
  count = 1
  IND_CKT <- c()
  for (i in 1:length(ED[[level]]$BE)){
    if (all(ED[[level]]$E[ED[[level]]$BE[i],1:2] %in% ED[[level]]$VERTICESINCKT)){
      IND_CKT <- c(IND_CKT,ED[[level]]$BE[i])
    }
  }
  #now find the minimum
  GRAPH_TILL_NOW <- ED[[level]]$E[RECONBUCKET,,drop = F]
  NODES_WITH_OUTGOING <- sort(unique(GRAPH_TILL_NOW[,1]))
  NODES_WITH_INCOMING <- sort(unique(GRAPH_TILL_NOW[,2]))

  IS_OUTTREE <- length(intersect(NODES_WITH_OUTGOING,ED[[level]]$VERTICESINCKT))
  IS_INTREE <- length(intersect(NODES_WITH_INCOMING,ED[[level]]$VERTICESINCKT))
  INCOMING_INTERFACEPOINT <- intersect(NODES_WITH_INCOMING,ED[[level]]$VERTICESINCKT)

  CKT <- ED[[level]]$E[IND_CKT,,drop = F]

  if (IS_INTREE > 0){
    b <- which(CKT[,2] == INCOMING_INTERFACEPOINT)

  } else {
    b <- which.min(CKT[,3])
    a <- min(CKT[,3])
  }
  TO_NOT_BE_INCLUDED <- IND_CKT[b]
  INDICES <- setdiff(IND_CKT,TO_NOT_BE_INCLUDED)
  return(INDICES)
}
#===============================================================================
#' Getting circuit with min dist.
#'
#' @param V Nodes of graphs.
#' @param E Edge of graphs. Containing weights.
#' @param rootNode Numeric. Rootnode of graphs.
#' @return circuit with min dist.
#-------------------------------------------------------------------------------
edmonds <- function(V,E,rootNode){
  edmonds_result <- list()
  ED <- list()
  ED$BV <- rootNode
  ED$BE <- numeric()
  ED$V <- V
  ED$E <- E
  CURRENT_i=1
  while (1){
    V <- ED$V
    E <- ED$E
    VERTICES_NOT_IN_BV <- setdiff(V,ED$BV)
    if(length(VERTICES_NOT_IN_BV) == 0){
      edmonds_result[CURRENT_i] <- list(ED)
      return(edmonds_result)
      break
    }
    VERTEX_ADDED <- VERTICES_NOT_IN_BV[1]
    ED$BV <- c(ED$BV,VERTEX_ADDED)
    EDGE_VALUE <- index_of_max_value_incoming_edge(G = E,VERTEX = VERTEX_ADDED)$EDGE_VALUE
    EDGE_ADDED <- index_of_max_value_incoming_edge(G = E,VERTEX = VERTEX_ADDED)$EDGE_ADDED
    if (EDGE_VALUE > 0){
      dist <- iscycle(G = E[ED$BE,,drop = F],S = E[EDGE_ADDED,1],D = E[EDGE_ADDED,2])$dist
      path <- iscycle(G = E[ED$BE,,drop = F],S = E[EDGE_ADDED,1],D = E[EDGE_ADDED,2])$path
      if (length(path) >= 2){
        ED$VERTICESINCKT <- path
      }
      ED$BE <- c(ED$BE,EDGE_ADDED)
      if (is.finite(dist)){
        tmp <- relabelgraph(ED)
        NEW_GRAPH <- tmp$NEW_GRAPH
        MAPVERT <- tmp$MAPVERT
        MAPEDGE <- tmp$MAPEDGE
        ED$MAPPINGVERT <- MAPVERT
        ED$MAPPINGEDGE <- MAPEDGE
        edmonds_result[CURRENT_i] <- list(ED)
        #print(CURRENT_i)
        CURRENT_i <- CURRENT_i + 1
        ED <- list()
        ED$BV <- NEW_GRAPH$BV
        ED$BE <- NEW_GRAPH$BE
        ED$V <- NEW_GRAPH$V
        ED$E <- NEW_GRAPH$E
      }
    }

  }
}
#===============================================================================
#' Getting gene lineage from score matrix.
#'
#' @param data Matrix. Score matrix of each site.
#' @param range Numeric. Extract the specified column of the score matrix for network construction.
#' @param thresh Numeric. Edges whose correlation is less than thresh will be filtered out.
#' @param normalization Logical. To normalize the data or not. Default:TRUE.
#' @param group Factor. Grouping information of each gene. Must contain the same length of nrow(data).
#' @param alpha significance level for a two-sided Wilcoxon rank-sum test of the distribution of scEnergy between two metagenes.Default:0.01.
#' @param theta1 Numeric. Select the genes which include theta1 energy in each cluster. Default:0.8.
#' @param rootNode Numeric. The root node (i.e., start state) for inferring gene family lineages.
#' @importFrom igraph add_edges
#' @importFrom igraph make_empty_graph
#' @importFrom igraph shortest_paths
#' @return a struct giving the gene lineage information.
#' @export
#-------------------------------------------------------------------------------
constructingNetwork <- function(data,range = NULL,thresh = 0.1,normalization = T,
                                group = NULL,alpha = 0.01,theta1 = 0.8,
                                rootNode = 1){
  if (!is.null(range)){
    data <- data[,range]
  }
  data_sd <- apply(data,2,sd)
  data <- data[,data_sd > 0]
  total_node <- ncol(data)
  data1 <- cor(data)
  deg <- rowSums(matrix(data1 > thresh,nrow(data1),ncol(data1)))
  data1 <- abs(upper.tri(data1))
  data2 <- matrix(data1 > thresh,nrow(data1),ncol(data1))
  data_edge <- which(data2,arr.ind = T)
  data_edge <- data_edge[order(data_edge[,1]),]
  num_edges <- nrow(data_edge)
  num_nodes <- length(unique(c(data_edge)))
  ID_select <- unique(c(data_edge))
  data_filter <- data[,ID_select]
  data_filter <- (data_filter - min(data_filter)) / (max(data_filter) - min(data_filter))
  ##load network information
  dataSumNeibor <- matrix(0,nrow(data_filter),ncol(data_filter))
  for (i in 1:ncol(dataSumNeibor)){
    a <- unique(c(data_edge[data_edge[,"row"] == i | data_edge[,"col"] == i,]))
    dataSumNeibor[,i] <- rowSums(data_filter[,a])
  }
  ##calculate the scEnergy for each cell
  scE = -data_filter * log(data_filter/dataSumNeibor,base = exp(1))
  scE[is.nan(scE)] <- 0
  scE[is.infinite(scE)] <- 0
  scEcell <- rowSums(scE)
  if (normalization){
    scEcell = ((scEcell/mean(scEcell))^2)/(1+(scEcell/mean(scEcell))^2)
  }
  scE.pca <- prcomp(scE,scale=T,rank=10,retx=T)
  ydata <- scE.pca$x[,1:2]
  y <- sort(unique(group))
  numCluster <- length(y)
  #(1) find the core cell in each cluster
  centroidCoreCell <- matrix(0,numCluster,2)
  FCoreCell <- matrix(0,numCluster,1)
  cell_dataframe <- data.frame(scE = scEcell,group = group,
                               PC1 = ydata[,1],PC2 = ydata[,2])
  cell_dataframe_sort <- cell_dataframe[order(cell_dataframe$scE),]
  cell_dataframe_list <- split(cell_dataframe_sort,cell_dataframe$group)
  cell_dataframe_list <- cell_dataframe_list[y]
  if (theta1 < 1){
    idxCoreCellCi <- lapply(cell_dataframe_list,
                            function(x)which((cumsum(x$scE) / sum(x$scE)) > theta1)[1])
    } else {
    idxCoreCellCi <- lapply(cell_dataframe_list,nrow)
    }
    idxCoreCellCi <- idxCoreCellCi[y]
    centroidCoreCell <- data.frame(group = y,
                                   PC1 = unlist(lapply(1:length(idxCoreCellCi),
                                                function(x)mean(cell_dataframe_list[[x]]$PC1[1:(idxCoreCellCi[[x]])]))),
                                   PC2 = unlist(lapply(1:length(idxCoreCellCi),
                                                function(x)mean(cell_dataframe_list[[x]]$PC2[1:(idxCoreCellCi[[x]])]))))
    #Tukey's trimean
    FCoreCell <- unlist(lapply(1:length(idxCoreCellCi),
                function(x)sum(0.25*quantile(cell_dataframe_list[[x]]$scE[1:(idxCoreCellCi[[x]])]
                                          ,c(0.25,0.5,0.75)) +
                            0.25*quantile(cell_dataframe_list[[x]]$scE[1:(idxCoreCellCi[[x]])]
                                          ,c(0.5)))))

    ##(2)calculate the Transition Probability (TP)
    ##(2.1)based on distance in low-dimension space
    D1 <- dist(as.matrix(centroidCoreCell[,2:3]),upper = T)
    epsilon <- 1/2*max(D1)
    P1 <- as.matrix(exp(-D1^2/epsilon^2))
    wx <- colSums(P1) #sum of all weights of the vertices
    wg <- sum(wx) #sum of all weights of the vertices
    sd <- wx/wg #stationary distribution
    #transition matrix
    P1 <- P1 / wx
    #symmetric transition matirx
    P1 <- P1 * sqrt(sd) / sqrt(rep(sd,each = ncol(P1)))
    #(2.2) based on the energy
    #calculate the probability of each matacell
    PCorecell <- exp(-FCoreCell)/sum(exp(-FCoreCell))
    P2 <- PCorecell
    #(2.3) define the likehood (probability) of transition by combine these two
    TP <- (diag(1-P2) %*% P1) + diag(P2)
    #(3) construct a probabilistic directed graph using the energy information
    C <- t(combn(numCluster,2))
    A <- matrix(0,numCluster,numCluster)
    for (i in 1:nrow(C)){
      p <- wilcox.test(cell_dataframe_list[[C[i,1]]]$scE,cell_dataframe_list[[C[i,2]]]$scE)$p.value
      if ((FCoreCell[C[i,1]] > FCoreCell[C[i,2]]) & (p < alpha)){
        A[C[i,1],C[i,2]] <- 1 - TP[C[i,1],C[i,2]]
      } else {
        A[C[i,2],C[i,1]] <- 1 - TP[C[i,2],C[i,1]]
        A[C[i,1],C[i,2]] <- 1 - TP[C[i,1],C[i,2]]
      }
      #print(i)
    }
    #(4)find the minimal directed spanning tree
    V = 1:numCluster
    E = data.frame(which(A > 0,arr.ind = T))
    E$A <- A[A > 0]
    E <- E[order(E$row),]
    rownames(E) <- NULL
    colnames(E) <- NULL
    E <- as.matrix(E)
    GT <- edmonds(V,E,rootNode)
    MDSTedgeOrder <- reconstruct_2(GT)
    MDST <- GT[[1]]$E[MDSTedgeOrder,,drop = F]
    MDST[MDST[,2] == rootNode,] <- MDST[MDST[,2] == rootNode,c(2,1,3)]
    idxS <- rootNode
    idxT <- setdiff(1:numCluster,MDST[,1])
    my_network <- add_edges(make_empty_graph(numCluster), t(MDST[,1:2]), weight=MDST[,3])
    path <- list()
    for (i in 1:length(idxT)){
      path[i] <- list(unlist(shortest_paths(my_network,idxS,idxT[i],weights = NA)$vpath))
    }
    MDST[,3] <- 1 - MDST[,3]
    MDST <- rbind(MDST,matrix(c(1:numCluster,1:numCluster,diag(TP)),ncol = 3))
    my_network <- add_edges(make_empty_graph(numCluster), t(MDST[,1:2]), weights=MDST[,3])
    lineageIfo <- list()
    lineageIfo$TP <- TP
    lineageIfo$MDST <- MDST
    lineageIfo$path <- path
    lineageIfo$PDG <- A
    lineageIfo$centroid <- centroidCoreCell
    lineageIfo$rootNode <- rootNode
    lineageIfo$group <- centroidCoreCell$group
    lineageIfo$net <- my_network
    return(lineageIfo)
}
#===============================================================================
#' Combine two alignment results.
#'
#' Solves (Needleman-Wunsch) global alignment problems on two alignment results.
#'
#' @param alignment1 Alignment result of group1.
#' @param alignment2 Alignment result of group2.
#' @param gapOpening Numeric. The cost for opening a gap in the alignment.
#' @param seqtype One of DNA, RNA and AA.
#' @param gapExtension Numeric. The incremental cost incurred along the length of the gap in the alignment.
#' @param substitutionMatrix Substitution matrix representing the fixed substitution scores for an alignment.
#' @param standard Logical. If TRUE, only keep the score of AA_STANDARD/DNA_BASES/RNA_BASES for protein/DNA/RNA sequences. Default:TRUE.
#' @param alignment_thres Numeric. Remove sites with amino acid ratio lower than alignment_thres. Default:0.
#' @return Combined multiple alignment results.
#' @export
#-------------------------------------------------------------------------------
##This function combine multiple alignment results of multiple groups
group_alignment <- function(alignment1,alignment2,seqtype = "AA",
                            gapOpening=-4, gapExtension=-1,
                            substitutionMatrix = "BLOSUM62",
                            standard = T,alignment_thres = 0
){
  width1 <- unique(width(alignment1))##group1 width
  width2 <- unique(width(alignment2))##group2 width
  if (length(width1) != 1 | length(width2) != 1){
    stop("alignment result must contain the same width!")
  }
  ##Convert all sequences to uppercase
  ##Convert all sequences to matrix format
  seq1 <- as.matrix(AAStringSet(toupper(as.character(alignment1))))
  seq2 <- as.matrix(AAStringSet(toupper(as.character(alignment2))))
  ##Generate amino acid substitution matrix
  if (length(substitutionMatrix) == 1 & (substitutionMatrix %in% c("BLOSUM45","BLOSUM50",
                                                                   "BLOSUM62","BLOSUM80",
                                                                   "BLOSUM100","PAM30",
                                                                   "PAM40","PAM70",
                                                                   "PAM120","PAM250"))){
    substitutionMatrix <- get(data(list = substitutionMatrix))
  }
  ##Extract standard amino acids
  if (seqtype == "AA"){
    if (standard){
      substitutionMatrix <- substitutionMatrix[AA_STANDARD,AA_STANDARD]
      my_alphabet <- AA_STANDARD
      seq_percent1 <- get_seq_percent(seq1,alphabet = my_alphabet)
      seq_percent2 <- get_seq_percent(seq2,alphabet = my_alphabet)
      ##Remove sites with low amino acid ratio
      seq_filter1 <- rowSums(seq_percent1) > alignment_thres
      seq_filter2 <- rowSums(seq_percent2) > alignment_thres
      seq1 <- seq1[,seq_filter1]
      seq2 <- seq2[,seq_filter2]
      seq_percent1 <- seq_percent1[seq_filter1,]
      seq_percent2 <- seq_percent2[seq_filter2,]
    } else {
    my_alphabet <- AA_ALPHABET
    seq_percent1 <- get_seq_percent(seq1,alphabet = my_alphabet)
    seq_percent2 <- get_seq_percent(seq2,alphabet = my_alphabet)
    }
  } else if (seqtype == "DNA"){
    if (standard){
      substitutionMatrix <- substitutionMatrix[DNA_BASES,DNA_BASES]
      my_alphabet <- DNA_BASES
      seq_percent1 <- get_seq_percent(seq1,alphabet = my_alphabet)
      seq_percent2 <- get_seq_percent(seq2,alphabet = my_alphabet)
      ##Remove sites with low amino acid ratio
      seq_filter1 <- rowSums(seq_percent1) > alignment_thres
      seq_filter2 <- rowSums(seq_percent2) > alignment_thres
      seq1 <- seq1[,seq_filter1]
      seq2 <- seq2[,seq_filter2]
      seq_percent1 <- seq_percent1[seq_filter1,]
      seq_percent2 <- seq_percent2[seq_filter2,]
    } else {
      my_alphabet <- DNA_ALPHABET
      seq_percent1 <- get_seq_percent(seq1,alphabet = my_alphabet)
      seq_percent2 <- get_seq_percent(seq2,alphabet = my_alphabet)
    }
  } else if (seqtype == "RNA"){
    if (standard){
      substitutionMatrix <- substitutionMatrix[RNA_BASES,RNA_BASES]
      my_alphabet <- RNA_BASES
      seq_percent1 <- get_seq_percent(seq1,alphabet = my_alphabet)
      seq_percent2 <- get_seq_percent(seq2,alphabet = my_alphabet)
      ##Remove sites with low amino acid ratio
      seq_filter1 <- rowSums(seq_percent1) > alignment_thres
      seq_filter2 <- rowSums(seq_percent2) > alignment_thres
      seq1 <- seq1[,seq_filter1]
      seq2 <- seq2[,seq_filter2]
      seq_percent1 <- seq_percent1[seq_filter1,]
      seq_percent2 <- seq_percent2[seq_filter2,]
    } else {
      my_alphabet <- RNA_ALPHABET
      seq_percent1 <- get_seq_percent(seq1,alphabet = my_alphabet)
      seq_percent2 <- get_seq_percent(seq2,alphabet = my_alphabet)
    }
  }
  ##Needleman-Wunsch algorithm calculates global alignment of multiple groups
  width1 <- nrow(seq_percent1)
  width2 <- nrow(seq_percent2)
  result <- matrix(0,width1+1,width2+1)
  path <- matrix(0,width1+1,width2+1)
  result[1,] <- c(0,seq(gapOpening,by = gapExtension,length.out = width2))
  result[,1] <- c(0,seq(gapOpening,by = gapExtension,length.out = width1))
  path[1,] <- 2
  path[,1] <- 1
  path[1,1] <- 0
  for (i in 2:(nrow(result))){
    for (j in 2:(ncol(result))){
      V1 <- ifelse((path[i,j-1] == 1),
                   result[i,j-1] + gapExtension,
                   result[i,j-1] + gapOpening)
      V2 <- ifelse((path[i-1,j] == 2),
                   result[i-1,j] + gapExtension,
                   result[i-1,j] + gapOpening)
      V3 <- result[i-1,j-1] + site_score(seq_percent1[i-1,],seq_percent2[j-1,],substitutionMatrix)
      path[i,j] <- which.max(c(V1,V2,V3))
      result[i,j] <- max(c(V1,V2,V3))
    }
  }
  ##Needleman-Wunsch algorithm backtracking to find an optimal solution
  alignment <- NULL
  n <- nrow(path);m <- ncol(path)
  while(1){
    if (path[n,m] == 3){
      tmp <- c(seq1[,n-1],seq2[,m-1])
      alignment <- cbind(alignment,tmp)
      n <- n - 1;m <- m - 1
    } else if (path[n,m] == 2){
      tmp <- c(rep("-",nrow(seq1)),seq2[,m-1])
      alignment <- cbind(alignment,tmp)
      m <- m - 1
    } else if (path[n,m] == 1){
      tmp <- c(seq1[,n-1],rep("-",nrow(seq2)))
      alignment <- cbind(alignment,tmp)
      n <- n - 1
    } else if (path[n,m] == 0){
      final_alignment <- AAStringSet(apply(alignment,1,function(x)paste0(rev(x),collapse = "")))
      names(final_alignment) <- c(names(alignment1),names(alignment2))
      seq_percent_all <-  get_seq_percent(alignment,alphabet = my_alphabet)
      ##return path, score_matrix, alignment_result, aa percent of each site
      final <- list(path = path,matrix = result,alignment = final_alignment,site_percent = seq_percent_all)
      return(final)
    }
  }
}
#===============================================================================
#' Calculate the match score between two sites.
#'
#' @param site1 score of site1.
#' @param site2 score of site2.
#' @param substitutionMatrix Substitution matrix representing the fixed substitution scores for an alignment.
#' @param quickly If TRUE, only keep the score of AA_STANDARD/DNA_BASES/RNA_BASES for protein/DNA/RNA sequences. Default:TRUE.
#' @return Combined multiple alignment results.
#' @export
#-------------------------------------------------------------------------------
site_score <- function(site1,site2,substitutionMatrix,quickly = F){
  if (quickly){
    site_score <- sum((site1 %*% t(site2)) * substitutionMatrix)
  } else {
    site_score_matrix <- matrix(0,nrow(substitutionMatrix),ncol(substitutionMatrix))
    colnames(site_score_matrix) <- rownames(site_score_matrix) <- rownames(substitutionMatrix)
    ##Calculate the matching probability of two sequences group
    for (i in rownames(site_score_matrix)){
      site_score_matrix[i,i] <- min(site1[i],site2[i])
    }
    mismatch <- 1 - sum(diag(site_score_matrix))
    ##Calculate the mismatching probability of two sequences group
    for (j in rownames(site_score_matrix)){
      for (k in colnames(site_score_matrix)){
        if (j != k){
          if ((site1[j] > site2[j]) & (site2[k] > site1[k])){
            site_score_matrix[j,k] <- (site1[j]-site2[j]) * (site2[k]-site1[k]) / mismatch
          } else {
            site_score_matrix[j,k] <- 0
          }
        } else if (j == k){
          site_score_matrix[j,k] <- site_score_matrix[j,k]
        }
      }
    }
    site_score <- sum(site_score_matrix * substitutionMatrix)
  }
  return(site_score)
}
#===============================================================================
#' Calculate the Percent matrix of each site.
#'
#' @param dat score matrix.
#' @param alphabet alphabet of dat.
#' @return Percent matrix of each site for each element in alphabet.
#-------------------------------------------------------------------------------
get_seq_percent <- function(dat,alphabet){

  result <- t(sapply(1:ncol(dat),
                     function(x)table(factor(dat[,x],
                                             levels = alphabet)))) / nrow(dat)
  return(result)
}
#===============================================================================
#' Combine multiple alignment results according to the distance of each alignment result.
#'
#' @param alignment_list A list containing multiple alignment results.
#' @param alignment_file A vector containing paths of multiple alignment results.
#' @param group_dist Dist. This parameter could be applied to import the calculated distance between groups. If NULL, distance between groups will be calculated automatically.
#' @param type Sequences type. One of DNA, RNA and AA.
#' @param tree Logical. If TRUE, also return hclust tree of each family. Default: FALSE.
#' @param ... Other parameter of group_alignment
#' @return multiple alignment results with hclust tree
#' @export
#-------------------------------------------------------------------------------
combine_alignment <- function(alignment_list = NULL,alignment_file = NULL,
                            group_dist = NULL,type = "AA",tree = F,...){
  if (is.null(alignment_list)){
    my_seq_alignment <- switch(type,DNA = lapply(alignment_file,Biostrings::readDNAStringSet),
                               RNA = lapply(alignment_file,Biostrings::readRNAStringSet),
                               AA = lapply(alignment_file,Biostrings::readAAStringSet))
  } else {
    my_seq_alignment <- alignment_list
  }
  k <- length(my_seq_alignment)
  if (is.null(group_dist)){
    group_result <- matrix(0,length(my_seq_alignment),length(my_seq_alignment))
    cat("\nCalculate the score between alignment results\n")
    pb1 <- txtProgressBar(style=3,title = "Calculate the score between alignment results")
    for (i in 1:nrow(group_result)){
      setTxtProgressBar(pb1, i/nrow(group_result),
                        title = "Calculate the score between alignment results")
      for (j in 1:ncol(group_result)){
        if (i <= j){
        alignment_result <- group_alignment(my_seq_alignment[[i]],my_seq_alignment[[j]],...)
        alignment_dims <- dim(alignment_result$matrix)
        alignment_score <- alignment_result$matrix[prod(dim(alignment_result$matrix))]
        alignment_length <- width(alignment_result$alignment)[1]
        group_result[i,j] <- alignment_score / alignment_length
        group_result[j,i] <- group_result[i,j]
        }
      }
    }
    close(pb1)
    cat("Calculate the dist between alignment results...")
    group_result <- scale(group_result)
    min_score <- min(group_result)
    diag(group_result) <- -10000
    max_score <- max(group_result)
    group_result <- (max_score - group_result) / (max_score - min_score)
    group_dist <- as.dist(group_result)
    cat("done\n")
  }
  group_hcl <- hclust(group_dist)
  group_tree <- cutree(group_hcl,1:(k-1))

  cat("Combine alignment results\n")
  pb2 <- txtProgressBar(style=3,label = "Combine alignment results")

  for (i in ncol(group_tree):1){
    setTxtProgressBar(pb2, (ncol(group_tree)-i+1)/ncol(group_tree),
                      label = "Combine alignment results")
    index <- which(duplicated(group_tree[,i]))
    group <- group_tree[index,i]
    combine <- which(group_tree[,i] == group)
    my_seq_alignment[[combine[1]]] <- group_alignment(my_seq_alignment[[combine[1]]],
                                                      my_seq_alignment[[combine[2]]],...)$alignment
    my_seq_alignment <- my_seq_alignment[-combine[2]]
    group_tree <- group_tree[-combine[2],]
  }
  close(pb2)
  if (tree){
    result <- list(combined_alignment = my_seq_alignment[[1]],tree = group_hcl)
    return(result)
  }
  return(my_seq_alignment[[1]])
}
#===============================================================================
##This function is in ggseqlogo package
makePFM <- function(seqs, seq_type='auto', namespace=NULL, keep_letter_mat=F){

  letter_mat = NA
  if(is.matrix(seqs)){
    # Process matrix
    if(is.null(rownames(seqs))) stop('Matrix must have letters for row names')

    num_pos = ncol(seqs)

    # Get namespace
    ns = findNamespace(rownames(seqs), seq_type, namespace)
    namespace = ns$namespace
    seq_type = ns$seq_type

    nseqs = NULL

    bg_prob = NA
    pfm_mat = seqs
    pfm_mat = apply(pfm_mat, 2, function(x) x / sum(x, na.rm=T))

    missing_rows = setdiff(namespace, rownames(pfm_mat))

    if(length(missing_rows) > 0){
      miss = matrix(rep(0, length(missing_rows) * ncol(pfm_mat)), nrow=length(missing_rows), dimnames = list(missing_rows))
      pfm_mat = rbind(pfm_mat, miss)
    }

    pfm_mat = pfm_mat[namespace,]

  }else{
    # Process sequences

    # Number of positions in alignment
    num_pos = nchar(seqs[1])
    # Number of sequences
    nseqs = length(seqs)
    # Letter matrix
    letter_mat = letterMatrix(seqs)


    # Get namespace
    ns = findNamespace(letter_mat, seq_type, namespace=namespace)
    namespace = ns$namespace
    seq_type = ns$seq_type

    # Construct PWM
    pfm_mat = apply(letter_mat, 2, function(pos.data){
      # Get frequencies
      t = table(pos.data)
      # Match to aa
      ind = match(namespace, names(t))
      # Create column
      col = t[ind]
      col[is.na(col)] = 0
      names(col) = namespace
      # Do relative frequencies
      col = col / sum(col)
      col
    })

    mat = matrix((letter_mat %in% namespace), nrow=nrow(letter_mat))
    attr(pfm_mat, 'nongapped') = apply(mat, 2, sum)
    attr(pfm_mat, 'nseqs') = nseqs
  }

  # Number of letters in ns
  N = length(namespace)

  # Assign seq type and namespace as attributes
  attr(pfm_mat, 'seq_type') = seq_type
  attr(pfm_mat, 'namespace') = namespace

  # Non-gapped columns
  if(seq_type == 'aa') namespace = c(namespace, 'X', 'B', 'Z')

  # Information content
  attr(pfm_mat, 'bits') = computeBits(pfm_mat, N, nseqs)

  # Assign AA names to rows/pos col
  rownames(pfm_mat) = namespace
  colnames(pfm_mat) = 1:num_pos

  if(keep_letter_mat) return(list(letter_mat = letter_mat, pfm=pfm_mat))

  return(pfm_mat)
}
#===============================================================================
##This function is in ggseqlogo package
letterMatrix <- function(input){
  # Ensure kmers are the same length characters
  seq.len = sapply(input, nchar)
  num_pos = seq.len[1]
  if(! all(seq.len == num_pos)) stop('Sequences in alignment must have identical lengths')

  # Construct matrix of letters
  split = unlist( sapply(input, function(seq){strsplit(seq, '')}) )

  t( matrix(split, seq.len, length(split)/num_pos) )
}
#===============================================================================
##This function is in ggseqlogo package
findNamespace <- function(letter_mat, seq_type, namespace){

  # Get all letters in our alignment
  sp = as.character(letter_mat)

  # Other namespace
  if(seq_type == "other"){
    if(is.null(namespace))
      stop('seq_type of "other" must have a defined namespace')

    namespace = as.character(namespace)
    # Get unique
    namespace = unique( unlist(strsplit(namespace, '')) )


    # Validate
    non_alphanumeric = grepl('[^a-zA-Z0-9\u03bb\u03b1\u03b2\u0393\u03b3\u0394\u03b4\u03b5\u03b6\u03b7\u03b8\u0398\u03b9\u03ba\u039b\u039b\u03bc\u039e\u03be\u03a0\u03c0\u03c1\u03c3\u03c4\u03c5\u03a6\u03c6\u03c7\u03c8\u03a8\u03a9\u03c9]', namespace)
    if( any( non_alphanumeric ) )
      stop('All letters in the namespace must be alphanumeric')

    # Ensure there is something in each column
    # apply(letter_mat, 2, function(column_letters){
    #   int = intersect(namespace, column_letters)
    #   if(length(int) == 0)
    #     stop('The alignment has no letters in namespace match aligned sequences in at least one column')
    # })

  }else{
    if(!is.null(namespace))
      stop('For custom namespaces please set seq_type to "other"')

    # Guess sequence type
    if(seq_type == "auto")
      seq_type = guessSeqType(sp)

    # Get predefined namespace
    namespace = get( sprintf('.%s_NAMESPACE', toupper(seq_type)) )()
  }

  return(list(seq_type = toupper(seq_type),
              namespace = namespace))
}

#===============================================================================
##This function is in ggseqlogo package
computeBits <- function(pwm, N=4, Nseqs=NULL){
  Nseqs = attr(pwm, 'nongapped')
  H_i = - apply(pwm, 2, function(col) sum(col * log2(col), na.rm=T))
  e_n = 0
  if(!is.null(Nseqs)) e_n = (1/logb(2)) * (N-1)/(2*Nseqs)

  R_i = log2(N) - (H_i  + e_n)
  # Set any negatives to 0
  R_i = pmax(R_i, 0)
  return(R_i)
}

#===============================================================================
