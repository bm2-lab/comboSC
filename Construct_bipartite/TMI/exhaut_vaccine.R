###################################################################
######## ------creat dataset and Initial clustering------##########
###################################################################

get_earliest_principal_node <- function(cds,
                                        colname = "seurat_clusters",
                                        time_bin = "CD8_ex_T_cells") {
#' set the root of trajectory
    cell_ids <- which(colData(cds)[, colname] == time_bin)
    closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
        (which.max(table(closest_vertex[cell_ids, ]))))]
    root_pr_nodes
}

pre_cluster_learn_order <- function(cds_subset, reduction_method = "UMAP",
                            cores = 32, num_dim = 100,
                            orderp = c("seurat_clusters", "CD8_ex_T_cells")) {
    #' intial operate on the cds
    cds_subset <- preprocess_cds(cds_subset, num_dim = 100)
    cat("Step 1: Normalize and pre-process the data  ^_^\n")
    cds_subset <- reduce_dimension(cds_subset,
        reduction_method = reduction_method, cores = cores
    )
    cat(" Step 3: Reduce the dimensions using UMAP ^_^ \n")
    cds_subset <- cluster_cells(cds_subset, reduction_method = reduction_method)
    cat(" Step 4: Cluster the cells ^_^ \n")
    cds_subset <- learn_graph(cds_subset)
    cat(" Step 5: Learn a graph ^_^ \n")
    ## Step 6: order cells
    if (length(orderp) == 2) {
        cds_subset <- order_cells(cds_subset,
            root_pr_nodes = get_earliest_principal_node(
                cds_subset,
                colname = orderp[1], time_bin = orderp[2]
            )
        )
    }
    cat(" Step 6: order cells ^_^ \n")
    return(cds_subset)
}


###################################################################
###################### -----find the path-----#####################
###################################################################
state_pseudo <- function(cds, reduction = "UMAP") {
#' coldata add the pseudotime and plot.
    cds$Pseudotime <- pseudotime(cds, reduction_method = reduction) %>% as.vector()
    cds_mete <- colData(cds) %>% data.frame %>% filter(., Pseudotime != Inf)
    CairoPDF(file = "./Resources/box_seutime20.32.06.pdf")
    ggplot(cds_mete, aes(x = seurat_clusters, y = Pseudotime,
                        fill = seurat_clusters), position = "dodge") +
                        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8))
    dev.off()
    return(cds_mete)
}

######################-----------------######################
exhau_interval <- function(cds_mete, group2, group1 = "CD8_ex_T_cells") {
    p1 <- cds_mete[cds_mete["seurat_clusters"] == group1, 5]
    p2 <- cds_mete[which(cds_mete["seurat_clusters"] == group2), 5]
    interval1 <- quantile(p1, c(0.25, 0.5, 0.75)) %>% unlist %>% as.vector
    interval2 <- quantile(p2, c(0.25, 0.5, 0.75)) %>% unlist %>% as.vector
    o <- order(c(interval1[1], interval2[1]))
    v1 <- list(interval1, interval2)[[o[2]]]
    v2 <- list(interval1, interval2)[[o[1]]]
    # p_order <- rank(c(interval1[[1]], interval1[[3]],
    #                 interval2[[1]], interval2[[3]]))
    if (v1[3] < v2[3]) {
        pseudotime1 <- v1[2:3]
        pseudotime2 <- v2[1:2]
    } else {
        pseudotime1 <- v2[c(1,3)]
        pseudotime2 <- v2[c(1,3)]
    }
    if (o[1] == 2) {
        ls <- pseudotime1
        pseudotime1 <- pseudotime2
        pseudotime2 <- ls
    }
    return(list(pseudotime1, pseudotime2))
}


######################-----------------######################

exhaust_cell <- function(cds_subset, cds_mete,
                    pseudotime, group, colname = "seurat_clusters") {

    cellid <- cds_mete[which(cds_mete$Pseudotime > pseudotime[1] &
                                cds_mete$Pseudotime < pseudotime[2] &
                                cds_mete[colname] == group), 1]
    if(length(cellid) < 50) {
        cellid <- cds_mete[which(cds_mete[colname] == group), 1]
    }
    return(cellid)
}

exhaust_path <- function(cds_subset, cds_mete,
                        pseudotime1, pseudotime2,
                        group2, group1 = "CD8_ex_T_cells",
                        colname = "seurat_clusters") {
    exhaust_cellid <- exhaust_cell(cds_subset, cds_mete, pseudotime1, group1)
    preexhaust_cellid <- exhaust_cell(cds_subset, cds_mete, pseudotime2, group2)
    
    act_exp <- assay(cds_subset)[, colnames(assay(cds_subset)) %in%
        union(exhaust_cellid, preexhaust_cellid)]
    indx <- factor(colnames(act_exp),
        levels = union(exhaust_cellid, preexhaust_cellid)
    )
    act_exp <- act_exp[, order(indx)] # 换到指定顺序
    return(list(act_exp, exhaust_cellid, preexhaust_cellid))
}

###################################################################
###################### ----de analysis--------#####################
###################################################################

de_analysis <- function(list_to_de, workers = 32) {
#' Differential expression analysis of cells on the way to exhaustion
    param <- MulticoreParam(workers = 6, progressbar = TRUE)
    register(param)
    if(dim(list_to_de[[1]])[2] < 50) {
        cat("warnning! Too few cells at either end of the trajectory！")
    }
    group <- factor(c(
        rep(1, length(list_to_de[[2]])),
        rep(2, length(list_to_de[[3]]))
    ))
    results <- DEsingle(
        counts = list_to_de[[1]] %>% as.matrix(), group = group,
        parallel = TRUE, BPPARAM = param
    )
    results_classified <- DEtype(
        results = results,
        threshold = 0.01
    )
    results_sig <- results_classified[which(results_classified$pvalue < 1e-2 &
                                    results_classified$norm_foldChange != "Inf"),
                                    ]
    de_gene <- data.frame("Gene" = results_sig %>% rownames(), 
                        foldchange = results_sig$norm_foldChange)
    return(de_gene)
}

de_analysis_Srurat <- function(list_to_de, p_val_adjust = FALSE) {
#' Differential expression analysis of cells on the way to exhaustion
    if(dim(list_to_de[[1]])[2] < 50) {
        cat("warnning! Too few cells at either end of the trajectory！")
    }
    pbmc <- CreateSeuratObject(counts = list_to_de[[1]])
    pbmc <- FindVariableFeatures(pbmc)   
    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 10)
    pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:70, nn.eps = 0.5)
    pbmc <- FindClusters(pbmc, n.start = 10, resolution = 2)

    groupA <- data.frame("names" = list_to_de[[2]], "group" = "A")
    groupB <- data.frame("names" = list_to_de[[3]], "group" = "B")
    group_meta <- rbind(groupA,groupB)
    g <- group_meta[,2] %>% as.vector
    names(g) <- unlist(group_meta[,1])
    pbmc$seurat_clusters <- g[unlist(rownames(pbmc@meta.data))]
    pbmc@active.ident <- as.factor(pbmc@meta.data$seurat_clusters)
    names(pbmc@active.ident) <- rownames(pbmc@meta.data)
    markers <- FindAllMarkers(object = pbmc, 
                            only.pos = TRUE,
                            logfc.threshold = 0.25)
    if (p_val_adjust) {
        markers <- markers[markers$p_val_adj < 0.05, ] 
        markers[which(markers$cluster == "B"),"avg_logFC"] <- -markers[which(markers$cluster == "B"),"avg_logFC"]
        markers$logFC <- exp(markers$avg_logFC)
        de_gene <- markers[c("gene","logFC")]
    } else{
        markers[which(markers$cluster == "B"),"avg_logFC"] <- -markers[which(markers$cluster == "B"),"avg_logFC"]
        markers$logFC <- exp(markers$avg_logFC)
        de_gene <- markers[c("gene","logFC")]
    }
    return(de_gene)
}

###################################################################
######################---------drug_match--------##################
###################################################################    
biograph_python <- function(men_delist, act_delist) {
    men_de_gene <<- de_analysis(men_delist) 
    write.csv(men_de_gene, "./Resources/men_de_gene.csv",
            row.names = FALSE, quote = FALSE)
    command <- "./Construct_bipartite/TMI/drug_match.py ./Resources/men_de_gene.csv ./Auxiliary/men_linshi" 
    system(command)
    df1 <<- fread("./Auxiliary/men_linshi",data.table = FALSE)
    df1 <- na.omit(df1)
    act_de_gene <<- de_analysis(act_delist) 
    write.csv(act_de_gene, "./Resources/act_de_gene.csv",
            row.names = FALSE, quote = FALSE)
    command <- "./Construct_bipartite/TMI/drug_match.py ./Resources/act_de_gene.csv ./Auxiliary/act_linshi" 
    system(command)
    df2 <<- fread("./Auxiliary/act_linshi",data.table = FALSE)
    df2 <- na.omit(df2)
    system("rm ./Auxiliary/*linshi")
    bio_matrix <- full_join(df1, df2, by = "pubchem_ID")
    bio_matrix <<- bio_matrix[, c(2, 3, 5)]
    bio_matrix[is.na(bio_matrix)] <- 0
    write.csv(bio_matrix, "./Resources/bio_matrix.csv",
            row.names = FALSE, quote = FALSE) 
    return(bio_matrix)
}

dfgene_drug <- function(df_gene) {
  gene_input <- as.data.frame(df_gene,stringsAsFactors = FALSE)
  gene_input[,2] <- as.numeric(gene_input[,2])
  gene_input[,1] <- unlist(gene_input[,1]) %>% as.character()
  gene_input <- na.omit(gene_input)
  gene_input[gene_input[,2] < 0.5,2] = "-1"
  gene_input[gene_input[,2] > 2,2] = "1"
  
  gene_input <- gene_input[which(gene_input[,2] %in% c("1","-1")), ]
  ESMID <- mapIds(org.Hs.eg.db, keys = gene_input[,1], keytype = "SYMBOL", column="ENSEMBL")
  ESMID <- na.omit(ESMID)
  gene_input <- gene_input[gene_input[,1] %in% names(ESMID),]
  rownames(gene_input) <- paste(ESMID[gene_input[,1]],"_at",sep = "")
  colnames(gene_input) <- c("feature","direction")

  res <- apply(drug.perturbation[,,c("tstat", "fdr")],
            2, function(x, HDAC){ 
            return(connectivityScore(x=x,
                                    y=HDAC[,2,drop=FALSE],
                                        method="fgsea", nperm=100))
            }, HDAC=gene_input)
  
  rownames(res) <- c("Connectivity", "P Value")

  rest <- t(res) %>% as.data.frame()
  rest <- rest[rest[,2]<0.05,]
  rest <- rest[order(rest[,1], decreasing=FALSE),]
  rest$drug <- rownames(rest)
  return(rest)
}

Pharm <- function(men_delist, act_delist) {
    men_de_gene <<- de_analysis(men_delist) 
    act_de_gene <<- de_analysis(act_delist) 
    act_drug <- dfgene_drug(act_de_gene)
    men_drug <- dfgene_drug(men_de_gene)
    drug_df <- full_join(act_drug,men_drug,by="drug")
    drug_df[is.na(drug_df)] <- 0
    write.csv(drug_df, "./Resources/bio_matrix.csv",
            row.names = FALSE, quote = FALSE) 
    return(drug_df)
}

dftam_drug <- function(df_gene) { 
  gene_input <- as.data.frame(df_gene,stringsAsFactors = FALSE)
  gene_input[,2] <- as.numeric(gene_input[,2])
  gene_input[,1] <- unlist(gene_input[,1]) %>% as.character()
  gene_input <- na.omit(gene_input)
  gene_input[gene_input[,2] < 0,2] = "-1"
  gene_input[gene_input[,2] > 0,2] = "1"
  
  gene_input <- gene_input[which(gene_input[,2] %in% c("1","-1")), ]
  ESMID <- mapIds(org.Hs.eg.db, keys = gene_input[,1], keytype = "SYMBOL", column="ENSEMBL")
  ESMID <- na.omit(ESMID)
  gene_input <- gene_input[gene_input[,1] %in% names(ESMID),]
  rownames(gene_input) <- paste(ESMID[gene_input[,1]],"_at",sep = "")

  res <- apply(drug.perturbation[,,c("tstat", "fdr")],
            2, function(x, HDAC){ 
            return(connectivityScore(x=x,
                                    y=HDAC[,2,drop=FALSE],
                                        method="fgsea", nperm=100))
            }, HDAC=gene_input)
  
  rownames(res) <- c("Connectivity", "P Value")

  rest <- t(res) %>% as.data.frame()
  rest <- rest[rest[,2]<0.05,]
  rest <- rest[order(rest[,1], decreasing=FALSE),]
  rest$drugs <- rownames(rest)
  return(rest)
}

Seu_tra <- function(pbmc){
    tra_path <- list(c("M1","M2"),c("CD8_act_T_cells","CD8_ex_T_cells"),c("CD8_mem_T_cells", "CD8_ex_T_cells"),c("mCAF","CAF"))
    pbmc_tra_path <- c(0,0,0,0)
    for(i in c(1:length(tra_path))){
        if(sum(tra_path[[i]] %in% unique(pbmc@meta.data$cluster)) == 2){
            pbmc_tra_path[i] <- 1    
        }
    }
    if(sum(pbmc_tra_path) == 0){
        print("Drug recommendation cannot be made because of the absence of immune microenvironment cells.")
    }
    if(sum(pbmc_tra_path) == 1){
        tra <- tra_path[pbmc_tra_path == 1]
        drug_df <- Drug_rep(pbmc,tra[[1]][1],tra[[1]][2])    
    }

    if(sum(pbmc_tra_path) == 2){
        tra <- tra_path[pbmc_tra_path == 1]
        drug_d1 <- Drug_rep(pbmc,tra[[1]][1],tra[[1]][2])    
        drug_d2 <- Drug_rep(pbmc,tra[[2]][1],tra[[2]][2]) 
        drug_df <- full_join(drug_d1,drug_d2,by="drugs")
        drug_df[is.na(drug_df)] <- 0
        drug_df$Score <- drug_df$Connectivity.x + drug_df$Connectivity.y
        drug_df <- drug_df[order(drug_df$Score, decreasing=FALSE),]
    }
    return(drug_df)
}

Drug_rep <- function(pbmc,group1,group2) {
    actexh.markers <- FindMarkers(pbmc, ident.1 = group2, ident.2 = group1,group.by = 'cluster',logfc.threshold = 0.75,min.pct = 0.25)
    actexh.markers <- actexh.markers[actexh.markers$p_val_adj<0.05,]
    actexh.markers %>% dim
    markers <- actexh.markers
    TAMgene <- data.frame(feature = rownames(markers),direction =  markers[,2],row.names = rownames(markers),stringsAsFactors = FALSE)
    tam_drug <- dftam_drug(TAMgene)
}


####################### main #######################

exhaust_caccine <- function(pbmc,
                            ciliated_genes = essential_genes,
                            group2 = "CD8_act_T_cells",
                            group1 = "CD8_ex_T_cells") {
    drug.perturbation <<- readRDS("./Resources/CMAP_signatures.rds")
    immu_pbmc <- subset(pbmc, subset = cluster %in% ciliated_cluster)    
    if(ncol(immu_pbmc) < 4500){
        return(Seu_tra(pbmc))
    }
    counts.data <- immu_pbmc@assays[[1]][]
    gene_metadata <- data.frame(
        "Gene" = rownames(counts.data),
        "gene_short_name" = rownames(counts.data), row.names = rownames(counts.data)
    )
    cell_metedata <- data.frame(
        "Barcode" = colnames(counts.data), "Sample" = rep("bcc", colnames(counts.data) %>% length()),
        row.names = colnames(counts.data)
    )

    cell_metedata["seurat_clusters"] <- immu_pbmc@meta.data$cluster
    immue_cds <- new_cell_data_set(counts.data,
                            cell_metadata = cell_metedata,
                            gene_metadata = gene_metadata
    )
    immue_cds <- pre_cluster_learn_order(immue_cds)
    cds_mete <<- state_pseudo(immue_cds)
    act_interval <- exhau_interval(cds_mete, group2 = "CD8_act_T_cells", group1 = "CD8_ex_T_cells")
    men_interval <- exhau_interval(cds_mete, group2 = "CD8_mem_T_cells", group1 = "CD8_ex_T_cells")
    men_delist <- exhaust_path(   
        immue_cds, cds_mete,
        men_interval[[1]], men_interval[[2]], "CD8_mem_T_cells"
    )
    act_delist <- exhaust_path(
        immue_cds, cds_mete,
        act_interval[[1]], act_interval[[2]], "CD8_act_T_cells"
    )
    # act_interval, men_interval
    bio_matrix <- Pharm(men_delist, act_delist)
  return(bio_matrix)
}