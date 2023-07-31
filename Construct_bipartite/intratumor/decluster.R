
################### ---tumor detail cluster----######################
tumor_scale <- function(pbmc) {
#' Typing cancer cells
    tumor_col <- c("Tumor_1", "Tumor_2", "Tumor", "Malignant")
    tumor_pbmc <- subset(pbmc, subset = cluster %in% tumor_col)
    tumor_pbmc <- FindVariableFeatures(tumor_pbmc)    
    tumor_pbmc <- RunPCA(tumor_pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 10)
    tumor_pbmc <- FindNeighbors(tumor_pbmc, reduction = "pca", dims = 1:70, nn.eps = 0.5)
    return(tumor_pbmc)
}



###### dimensional reductio and clustering
dimrandcluster <- function(pbmc) {
#' dimensional reductio and clustering, if the result is not ideal,
#' you can try to change the parameters inside
    pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 10)
    pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:70, nn.eps = 0.5)
    pbmc <- FindClusters(pbmc, n.start = 10, resolution = 0.8)
    return(pbmc)
}


###################################################################
############# ---blast with reference cell line---#################
###################################################################












## Reasonable features can make the similarity algorithm match better

add_modu <- function(module_gene, matrix, rfgene) {
    #' improve the weight of the import gene.
    local_gene <- Reduce(
        intersect,
        list(
            v1 = module_gene,
            v2 = rownames(rfgene),
            v3 = rownames(matrix)
        )
    )
    imp_genes_mat1 <- matrix[rownames(matrix) %in% local_gene, ]
    imp_genes_mat2 <- rfgene[rownames(rfgene) %in% local_gene, ]
    imp_genes_mat1 <- imp_genes_mat1[order(row.names(imp_genes_mat1)), ]
    imp_genes_mat2 <- imp_genes_mat2[order(row.names(imp_genes_mat2)), ]
    ls1 <- rbind(imp_genes_mat1, imp_genes_mat1)
    ls2 <- rbind(imp_genes_mat2, imp_genes_mat2)
    rownames(ls1) <- c(1:nrow(ls1))
    rownames(ls2) <- c(1:nrow(ls2))
    new_ma <- rbind(matrix, ls1)
    new_rf <- rbind(rfgene, ls2)
    return(list(new_ma, new_rf))
}

only_modulegene <- function(expression.profile, module_gene) {
    #' Only select genes associated with cancer (
    #' exp:only_modulegene(me,essential_genes)
    deducted <- expression.profile[
        rownames(expression.profile) %in% module_gene,
    ]
    return(deducted)
}

mapfrom_Addmodule <- function(n, rfgene) {
    #' append rf gene in the rfgene
    new.pro <- add_modu(module_gene, pbmc.matrix, rfgene)
    pbmc.matrix <- new.pro[[1]]
    mj <- PAIT(n)
    sce <- define_sce(mj[[1]], mj[[2]])
    map(new.pro[[2]])
}

Sum_row <- function(df, n) {
    #' sum df to exract
    dis_row <- c()
    for (i in c(2:nrow(df))) {
        sum_row <- sum(df[i, n:ncol(df)])
        if (sum_row == 0) {
            dis_row <- c(dis_row, i)
        }
    }
    dis_row <- unlist(dis_row)
    del_len <- length(dis_row)
    for (j in c(1:del_len - 1)) {
        de <- dis_row[del_len - j]
        df <- df[-de, ]
    }
    muti <- list(df, dis_row)
    return(df)
}

define_sce <- function(df, coldata = TRUE) {

    #' define our database
    if (coldata == TRUE) {
        coldata <- data.frame(cell = colnames(df), cell_type1 = colnames(df))
    }
    e <- SingleCellExperiment(
        assays = list(normcounts = as.matrix(df)),
        colData = coldata
    )
    logcounts(e) <- log2(normcounts(e) + 1)
    rowData(e)$feature_symbol <- rownames(e)
    # isSpike(e, "ERCC") <- grepl("^ERCC-", rownames(e))
    is.spike <- grepl("^ERCC", rownames(e))
    e <- splitAltExps(e, ifelse(is.spike, "ERCC", "gene"))
    e <- e[!duplicated(rownames(e)), ]
    return(e)
}

cell_coverage <- function(pbmc, scmapcluster_results2) {
#' According to scmapcluster_results2, get how many cells
#  are matched(not cluster);now. cluster also
#' Example:cell_matched <- cell_coverage(pbmc, scmapcluster_results2)
    matched_cluster <- unique(scmapcluster_results2$scmap_cluster_labs) %>%
        as.vector() %>%
        .[-which(. == "unassigned")]
    cluster_num <- table(pbmc$seurat_clusters) %>% as.data.frame
    cell_matched <- sapply(matched_cluster, function(x) {
        cluster_num[which(cluster_num[, 1] == x), 2]
    }) %>% unlist %>% sum
    cell_coverage <- cell_matched / sum(cluster_num$Freq)
    cluster_coverage <- (length(df_map) - 1) / length(unique(identiy[ ,1]))
    general_coverage <- (cell_coverage + cluster_coverage) / 2
    return(general_coverage)
}

map <- function(sce = sce, rfgene = rfgene_rpkm, threshold = 0.5) {
    #' begin scmap,threshold belong to similarity
    #' exp: map(rfgene_count,0.5)
    e <- define_sce(rfgene)
    sce <- selectFeatures(sce)
    # if need pic ,please add , suppress_plot = FALSE;
    sce <- indexCluster(sce)
    # if heatmap is needed
    # MAY:heatmap(as.matrix(metadata(sce)$scmap_cluster_index))
    scmapcluster_results2 <<- scmapCluster(
        projection = e,
        index_list = list(yan = metadata(sce)$scmap_cluster_index),
        threshold = threshold
    )
    df_map <<- table(scmapcluster_results2$scmap_cluster_labs)
    cell_matched <- cell_coverage(pbmc, scmapcluster_results2)
    return(cell_matched)
}

#map(sce, ccleandgdsc)


############################################################
###################### -----opt_cluster------###############
############################################################

res_cov <- function(pbmc, restocluster,res, threshold, rfgene = rfgene_rpkm) {
    #' Under a certain resolution get the Coverage, print specific match case
    #' exp:res_cov(0.5)
    identiy <<- as.data.frame(unlist(restocluster[as.character(res)]))
    colnames(identiy) <- "cell_type1"
    identiy["cell"] <- rownames(identiy)
    sce <<- define_sce(only_modulegene(pbmc@assays[[1]][],
                                    essential_genes), identiy)
    result_map <- map(sce, only_modulegene(rfgene, essential_genes),threshold = threshold)

    return(result_map)
}

res_opt <- function(res_rank, threshold = 0.5, rfgene = rfgene_rpkm, pbmc. = pbmc) {
    #' Select the most appropriate resolution in the resolution of an interval
    #' exp:llss=res_opt(res_rank,0.3)

    cl <- makeCluster(4, type="FORK")
    #registerDoParallel(cl) 
    re <- function(res) as.vector(FindClusters(pbmc, n.start = 10, resolution = res)@meta.data[paste("RNA_snn_res.",as.character(res),sep = "")])
    restocluster <- foreach(1, .combine="cbind", .packages = c("Seurat")) %dopar% re(res_rank)

    names(restocluster) <- res_rank
    coverage <- sapply(res_rank, res_cov, restocluster = restocluster, threshold = threshold,pbmc = pbmc)
    stopCluster(cl)
    df <- data.frame(resolution = res_rank, coverage = coverage)

    resest <- res_rank[which.max(coverage)]
    coverest <- res_cov(pbmc, resest, restocluster = restocluster,threshold = threshold)



    return(resest)
}



main_optcluster <- function(res_rank = seq(0.4, 5, 0.2),
                            threshold = 0.5,
                            rfgene = rfgene_rpkm) {
    best_res <- res_opt(res_rank, threshold, rfgene)
    return(best_res)
}

###################################################################
###################### -------matrix--------######################
###################################################################

map_topn <- function(scmapcluster_results, num, rfgene = rfgene_rpkm) {
    #' print the most matched cell line of the Standard cell lines
    #' tips:most_matched(scmapcluster_results,3)
    x <- as.vector(unique(scmapcluster_results$scmap_cluster_labs))
    x <- x[which(x != "unassigned")]
    mml <- c()
    for (i in x) {
        exist <- scmapcluster_results$scmap_cluster_labs == i
        mx <- sort(scmapcluster_results$scmap_cluster_siml[exist],
                    decreasing = TRUE)[1:min(num, sum(exist == TRUE))]
        a <- unlist(lapply(mx, function(x) {
            colnames(rfgene)[
                which(scmapcluster_results$scmap_cluster_siml == x)]
        }))
        mml <- c(mml, a)
    }
    mml <- as.data.frame(array(mml, dim = c(num, length(x))))
    colnames(mml) <- x
    return(mml)
}

creat_use_cell <- function(mel_top, rfdrug, quan = "LN_IC50") {
    #' define target ,creat useful celllines and drug
    group <- list()
    col_n <- 1
    while (col_n <= ncol(mel_top)) {
        row_n <- 1
        while (row_n <= nrow(mel_top)) {
            a <- mel_top[row_n, col_n]
            if (is.na(a)) {
                row_n <- row_n + 1
                break
            }
            n <- str_extract_all(a, "[0-9]")[[1]][1]
            le <- substr(str_extract_all(a, "[A-Z]")[[1]][1], 1, 2)
            if (is.na(n)) {
                st1 <- rfdrug
            } else {
                st1 <- rfdrug[grep(n, rfdrug$CELL_LINE_NAME,
                    ignore.case = TRUE
                ), ]
            }
            if (length(le) > 0) {
                st1 <- st1[grep(le, st1$CELL_LINE_NAME, ignore.case = TRUE), ]
            } # prevent null to useless
            st1 <- filt_drug(st1, 0, quan = quan) # do nothing to z-score
            if (length(st1[, 1]) == 0) {
                row_n <- row_n + 1
            } else {
                eval(parse(text = paste("group$", a, "=st1", sep = "")))
                break
            }
        }
        col_n <- col_n + 1
    }
    return(group)
}



pro_drug_df <- function(mel_top, rfdrug, not_write = TRUE, quan = "LN_IC50") {
#' construct the connection matrixs,mel_top is the topn matrix of the cell lines
    group <- creat_use_cell(mel_top, rfdrug, quan = quan)
    drug_list <- unique(unlist(lapply(group, rownames)))
    drug_df <- data.frame()
    for (i in c(1:length(group))) {
        for (j in drug_list) {
            if (j %in% rownames(group[[i]])) {
                drug_df[j, names(group)[i]] <-
                    group[[i]][rownames(group[[i]]) == j, 2]
            } else {
                drug_df[j, names(group[i])] <- 0
            }
        }
    }
    if (not_write != TRUE) {
        write.csv(drug_df, not_write, row.names = TRUE)
    }
    return(drug_df)
}

###################################################################
########################--------main--------#######################
###################################################################
###################################################################

######################-----------------######################
decluster <- function(pbmc,
                    res_rank = seq(0.4, 3, 0.2),
                    threshold = threshold,
                    rfgene = rfgene_rpkm,
                    essential_genes = essential_genes) {
    cat('###################################################\nStep3.1 is start!\n#################################################\n')
    pbmc <- tumor_scale(pbmc)
    pbmc <<- FindClusters(pbmc, n.start = 10)
    cat('###################################################\nStep3.1 is OK!\n#################################################\n')
    identiy <- as.data.frame(pbmc$seurat_clusters)
    colnames(identiy) <- "cell_type1"
    identiy["cell"] <- rownames(identiy)
    sce <<- define_sce(only_modulegene(pbmc@assays[[1]][], essential_genes),identiy)
    cat('###################################################\nStep3.2 is OK!\n#################################################\n')
    #best_res <- main_optcluster(res_rank = res_rank, threshold = threshold)
    best_res <- res_opt(pbmc = pbmc, res_rank = res_rank, threshold = threshold)
    cat('###################################################\nStep3.3 is OK!\n#################################################\n')
    mel_top <<- map_topn(scmapcluster_results2, 10, rfgene_rpkm)
    drug_df <- pro_drug_df(mel_top, rfdrug, "./Resources/bio_matrix.csv", quan = "LN_IC50")
    return(drug_df)
}



