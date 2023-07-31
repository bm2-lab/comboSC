######################----filter cells--------#####################
cell_filt <- function(exp, metadata) {
    pbmc <- CreateSeuratObject(counts = exp, meta.data = metadata)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
    pbmc <- NormalizeData(object = pbmc,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)
    pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 &
                        nFeature_RNA < 25000 & percent.mt < 5)
    pbmc <- FindVariableFeatures(pbmc)
    pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
    pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 10)
    pbmc <- FindNeighbors(pbmc, reduction = "pca",nn.eps = 0.5)
    pbmc <- FindClusters(pbmc, n.start = 10)
    return(pbmc)
}

pca_cluster <- function(pbmc) {
    pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 10)
    pbmc <- FindNeighbors(pbmc, reduction = "pca",nn.eps = 0.5)
    pbmc <- FindClusters(pbmc, n.start = 10)
    pbmc <- RunUMAP(pbmc, dims = 1:60, nthreads = 4, max_iter = 2000)
    pbmc <- RunTSNE(pbmc, dims = 1:60, nthreads = 4, max_iter = 2000)
    return(pbmc)
}

###################----Remove batch effects----####################
rm_batch <- function(exp, metadata = bcc_metadata) {
    if ("batch" %in% colnames(metadata)) {
        exp <- align_cds(exp, alignment_group = "batch")
    }
    return(exp)    
}

######################----remove drop out----######################
rm_drop <- function(pbmc) {
    bcc_exped <- saver(pbmc@assays[[1]][], ncores = 8)
    return(bcc_exped)
}


###################################################################
######################---cluster and visual---#####################
###################################################################
Ana <- function(Pbmc){
    lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
    source("./Resources/sc-type-master/R/gene_sets_prepare.R")
    source("./Resources/sc-type-master/R/sctype_score_.R")
    db_ = "./Resources/sc-type-master/ScTypeDB_full.xlsx"
    tissue = "Immune system" 
    gs_list = gene_sets_prepare(db_, tissue)
    # pbmc <- readRDS("/home/tangchen/scRNA/RDS/bcc.rds")
    pbmc <- Pbmc
    es.max = sctype_score(scRNAseqData = pbmc[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 


    cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
    }))

    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    sctype_ident <- data.frame(sctype_scores[!duplicated(unlist(sctype_scores[,1])),c(1,2)])
    si_dict <- sctype_ident[,2]
    names(si_dict) <- sctype_ident[,1]
    pbmc@meta.data$cluster <- unlist(si_dict[pbmc@meta.data$seurat_cluster])
    return(pbmc)
}
###################################################################
########################--------main--------#######################
###################################################################
scrna_preprocess <- function(exp, metadata, batch = FALSE, core = 8) {
    pbmc <- cell_filt(exp, metadata)
    if(!("cluster" %in% colnames(metadata))){
        pbmc <- Ana(pbmc)
    }
    pbmc <- pca_cluster(pbmc)
    if(batch) {
        rm_batch(cds) 
    }
    return(pbmc)
}

