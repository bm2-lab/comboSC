setwd("/home/tangchen/code/precision_medicine/Sctc")
source("./resources/dep_data.R")
source("./Data_preprocessing/scprocess.R")
source("./immunity_evaluation/immuscore.R")
source("./intelligent_optimization/graph_optimization.R")
source("./Construct_bipartite/intratumor/decluster.R")
source("./Construct_bipartite/TMI/exhaut_vaccine.R")
source("./Construct_bipartite/construct_bipartite.R")

ciliated_cluster <- "Macrophages"

load("/home/tangchen/vde/ls.cafmar.Rdata", envir = parent.frame(), verbose = FALSE)
###################################################################
######################---------CAFS--------######################
##################################################################
ciliated_cluster <- "CAFs"
TMI_metadata <- bcc_metadata[bcc_metadata$cluster %in% ciliated_cluster, ]

TMI_exp <- bcc_exp[, colnames(bcc_exp) %in%
                        unlist(as.vector(TMI_metadata[, 1]))]
indx <- factor(TMI_metadata$cell.id, levels = colnames(TMI_exp))
TMI_metadata <- TMI_metadata[order(indx), ]
gene_metadata <- data.frame(
    "Gene" = rownames(TMI_exp),
    "gene_short_name" = rownames(TMI_exp), row.names = rownames(TMI_exp)
) # gene names - only have hgnc_symbol
cell_metedata <- data.frame(
    "Barcode" = colnames(TMI_exp), 
    "Sample" = rep("bcc", TMI_metadata[, 1] %>% unlist %>% length), row.names = colnames(TMI_exp)
)
TMI_cds <- new_cell_data_set(TMI_exp %>% as.matrix(),
                            cell_metadata = cell_metedata,
                            gene_metadata = gene_metadata
)

TMI_cds <- pre_cluster_learn_order(TMI_cds,orderp = c("seurat_clusters", ciliated_cluster))

CairoPDF(file = "/home/tangchen/public/png/tmi.6.9.pdf")
plot_cells(TMI_cds,
    color_cells_by = "pseudotime",
    label_cell_groups = TRUE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5
)
dev.off()

CairoPDF(file = "/home/tangchen/public/png/tmi.6.90.pdf")
plot_cells(TMI_cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5
)
dev.off()

#patient005最多的CAF
cds_subset <- TMI_cds
## Step 1: Normalize and pre-process the data
cds_subset <- preprocess_cds(cds_subset, num_dim = 100)
## Step 3: Reduce the dimensions using UMAP
cds_subset <- reduce_dimension(cds_subset, reduction_method = "UMAP", cores = 32)
## Step 4: Cluster the cells
cds_subset <- cluster_cells(cds_subset, reduction_method = "UMAP", num_iter = 100, resolution=1e-10)
## Step 5: Learn a graph
cds_subset <- learn_graph(cds_subset)
## Step 6: order cells
orderp = c("seurat_clusters", "CAFs")
# cds_subset <- order_cells(cds_subset,
#             root_pr_nodes = get_earliest_principal_node(
#                 cds_subset,
#                 colname = orderp[1], time_bin = orderp[2]
#             )
# )

cds_subset <- order_cells(cds_subset)

CairoPDF(file = "/home/tangchen/public/png/3ls.pdf")
plot_cells(cds_subset, color_cells_by="cluster",
                label_cell_groups = TRUE,
                label_leaves = FALSE,
                label_branch_points = FALSE,
                graph_label_size = 1.5)
dev.off()

marker_test_res <- top_markers(cds_subset, group_cells_by="cluster", 
                        reference_cells=1000, cores=64)

marker_test_res[grep("ID", marker_test_res$gene_id),c(1,10)]

vCAF  <- c("Esam","Gng11","Higd1b","Cox4i2","Cygb","Gja4","Eng")
mCAF  <- c("Dcn","Col12a1","Mmp2","Lum","Mrc2","Bicc1","Lrrc15","Mfap5","Col3A1","Mmp14","Spon1","Pdgfrl","Serpinf1","Lrp1","Gfpt2","Ctsk","Cdh11","Itgbl1","Col6a2","Postn","Ccdc80","Lox","Vcan","Col1a1","Fbn1","Col1a2","Pdpn","Col6a1","Fstl1","Col5a2","Aebp1")
dCAF  <- c("Tspan2","Reck")

caf <- dCAF
for( g in caf){
    print(marker_test_res[grep("g", marker_test_res$gene_id), ])
}


###################################################################
######################---------MAR--------######################
###################################################################

ciliated_cluster <- "Macrophages"
TMI_metadata <- bcc_metadata[bcc_metadata$cluster %in% ciliated_cluster, ]

TMI_exp <- bcc_exp[, colnames(bcc_exp) %in%
                        unlist(as.vector(TMI_metadata[, 1]))]
indx <- factor(TMI_metadata$cell.id, levels = colnames(TMI_exp))
TMI_metadata <- TMI_metadata[order(indx), ]
gene_metadata <- data.frame(
    "Gene" = rownames(TMI_exp),
    "gene_short_name" = rownames(TMI_exp), row.names = rownames(TMI_exp)
) # gene names - only have hgnc_symbol
cell_metedata <- data.frame(
    "Barcode" = colnames(TMI_exp), 
    "Sample" = rep("bcc", TMI_metadata[, 1] %>% unlist %>% length), row.names = colnames(TMI_exp)
)
cell_metedata$seurat_clusters <- TMI_metadata$cluster
TMI_cds <- new_cell_data_set(TMI_exp %>% as.matrix(),
                            cell_metadata = cell_metedata,
                            gene_metadata = gene_metadata
)

orderp = c("seurat_clusters", "Macrophages")
#patient005最多的CAF
cds_subset <- TMI_cds
## Step 1: Normalize and pre-process the data
cds_subset <- preprocess_cds(cds_subset, num_dim = 100)
## Step 3: Reduce the dimensions using UMAP
cds_subset <- reduce_dimension(cds_subset, reduction_method = "UMAP", cores = 32)
## Step 4: Cluster the cells
cds_subset <- cluster_cells(cds_subset, reduction_method = "UMAP", num_iter = 100, resolution=1e-7)
table(as.character(clusters(cds_subset)[colnames(cds_subset)]))
## Step 5: Learn a graph
cds_subset <- learn_graph(cds_subset)
## Step 6: order cells
cds_subset <- order_cells(cds_subset,
            root_pr_nodes = get_earliest_principal_node(
                cds_subset,
                colname = orderp[1], time_bin = orderp[2]
            )
)

CairoPDF(file = "/home/tangchen/public/png/1ls.pdf")
plot_cells(cds_subset, color_cells_by="cluster",
                label_cell_groups = TRUE,
                label_leaves = FALSE,
                label_branch_points = FALSE,
                graph_label_size = 1.5)
dev.off()

marker_test_res <- top_markers(cds_subset, group_cells_by="cluster", 
                        reference_cells=1000, cores=64)


###################################################################
######################---------lung--------######################
###################################################################

mt <- fread("/home/tangchen/scRNA/lung/lungGSE127465/GSE127465_human_cell_metadata_54773x25.tsv") 


###################################################################
######################---------scCATCH--------#####################
###################################################################

pbmc <- CreateSeuratObject(as.matrix(pt1_exp), meta.data = pt1_metadata)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc <- NormalizeData(
    object = pbmc,
    normalization.method = "LogNormalize", scale.factor = 10000
)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
pbmc <- RunPCA(pbmc, npcs = 100, ndims.print = 1:5, nfeatures.print = 10)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:70, nn.eps = 0.5)
pbmc <- FindClusters(pbmc, n.start = 10, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:60, nthreads = 4, max_iter = 2000)


clu_markers <- findmarkergenes(object = pbmc, species = 'Human', cluster = 'All', match_CellMatch = FALSE,cancer = "Lung",tissue = "Skin-related",cell_min_pct = 0.25,logfc = 0.25, pvalue = 0.05)

clu_ann <- scCATCH(object = clu_markers$clu_markers, species = 'Human', cancer = "Lung", tissue = 'Skin-related')


###################################################################
######################-------garneet----------#####################
###################################################################
TMI_classifier <- train_cell_classifier(cds = cds_subset,
                                        marker_file = "/home/tangchen/public/cell_mark/ls.txt",
                                        db=org.Ce.eg.db::org.Ce.eg.db,
                                        cds_gene_id_type = "ENSEMBL",
                                        num_unknown = 10,   
                                        marker_file_gene_id_type = "ENSEMBL",
                                        cores = 32)

pbmc_cds <- classify_cells(cds_subset, TMI_classifier,
                        db = org.Ce.eg.db::org.Ce.eg.db,
                        cluster_extend = TRUE,
                        cds_gene_id_type = "ENSEMBL")


assigned_type_marker_test_res <- top_markers(cds_subset,
    reference_cells = 1000,
    cores = 32
)

garnett_markers <- assigned_type_marker_test_res %>%
    filter(marker_test_q_value < 0.01 & specificity >= 0.27) %>%
    group_by(cell_group) %>%
    top_n(20, marker_score)
# Exclude genes that are good markers for more than one cell type:

garnett_markers <- garnett_markers %>%
    group_by(gene_short_name) %>%
    filter(n() == 1)
