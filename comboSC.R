#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

ex <- args[1]
meta <- args[2]
resolution_rank <- args[3]
threshold <- args[4]
rfgene <- args[5]
Sample_id <- args[7]

# library(getopt)
# 首先将第一个参数的具体信息进行描述
# 每行五个，第五个可选，也就是说第五列可以不写
# byrow 按行填充矩阵的元素
# ncol  每行填充五个元素
# spec <- matrix(
#   c("first",  "f", 2, "integer", "This is first!",
#     "second", "s", 1, "character",  "This is second!",
#     "third",  "t", 2, "double",  "This is third!",
#     "help",   "h", 0, "logical",  "This is Help!"),
#   byrow=TRUE, ncol=5)
# # 使用getopt方法
# opt <- getopt(spec=spec)
# # opt实际上就是一个列表，直接使用$来索引到对应的参数的值
# print(opt$first)
# print(opt$second)
# print(opt$third)

setwd("/home/tangchen/code/comboSC/")
source("./resources/Allpackage.R")
source("./resources/dep_data.R")
source("./Data_preprocessing/scprocess.R")
source("./immunity_evaluation/immuscore.R")
source("./intelligent_optimization/graph_optimization.R")
source("./Construct_bipartite/intratumor/decluster.R")
source("./Construct_bipartite/TMI/exhaut_vaccine.R")
source("./Construct_bipartite/construct_bipartite.R")

ex <- "/home/tangchen/code/precision_medicine/Sctc/resources/bcc_exp_006.csv.gz"
meta <- "/home/tangchen/code/precision_medicine/Sctc/resources/bcc_meta_006.csv.gz"



exp <-  fread(ex,data.table = FALSE)
metadata <- fread(meta,data.table = FALSE)
rownames(exp) =unlist(exp[, 1])
exp = exp[, -1]
metadata <- metadata[, -1]
rownames(metadata) =unlist(metadata[, 1])

pbmc <- readRDS("/home/tangchen/scRNA/RDS/bcc.rds")
# pbmc <- subset(pbmc, subset = patient == "su001")

COMBOSC <- function(exp, metadata,
                res_rank = seq(0.4, 3, 0.2),
                there = 0.5,
                rfgene = rfgene_rpkm,
                Output_line = 50,
                Sample_id = "006") {
                # essential_genes = essential_genes) 
    pbmc <- scrna_preprocess(exp, metadata)
    Patient_class <- Tre_pre(pbmc,Sample_id)           
    drug_df <- construct_bipartite(patient_class = Patient_class, pbmc, res_rank , there, rfgene , essential_genes)
    solution_recommended <- graph_optimization(mel_top, sin_res, patient_class = Patient_class, output_line = Output_line, sample_id = Sample_id)
}  

solution_recommended <- COMBOSC(exp, metadata,
                        res_rank = seq(0.4, 3, 0.2),
                        there = threshold,
                        rfgene = rfgene_rpkm,
                        Output_line = 100,
                        Sample_id = Sample_id)

