#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

ex <- args[1]
meta <- args[2]
Sample_id <- args[3]
threshold <- args[4]


source("./Resources/Allpackage.R")
source("./Data_preprocessing/scprocess.R")
source("./Immunity_evaluation/Immuscore_Tres.R")
source("./Intelligent_optimization/graph_optimization.R")
source("./Construct_bipartite/intratumor/decluster.R")
source("./Construct_bipartite/TMI/exhaut_vaccine.R")
source("./Construct_bipartite/construct_bipartite.R")
load("./Resources/data_comboSC.rdata")


exp <-  fread(ex,data.table = FALSE)
metadata <- fread(meta,data.table = FALSE)
rownames(exp) =unlist(exp[, 1])
exp = exp[, -1]
metadata <- metadata[, -1]
rownames(metadata) =unlist(metadata[, 1])


COMBOSC <- function(exp, metadata,
                res_rank = seq(0.4, 3, 0.2),
                there = 0.5,
                rfgene = rfgene_rpkm,
                Output_line = 50,
                Sample_id = "006") {
    pbmc <- scrna_preprocess(exp, metadata)
    Patient_class <- Tre_pre(pbmc,Sample_id)           
    drug_df <- construct_bipartite(patient_class = Patient_class, pbmc, res_rank , there, rfgene , essential_genes)
    solution_recommended <- graph_optimization(mel_top, sin_res, patient_class = Patient_class, output_line = Output_line, sample_id = Sample_id)
}  

solution_recommended <- COMBOSC(exp, metadata,
                        res_rank = seq(0.4, 3, 0.2),
                        there = threshold,
                        rfgene = rfgene_rpkm,
                        Output_line = 50,
                        Sample_id = Sample_id)

