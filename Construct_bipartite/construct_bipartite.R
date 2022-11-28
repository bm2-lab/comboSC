# Author:Tang Chen
# Email: chentang@tongji.edu.cn
# date: 2020-02-24 10:20:47
# --------------
# About project:

# save.image('./resources/rxhaut_vaccine.Rdata')
# load('./resources/rxhaut_vaccine.Rdata')
construct_bipartite <- function(patient_class,
                                pbmc,
                                res_rank = seq(0.4, 3, 0.2),
                                there = 0.5,
                                rfgene = rfgene_rpkm,
                                essential_genes = essential_genes,
                                group2 = "CD8_act_T_cells",
                                group1 = "CD8_ex_T_cells") {
    if(patient_class == 2){
        drug_df <- c()
    } 
    if (patient_class == 1) {
        cat("Drug combination paired with ICB therapy is selected.\n")
        drug_df <- exhaust_caccine(pbmc = pbmc,
                            group2 = "CD8_act_T_cells",
                            group1 = "CD8_ex_T_cells")
    }
    if(patient_class == 0) {
        cat("Drug combination for targeting therapy is selected.\n")
        drug_df <- decluster(pbmc = pbmc,
                    res_rank = res_rank,
                    there = there,
                    rfgene = rfgene_rpkm,
                    essential_genes = essential_genes)
    }

    cat('###################################################\nStep3 is OK!\n#################################################\n')
    return(drug_df)
}


###################################################################
######################--------example---------######################
###################################################################

# drug_df <- construct_bipartite(patient_class = 1,
#                                 exp = pexp,
#                                 metadata = pt1_metadata,
#                                 res_rank = seq(0.4, 5, 0.2),
#                                 there = 0.5,
#                                 rfgene = rfgene_rpkm,
#                                 essential_genes = essential_genes,
#                                 group2 = "CD8_act_T_cells",
#                                 group1 = "CD8_ex_T_cells")
