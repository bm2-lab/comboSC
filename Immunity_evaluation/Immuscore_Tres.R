
countToTpm <- function(count,effLen){
    rate <- log(count) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}


classified <- function(eva,Sample_id = "su001") {
    sl <- eva[,2]
    score <- eva[eva[,1] == Sample_id,2]
    line <- sl[round(length(sl)*0.9)]
    if(score > line){
        Patient_class <- 2
    }else if (score > 0) {
        Patient_class <- 1
    }else if (score < 0) {
         Patient_class <- 0
    }
    return(Patient_class)
}


Tre_pre <- function(pbmc,cancertype = "BCC",Sample_id = "008") {
    tumor_col <- c("CD8_mem_T_cells","CD8_act_T_cells","CD8_ex_T_cells","CD4_T_cells","Tcell_prolif","Tregs")
    cancer_type <- c("BCC" ,"BRCA" ,"CRC" ,"HNSC" ,"NSCLC" ,"PAAD" ,"SKCM" ,"LIHC" ,"AML" ,"UCEC")
    pbmc <- subset(pbmc, subset = cluster %in% tumor_col)
    gen <- read.table("./Resources/All_hg19gene_len.txt",header = T,sep="\t")
    rownames(gen) <- gen$Gene
    if(length(unique(pbmc@meta.data$patient))== 1){
        if(pbmc@meta.data$cancertype[2] %in% cancer_type){
            Texp_path <- paste0("./Immunity_evaluation/Reference/",pbmc@meta.data$cancertype[2],"_texp.txt")
            Texp = fread(Texp_path, data.table = FALSE)
        }else{
            Texp = fread("./Immunity_evaluation/Reference/ref_gen.txt", data.table = FALSE)
        }
        gene<-intersect(rownames(pbmc@assays$RNA@counts),Texp[,1])
        expr <- pbmc@assays$RNA@counts[gene,]
        len<-gen[gene,]
        tpms <- as.data.frame(expr)
        for(i in c(1:dim(expr)[2])){
            tpms[,i] <- countToTpm(expr[,i],effLen = len$Length)
        }
        nexp <- log2(tpms/10+1)
        nex <- apply(nexp, 1, mean)
        Texp[Sample_id] <- nex
     } else {
    sn = table(pbmc@meta.data$patient) > 20
    print(sn)
    plist <- unique(pbmc@meta.data$patient)[sn]
    all_exp <- pbmc@assays$RNA@counts
    gene<-intersect(rownames(all_exp),rownames(gen))
    all_exp<-all_exp[gene,]
    len<-gen[gene,]
    Texp <- matrix(data = 0, nrow = dim(all_exp)[1], ncol = plist %>% length() + 1, byrow = FALSE,
                dimnames = list(rownames(all_exp),c("Gene",plist))) %>% as.data.frame()

    Texp["Gene"] <- rownames(all_exp)
    for(p in plist){
        print(p)
        su <- subset(pbmc, subset = patient == p)
        expr <- su@assays$RNA@counts[gene,]
        tpms <- as.data.frame(apply(expr, 2, countToTpm, effLen = len$Length))
        nexp <- log2(tpms/10+1)
        nex <- apply(nexp, 1, mean)
        Texp[p] <- nex
    }}
    output <- "./Immunity_evaluation/Input/ls.txt"
    write.table(Texp,output,sep = "\t",row.names = FALSE,quote = FALSE)
    output2 <- "./Immunity_evaluation/Output/ls.txt"
    command = paste("Tres.py -m predict -i ",output," -n 1 -o ",output2,sep = "")
    system(command)
    eva <-  fread(output2, data.table = FALSE)
    return(eva)
}

Immune_score <- function(pbmc = pbmc, Sample_id = "su008"){
    score_df <- Tre_pre(pbmc,Sample_id = Sample_id)
    Patient_class <- classified(score_df, Sample_id)
}