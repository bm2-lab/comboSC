
countToTpm <- function(count,effLen){
    rate <- log(count) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

Tre_pre <- function(pbmc,patient_id = "su001") {
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
        su <- subset(pbmc, subset = patient == p)
        expr <- su@assays$RNA@counts[gene,]
        tpms <- as.data.frame(apply(expr, 2, countToTpm, effLen = len$Length))
        nexp <- log2(tpms/10+1)
        nex <- apply(nexp, 1, mean)
        Texp[p] <- nex
    }
    # output <- paste("/home/tangchen/code/comboSC/immunity_evaluation/Input/",strsplit(ex,"_")[[1]][2],"_texp.txt",sep = "")
    output <- "./immunity_evaluation/Input/ls.txt"
    write.table(Texp,output,sep = "\t",row.names = FALSE,quote = FALSE)
    # output2 = paste("/home/tangchen/code/comboSC/immunity_evaluation/Output/tre_",basename(output),sep = "")
    output2 <- "./immunity_evaluation/Output/ls.txt"
    command = paste("./immunity_evaluation/Tres.py -m predict -i ",output," -n 1 -o ",output2,sep = "")
    system(command)
    eva <-  fread(output2, data.table = FALSE)
    return(eva)
}

classified <- function(eva,patient_id = "su001") {
    sl <- eva[,2]
    score <- eva[eva[,1] == patient_id,2]
    line <- sl[round(length(sl)*0.9)]
    if(score > line){
        Patient_class <- 2
    }else if (score > 0) {
        Patient_class <- 1
    }else{
        Patient_class <- 1
    }
    return(Patient_class)
}

Immune_score <- function(pbmc, patient_id = "su001"){
    gen <- read.table("./resources/All_hg19gene_len.txt",header = T,sep="\t")
    rownames(gen)<-gen$Gene
    pbmc <- Tre_pre(pbmc)
    Patient_class <- classified(pbmc, patient_id)
}
# ###single cell
# pbmc.list <- c(pbmc, pbmc2)
# merged.anchors <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:30)

# pbmc <- IntegrateData(anchorset = merged.anchors, dims = 1:30)     #进行数据集整合
# Tre_pre(pbmc)

# pbmc <- readRDS("/home/tangchen/scRNA/RDS/bcc.rds")
# pbmc <- subset(pbmc, subset = patient == "su008")
# pbmc@meta.data$Patient <- "008"
# # pbmc <- subset(pbmc, subset = treatment == "pre")
# pbmc@meta.data$Patient <- pbmc@meta.data$patient

# Tre_pre(pbmc,Sample_id)
