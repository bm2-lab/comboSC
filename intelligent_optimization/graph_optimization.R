get_cluster <- function(identiy, mel_top, cell_opted,velocity = FALSE) {
    #' get the weight and the number of the every clusters
    cluster_number <- data.frame(table(identiy[, 1]))
    cluster_number$freq <- cluster_number$Freq / sum(cluster_number$Freq)
    colnames(cluster_number) <- c("cluster", "count", "freq")
    cluster_number_freq <- cluster_number$freq
    names(cluster_number_freq) <- cluster_number[, 1]
    cell_opted_weight <- unlist(lapply(colnames(mel_top), function(x) {
        cluster_number_freq[x]
    }))
    if(velocity){
        cell_opted_weight <- velocity * cell_opted_weight
    }
    names(cell_opted_weight) <- cell_opted
    return(cell_opted_weight)
}

drug_filter <- function(cell_opted, sen_there, sin_res = sin_res) {
    # select cell lines, remove low activity drugs
    # Screening for low reactivity drugs
    m <- sin_res[abs(drugsing_res$S) < sen_there, ]
    n <- data.frame(m[m$cell_line_name %in% cell_opted, ])
    drug_used <- data.frame(
        "drug" = n$drug_comb,
        "cell_line" = n$cell_line_name, "css" = n$css,
        stringsAsFactors=FALSE
    )
    cluster <- unique(drug_used$cell_line)
    drug <- unique(drug_used$drug)
    bipart <- matrix(0, length(drug), length(cluster), dimnames = list(drug,cluster))
    for(i in c(1:length(drug))){
        bipart[as.character(drug_used[i,1]),as.character(drug_used[i,2])] <- drug_used[i,3]
    }
    bipart01 <- bipart
    bipart01[bipart01 > 0] <- 1
    bipart01[bipart01 <= 0] <- 0
    bipart01 <- bipart01[rowSums(bipart01[,]) > round(length(drug)*0.3),]
    drug_used <- drug_used[drug_used$drug %in% rownames(bipart01), ]
    tox <- fread("./resources/drugtox_names.csv",data.table = FALSE)
    drug_used[drug_used[,1] %in% tox[,2],"css"] <- drug_used[drug_used[,1] %in% tox[,2],"css"] * 0.9
    return(drug_used)
}

sum_max <- function(drug01, drug02, cell_opted, cell.opted.weight) {
    drug01_line <- drug_used[drug_used$drug == drug01, ]
    drug02_line <- drug_used[drug_used$drug == drug02, ]
    value <- 0
    for (CellLine in cell_opted) {
        v1 <- if (CellLine %in% drug01_line$cell_line) drug01_line[drug01_line$cell_line == CellLine, 3] else 0
        v2 <- if (CellLine %in% drug02_line$cell_line) drug02_line[drug02_line$cell_line == CellLine, 3] else 0
        value <- value + cell.opted.weight[CellLine]*(v1 + v2)
    }

    return(value)
}


cluster_weight <- function(df, cell_opted_weight) {
    #' put weight on every clusters
    cell_opted_weight <- na.omit(cell_opted_weight)
    for (i in c(1:length(df))) {
        df[i, 3] <- df[i, 3] * cell_opted_weight[as.character(df[i, 2])]
    }
    df <- na.omit(df)
    df <- df[!duplicated(df[, c(1, 2)]), ]
    return(df)
}




drug_scorepre <- function(drug_used) {
    all_drug <- as.vector(unique(drug_used$drug))
    drug_score <- rep(0, length(all_drug))
    names(drug_score) <- all_drug
    for (dg in all_drug) {
        cluster_unweighted <- as.vector(drug_used[drug_used$drug == dg, ])
        cluster_weighted <- cluster_weight(cluster_unweighted,
                                            cell_opted_weight)
        drug_score[dg] <- sum(cluster_weighted[, 3])
    }
    return(drug_score)
}

cut_step <- function(drug_score, cut_point, cut_point_start = 1) {
#' sort the data and get the corresponding interval
#    all_drug <- as.vector(unique(drug_used$drug))
    all_drug <- names(drug_score)
    drug.order <- order(drug_score, decreasing = TRUE)
    lefted_num <- drug.order[cut_point_start :
                        (cut_point_start + cut_point * length(drug.order))]
    drug_lefted <- all_drug[lefted_num]
    return(drug_lefted)
}

cell3to1 <- function(mel_top,norm_drug) {
#' Select the most suitable celllines from top3
    group <- c()
    Col <- 1
    while (Col <= ncol(mel_top)) {
        Row <- 1
        while (Row <= nrow(mel_top)) {
            #print(c(Col,Row))
                a <- mel_top[Row, Col]
                rd <- eval(parse(text = paste("norm_drug$ ", a, sep = "")))
                #print(paste("norm_drug$", a, sep = ""))
                if (length(rd) > 0){
                    group <- c(group, rd)
                    break
                }
                Row <- Row + 1 
            }
        
        Col <- Col + 1
    }
    return(group)
}


sin_scores <- function(drug_lefted01, drug_used, db_res = db_res ,cell_opted_weight = cell_opted_weight, full_ranking = "/home/tangchen/vde/full_ranking.csv") {
#' Greedy algorithm. 
    # rest_drug <<- as.vector(unique(drug_used$drug))
    mid_score <- data.frame("drugcomb" = "test", "score" = 1,stringsAsFactors =FALSE)
    rest_drug <- drug_lefted01
    for (drug1 in drug_lefted01) {
        rest_drug <- rest_drug[-which(rest_drug == drug1)]
        if (length(rest_drug) == 0) {
                break
            }
        for (drug2 in rest_drug) {
            double_drug <- paste(drug1, "+", drug2) 
            na_re<- paste(drug1, " & ", drug2)
            an_er<- paste(drug2, " & ", drug1)

            if (na_re %in% db_res$drug_comb | an_er %in% db_res$drug_comb) {
                css <- 0
                df <- db_res[db_res$drug_comb == double_drug, c(3, 4, 6)]
                df <- cluster_weight(df, cell_opted_weight)
                if (length(df) < 5) {
                    sin_drug <- setdiff(cell_opted, as.vector(df[, 2]))
                    css <- css + sum_max(drug1,
                                        drug2, sin_drug,cell_opted_weight)
                }
                css <- df[, 3] %>% as.vector %>% as.numeric %>% sum
                cat(c(drug1, drug2),'\n')
            } else {
                css <- sum_max(drug1,
                            drug2, cell_opted = cell_opted, cell.opted.weight = cell_opted_weight)
            }
            mid_score <- rbind(mid_score, c(double_drug, css))       
        }
    }
    sin_score <- mid_score[order(mid_score[, 2] %>% as.numeric), ]
    dp <- sin_score[duplicated(sin_score[,2]),2]
    sin_score <- sin_score[!sin_score[,2] %in% dp, ]
    if(full_ranking!=FALSE){
        write.csv(sin_score, full_ranking)
    
    }
    return(sin_score)
}

TMI_opt <- function(bio_matrix,Output_line = 30,sample_id = sample_id, patient_class = patient_class) {
    bio_matrix <- data.frame(bio_matrix, stringsAsFactors = FALSE)
    bio_matrix$`Score value` <- bio_matrix$Connectivity.x + bio_matrix$Connectivity.y
    bio_matrix <- bio_matrix[order(bio_matrix$Score, decreasing=FALSE),]
    bio_matrix$`Drug combination` <- paste0(bio_matrix$gene," + ICI")
    bio_matrix$`Personalized score level` <- "Middle"
    solution_recommended <- bio_matrix[order(bio_matrix$`Score value`, decreasing = FALSE)[1:Output_line],]
    solution_recommended <- solution_recommended[c("Drug combination","Personalized score level", "Score value")]
    result_path <- paste("./auxiliary/",sample_id, ".", patient_class,".solution_recommended.csv", sep = "")
    write.csv(solution_recommended, result_path, row.names = FALSE)
    return(solution_recommended)
}


######################------main--------######################
graph_optimization  <- function(mel_top = mel_top, sin_res = sin_res, 
                                patient_class = 1, output_line = 2 ,
                                cut_point = 1,sample_id = "tumor") {
    if(patient_class == 2){
        result_path <- paste("./auxiliary/",sample_id, ".", 
                                patient_class,".solution_recommended.csv", sep = "")

        solution_recommended <- data.frame("Drug combination" = c("Na"), "Score value" = c("Na"), "Personalized                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               score level" = c("Na"))
        write.csv(solution_recommended, result_path, row.names = FALSE)
    }                                                                                                        
    if(patient_class == 0) {
        solution_recommended <- TMI_opt(bio_matrix,Output_line = output_line,
                                sample_id = sample_id, patient_class = patient_class) 
    } 
    if(patient_class == 1) {
        cell_opted <<- cell3to1(mel_top = mel_top, norm_drug)
        drug_used <- drug_filter(cell_opted, 1, sin_res)
        drug_used <<-  na.omit(drug_used[!duplicated(drug_used[,c(1,2)]),])
        cell_opted_weight <<- get_cluster(identiy, mel_top, cell_opted)
        drug_score <<- drug_scorepre(drug_used)
        drug_lefted01 <<- cut_step(drug_score, cut_point = 1)  
        drug_lefted01 <- drug_lefted01[which(drug_lefted01 != "NA")]
        cl <- makeCluster(cores, type="FORK")
        cat("\nLast step!\n")
        sin_score <- sin_scores(drug_lefted01, drug_used, db_res,cell_opted_weight,full_ranking = "/home/tangchen/vde/full_ranking.csv") ####问题
        stopCluster(cl)
        sin_score_topn <- solution_ndrugs(sin_score)
        solution_recommended <- sin_score_topn
        Score_list <- c("Middle", "Low")
        solution_recommended[,3] <- Score_list[2]
        colnames(solution_recommended) <- c("DrugCombination","Score value", "Personalized score level")
        solution_recommended <- solution_recommended[c("DrugCombination","Personalized score level", "Score value")]
        solution_recommended[, 3] <- as.character(solution_recommended[, 3]) 
        result_path <- paste("./auxiliary/",sample_id, ".", patient_class,".solution_recommended.csv", sep = "")
        write.csv(solution_recommended, result_path, row.names = FALSE)
        print(command)
        system(command, intern = TRUE)
    }
    cat('###################################################\nLast step is OK!\n#################################################\n')
    return(solution_recommended)
}
