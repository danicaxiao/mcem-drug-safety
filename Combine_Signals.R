rm(list = ls())
library("pROC")
library("RankAggreg")
rank.stat = 'Q_0.05(lambda)'

get_score = function(path, baseline_file){
  #setwd(path)
  path2 <- paste(path,"/",baseline_file,sep="")
  load(path2)
  gpsScore = mgps$ALLSIGNALS
  colnames(gpsScore)[colnames(gpsScore)=="Q_0.05(lambda)"]=rank.stat
  
  path3 <- paste(path, "/","accum.RData",sep="")
  load(path3)
  #next.mgps <- solo.mgps
  
  #files = list.files()
  #print(files)
  #for(i in 1:length(files)){
  #  file = files[i]
   # if(length(grep('accum.RData', file))!=0){
    #  file2 = paste(path,"/",file,sep="")
    #  cmd = paste('load("',file2,'")',sep="")
    #  print(cmd)
    #  eval(parse(text=cmd))
    #}
  
  gpsScore1 = next.mgps$ALLSIGNALS
  
  gpsScore$id = paste(gpsScore$drug,gpsScore$event,sep="_")
  o = order(gpsScore$id)
  gpsScore = gpsScore[o,]
  
  gpsScore1$id = paste(gpsScore1$drug,gpsScore1$event,sep="_")
  o = order(gpsScore1$id)
  gpsScore1 = gpsScore1[o,]
  
  
  ##gold standard
  ref = read.table("/home/cxiao/workspace/omop_adr_drug_concept_cut.tsv", sep="\t",comment="",header=TRUE)

  colnames(ref) = c("drug_id","drugvalue","adr","adr_specific","adr_id","ref")
  ref$adr_id = as.character(ref$adr_id)
  
  ref$adr_id[ref$adr_id=="ACUTE LIVER INJURY"] = "99999999"
  ref$adr_id[ref$adr_id=="GI BLEED"] = "99999998"
  ref$adr_id[ref$adr_id=="ACUTE MYOCARDIAL INFARCTION"] = "99999997"
  ref$adr_id[ref$adr_id=="ACUTE KIDNEY INJURY"] = "99999996"
  ref$id = paste(ref$drug_id,ref$adr_id,sep="_")
  
  subGPSScore = gpsScore[gpsScore$id%in%ref$id,]
  subGPSScore1 = gpsScore1[gpsScore1$id%in%ref$id,]
  
  refGPS = merge(ref,subGPSScore)
  refGPS = unique(refGPS[,c("id","adr_id","ref",rank.stat,"Se")])
  
  refGPS1 = merge(ref,subGPSScore1)
  refGPS1 = unique(refGPS1[,c("id","adr_id","ref",rank.stat,"Se")])
  
  
  return_list <- list("refGPS1"=refGPS1, "refGPS"=refGPS)
  return(return_list)
}

get_signal_info = function(path){
  baseline_file <- "faers.RData"
  return_list <- get_score(path, baseline_file)
  refGPS <- return_list$refGPS1
  refGPS_base <- return_list$refGPS
  score_base <- refGPS_base$`Q_0.05(lambda)`
  var_base <- refGPS_base$Se
  score <- refGPS$`Q_0.05(lambda)`
  var <- refGPS$Se
  id_base <- refGPS_base$id
  id <- refGPS$id
  return_list <- list("ref_base"=refGPS_base$ref, "score_base"=score_base, "var_base"=var_base, "id_base"=id_base, "ref"=refGPS$ref, "score"=score, "var"=var, "id"=id)
  return(return_list)
}

get_auc_individual = function(path_list){
  path_base <- "/home/cxiao/workspace/cloud"
  for(i in 1:length(path_list)){
    print(i)
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    ref_base <- return_list$ref_base
    score_base <- return_list$score_base
    ref <- return_list$ref
    score <- return_list$score
    print(roc(ref_base, score_base, percent=TRUE, plot=TRUE))
    print(roc(ref, score, percent=TRUE, plot=TRUE))
  }
}


combine_signals_kdd13 = function(signals, vars, selection){
  num_signals <- length(signals)
  denominator <- 0
  nominator <- 0
  if(selection == 0){
    for(i in 1:num_signals){
      denominator <- denominator + 1/as.numeric(unlist(vars[i]))
      nominator <- nominator + as.numeric(unlist(signals[i])) / as.numeric(unlist(vars[i]))
    }
  }else{
    for(i in 1:num_signals){
      nominator <- nominator + as.numeric(unlist(signals[i]))
    }
    denominator <- num_signals
  }
  yj_bar <- nominator / denominator
  sj <- var(yj_bar)
  J <- length(yj_bar)
  
  theta <- sum(yj_bar) / J
  diff <- (yj_bar - theta)^2 - sj
  tau <- sum(diff) / J
  
  Bj <- tau / (tau + sj)
  
  mu <- Bj * yj_bar + (1-Bj) * theta
  
  return(mu)
}

combine_signals_avg = function(signals, vars, selection){
  num_signals <- length(signals)
  denominator <- 0
  nominator <- 0
  if(selection == 0){
    for(i in 1:num_signals){
      denominator <- denominator + 1/as.numeric(unlist(vars[i]))
      nominator <- nominator + as.numeric(unlist(signals[i])) / as.numeric(unlist(vars[i]))
    }
  }else{
    for(i in 1:num_signals){
      nominator <- nominator + as.numeric(unlist(signals[i]))
    }
    denominator <- num_signals
  }
  mu <- nominator / denominator
  
  return(mu)
}

combine_signals_ours = function(signals1, vars1, signals2, vars2){
  num_signals1 <- length(signals1)
  num_signals2 <- length(signals2)
  denominator1 <- 0
  denominator2 <- 0
  nominator1 <- 0
  nominator2 <- 0
  median1 <- numeric(num_signals1)
  median2 <- numeric(num_signals2)
  for(i in 1:num_signals1){
    nominator1 <- nominator1 + as.numeric(unlist(signals1[i])) / as.numeric(unlist(vars1[i]))
    denominator1 <- denominator1 + 1/as.numeric(unlist(vars1[i]))
    median1[i] <- mean(as.numeric(unlist(vars1[i])))#median(as.numeric(unlist(vars1[i]))) 
  }  
  for(i in 1:num_signals2){
    nominator2 <- nominator2 + as.numeric(unlist(signals2[i])) / as.numeric(unlist(vars2[i]))
    denominator2 <- denominator2 + 1/as.numeric(unlist(vars2[i]))
    median2[i] <- mean(as.numeric(unlist(vars2[i])))#median(as.numeric(unlist(vars2[i])))
  }
  
  yj_bar1 <- nominator1 / denominator1
  yj_bar2 <- nominator2 / denominator2
  yj_bar <- yj_bar1 / mean(median1) + yj_bar2 / mean(median2)
  
  
  sj <- var(yj_bar)
  J <- length(yj_bar)
  
  theta <- sum(yj_bar) / J
  diff <- (yj_bar - theta)^2 - sj
  tau <- sum(diff) / J
  
  Bj <- tau / (tau + sj)
  
  mu <- Bj * yj_bar + (1-Bj) * theta
  
}

get_auc_combined = function(path_list){
  path_base <- "/home/cxiao/workspace/cloud"
  combine_list_base <- list()
  combine_list <- list()
  for (i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    combine_list_base[[i]] = return_list$id_base
    combine_list[[i]] = return_list$id
  }
  
  signal_intersect_base <- Reduce(intersect,combine_list_base)
  signal_intersect <- Reduce(intersect,combine_list)
  signals_base_list <- list()
  vars_base_list <- list()
  signals_list <- list()
  vars_list <- list()
  
  for(i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    ref_base <- return_list$ref_base
    score_base <- return_list$score_base
    var_base <- return_list$var_base
    id_base <- return_list$id_base
    ref <- return_list$ref
    score <- return_list$score
    var <- return_list$var
    id <- return_list$id
    score_base <- score_base[id_base %in% signal_intersect_base]
    var_base <- var_base[id_base %in% signal_intersect_base]
    score <- score[id %in% signal_intersect]
    var <- var[id %in% signal_intersect]
    signals_base_list[[i]] <- score_base
    signals_list[[i]] <- score
    vars_base_list[[i]] <- var_base
    vars_list[[i]] <- var
  }
  ref_base <- ref_base[id_base %in% signal_intersect_base]
  ref <- ref[id %in% signal_intersect]
  combined_score_base_kdd13_avg <- combine_signals_kdd13(signals_base_list, vars_base_list,1)
  combined_score_base_kdd13_wavg <- combine_signals_kdd13(signals_base_list, vars_base_list,0)
  combined_score_base_avg <- combine_signals_avg(signals_base_list, vars_base_list,1)
  combined_score_base_wavg <- combine_signals_avg(signals_base_list, vars_base_list,0)
  
  combined_score_kdd13_avg <- combine_signals_kdd13(signals_list, vars_list,1)
  combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list, vars_list,0)
  combined_score_avg <- combine_signals_avg(signals_list, vars_list,1)
  combined_score_wavg <- combine_signals_avg(signals_list, vars_list,0)
  
  print("kdd13_avg")
  print(roc(ref_base, combined_score_base_kdd13_avg, percent=TRUE, plot=TRUE))
  print(roc(ref, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  
  print("kdd13_wavg")
  print(roc(ref_base, combined_score_base_kdd13_wavg, percent=TRUE, plot=TRUE))
  print(roc(ref, combined_score_kdd13_wavg, percent=TRUE, plot=TRUE))
  
  print("avg")
  print(roc(ref_base, combined_score_base_avg, percent=TRUE, plot=TRUE))
  print(roc(ref, combined_score_avg, percent=TRUE, plot=TRUE))
  
  print("wavg")
  print(roc(ref_base, combined_score_base_wavg, percent=TRUE, plot=TRUE))
  print(roc(ref, combined_score_wavg, percent=TRUE, plot=TRUE))

}

get_auc_combined_claims = function(path_list, claims_path){
  path_base <- "/home/cxiao/workspace/cloud"
  combine_list <- list()
  for (i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    combine_list[[i]] = return_list$id
  }
  setwd(claims_path)
  claims.data <- "claims_403" 
  data <- read.table(claims.data,header=TRUE, sep="|", comment="")
  data <- as.data.frame(data)
  
  claims.data2 <- "claims_404" 
  data2 <- read.table(claims.data2,header=TRUE, sep="|", comment="")
  data2 <- as.data.frame(data2)
  
  data$CONDITION_CONCEPT_NAME <- as.character(data$CONDITION_CONCEPT_NAME)
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Liver Failure", data$CONDITION_CONCEPT_NAME)] <- "99999999"
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Upper GI Ulcer Hospitalization", data$CONDITION_CONCEPT_NAME)] <- "99999998"
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Acute myocardial Infarction", data$CONDITION_CONCEPT_NAME)] <- "99999997"
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Renal Failure", data$CONDITION_CONCEPT_NAME)] <- "99999996"
  
  data2$CONDITION_CONCEPT_NAME <- as.character(data2$CONDITION_CONCEPT_NAME)
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Liver Failure", data2$CONDITION_CONCEPT_NAME)] <- "99999999"
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Upper GI Ulcer Hospitalization", data2$CONDITION_CONCEPT_NAME)] <- "99999998"
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Acute myocardial Infarction", data2$CONDITION_CONCEPT_NAME)] <- "99999997"
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Renal Failure", data2$CONDITION_CONCEPT_NAME)] <- "99999996"
  
  data$id = paste(data$DRUG_CONCEPT_ID,data$CONDITION_CONCEPT_NAME,sep="_")
  o = order(data$id)
  data = data[o,]
  var_claims <- data$SELOGRR
  score_claims_LB95RR <- data$LB95RR
  ref_claims <- data$GROUND_TRUTH
  
  data2$id = paste(data2$DRUG_CONCEPT_ID,data2$CONDITION_CONCEPT_NAME,sep="_")
  o2 = order(data2$id)
  data2 = data2[o2,]
  var_claims2 <- data2$SELOGRR
  score_claims_LB95RR2 <- data2$LB95RR
  ref_claims2 <- data2$GROUND_TRUTH
  
  combine_list1 <- combine_list
  combine_list2 <- combine_list
  combine_list3 <- combine_list
  
  combine_list1[[length(path_list)+1]] <- data$id
  combine_list2[[length(path_list)+1]] <- data2$id
  combine_list3[[length(path_list)+1]] <- data$id
  combine_list3[[length(path_list)+2]] <- data2$id
  
  signal_intersect1 <- Reduce(intersect,combine_list1)
  signal_intersect2 <- Reduce(intersect,combine_list2)
  signal_intersect3 <- Reduce(intersect,combine_list3)
  
  signals_list1 <- list()
  vars_list1 <- list()
  signals_list2 <- list()
  vars_list2 <- list()
  signals_list3 <- list()
  vars_list3 <- list()
  
  for(i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    ref <- return_list$ref
    score <- return_list$score
    var <- return_list$var
    id <- return_list$id
    
    score1 <- score[id %in% signal_intersect1]
    var1 <- var[id %in% signal_intersect1]
    signals_list1[[i]] <- score1
    vars_list1[[i]] <- var1
    
    score2 <- score[id %in% signal_intersect2]
    var2 <- var[id %in% signal_intersect2]
    signals_list2[[i]] <- score2
    vars_list2[[i]] <- var2
    
    score3 <- score[id %in% signal_intersect3]
    var3 <- var[id %in% signal_intersect3]
    signals_list3[[i]] <- score3
    vars_list3[[i]] <- var3
  }
  ref1 <- ref[id %in% signal_intersect1]
  ref2 <- ref[id %in% signal_intersect2]
  ref3 <- ref[id %in% signal_intersect3]
  
  score_claims_LB95RR_1 <- score_claims_LB95RR[data$id %in% signal_intersect1]
  var_claims_1 <- var_claims[data$id %in% signal_intersect1]
  
  score_claims_LB95RR2_1 <- score_claims_LB95RR2[data2$id %in% signal_intersect2]
  var_claims2_1 <- var_claims2[data2$id %in% signal_intersect2]
  
  score_claims_LB95RR_2 <- score_claims_LB95RR[data$id %in% signal_intersect3]
  var_claims_2 <- var_claims[data$id %in% signal_intersect3]
  
  score_claims_LB95RR2_2 <- score_claims_LB95RR2[data2$id %in% signal_intersect3]
  var_claims2_2 <- var_claims2[data2$id %in% signal_intersect3]
  
  #Signal combination
  claims_signal_list1 <- list(score_claims_LB95RR_1)
  claims_var_list1 <- list(var_claims_1)
  claims_signal_list2 <- list(score_claims_LB95RR2_1)
  claims_var_list2 <- list(var_claims2_1)
  claims_signal_list3 <- list(score_claims_LB95RR_2, score_claims_LB95RR2_2)
  claims_var_list3 <- list(var_claims_2, var_claims2_2)
  
  combined_score_ours1 <- combine_signals_ours(signals_list1, vars_list1, claims_signal_list1, claims_var_list1)
  print("ours")
  print("Faers + 403")
  print(roc(ref1, combined_score_ours1, percent=TRUE, plot=TRUE))
  
  combined_score_ours2 <- combine_signals_ours(signals_list2, vars_list2, claims_signal_list2, claims_var_list2)
  print("Faers + 404")
  print(roc(ref2, combined_score_ours2, percent=TRUE, plot=TRUE))
  
  combined_score_ours3 <- combine_signals_ours(signals_list3, vars_list3, claims_signal_list3, claims_var_list3)
  print("Faers + 404 + 403")
  print(roc(ref3, combined_score_ours3, percent=TRUE, plot=TRUE))
  

  signals_list1[[length(path_list)+1]] <- score_claims_LB95RR_1
  vars_list1[[length(path_list)+1]] <- var_claims_1
  
  signals_list2[[length(path_list)+1]] <- score_claims_LB95RR2_1
  vars_list2[[length(path_list)+1]] <- var_claims2_1
  
  signals_list3[[length(path_list)+1]] <- score_claims_LB95RR_2
  signals_list3[[length(path_list)+2]] <- score_claims_LB95RR2_2
  vars_list3[[length(path_list)+1]] <- var_claims_2
  vars_list3[[length(path_list)+2]] <- var_claims2_2
  
  
  combined_score_kdd13_avg <- combine_signals_kdd13(signals_list1, vars_list1,1)
  combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list1, vars_list1,0)
  combined_score_avg <- combine_signals_avg(signals_list1, vars_list1,1)
  combined_score_wavg <- combine_signals_avg(signals_list1, vars_list1,0)
  print("Faers + 403")
  print("kdd13_avg")
  print(roc(ref1, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  print("kdd13_wavg")
  print(roc(ref1, combined_score_kdd13_wavg, percent=TRUE, plot=TRUE))
  print("avg")
  print(roc(ref1, combined_score_avg, percent=TRUE, plot=TRUE))
  print("wavg")
  print(roc(ref1, combined_score_wavg, percent=TRUE, plot=TRUE))
  
  combined_score_kdd13_avg <- combine_signals_kdd13(signals_list2, vars_list2,1)
  combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list2, vars_list2,0)
  combined_score_avg <- combine_signals_avg(signals_list2, vars_list2,1)
  combined_score_wavg <- combine_signals_avg(signals_list2, vars_list2,0)
  print("Faers + 404")
  print("kdd13_avg")
  print(roc(ref2, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  print("kdd13_wavg")
  print(roc(ref2, combined_score_kdd13_wavg, percent=TRUE, plot=TRUE))
  print("avg")
  print(roc(ref2, combined_score_avg, percent=TRUE, plot=TRUE))
  print("wavg")
  print(roc(ref2, combined_score_wavg, percent=TRUE, plot=TRUE))
  
  combined_score_kdd13_avg <- combine_signals_kdd13(signals_list3, vars_list3,1)
  combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list3, vars_list3,0)
  combined_score_avg <- combine_signals_avg(signals_list3, vars_list3,1)
  combined_score_wavg <- combine_signals_avg(signals_list3, vars_list3,0)
  print("Faers + 403 + 404")
  print("kdd13_avg")
  print(roc(ref3, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  print("kdd13_wavg")
  print(roc(ref3, combined_score_kdd13_wavg, percent=TRUE, plot=TRUE))
  print("avg")
  print(roc(ref3, combined_score_avg, percent=TRUE, plot=TRUE))
  print("wavg")
  print(roc(ref3, combined_score_wavg, percent=TRUE, plot=TRUE))
  
}

################################ FAERS #######################################

path_list <- c(2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014)
claims_path <- "/home/cxiao/workspace/"

get_auc_individual(path_list)

get_auc_combined(path_list)

get_auc_combined_claims(path_list, claims_path)



rank_aggre_function = function(path_list){
  path_base <- "/home/cxiao/workspace/cloud"
  combine_list <- list()
  for (i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    combine_list[[i]] = return_list$id
  }
  setwd(claims_path)
  claims.data <- "claims_403" 
  data <- read.table(claims.data,header=TRUE, sep="|", comment="")
  data <- as.data.frame(data)
  
  claims.data2 <- "claims_404" 
  data2 <- read.table(claims.data2,header=TRUE, sep="|", comment="")
  data2 <- as.data.frame(data2)
  
  data$CONDITION_CONCEPT_NAME <- as.character(data$CONDITION_CONCEPT_NAME)
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Liver Failure", data$CONDITION_CONCEPT_NAME)] <- "99999999"
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Upper GI Ulcer Hospitalization", data$CONDITION_CONCEPT_NAME)] <- "99999998"
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Acute myocardial Infarction", data$CONDITION_CONCEPT_NAME)] <- "99999997"
  data$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Renal Failure", data$CONDITION_CONCEPT_NAME)] <- "99999996"
  
  data2$CONDITION_CONCEPT_NAME <- as.character(data2$CONDITION_CONCEPT_NAME)
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Liver Failure", data2$CONDITION_CONCEPT_NAME)] <- "99999999"
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Upper GI Ulcer Hospitalization", data2$CONDITION_CONCEPT_NAME)] <- "99999998"
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Acute myocardial Infarction", data2$CONDITION_CONCEPT_NAME)] <- "99999997"
  data2$CONDITION_CONCEPT_NAME[grepl("OMOP Acute Renal Failure", data2$CONDITION_CONCEPT_NAME)] <- "99999996"
  
  data$id = paste(data$DRUG_CONCEPT_ID,data$CONDITION_CONCEPT_NAME,sep="_")
  o = order(data$id)
  data = data[o,]
  var_claims <- data$SELOGRR
  score_claims_LB95RR <- data$LB95RR
  ref_claims <- data$GROUND_TRUTH
  
  data2$id = paste(data2$DRUG_CONCEPT_ID,data2$CONDITION_CONCEPT_NAME,sep="_")
  o2 = order(data2$id)
  data2 = data2[o2,]
  var_claims2 <- data2$SELOGRR
  score_claims_LB95RR2 <- data2$LB95RR
  ref_claims2 <- data2$GROUND_TRUTH
  
  combine_list3 <- combine_list
  
  combine_list3[[length(path_list)+1]] <- data$id
  combine_list3[[length(path_list)+2]] <- data2$id
  
  signal_intersect3 <- Reduce(intersect,combine_list3)
  
  scores_list3 <- list()
  signals_list3 <- list()
  ref_list3 <- list()
  
  signal_matrix <- matrix(, nrow=13, ncol=194)
  
  for(i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    ref <- return_list$ref
    score <- return_list$score
    var <- return_list$var
    id <- return_list$id
    
    score3 <- score[id %in% signal_intersect3]
    ref3 <- ref[id %in% signal_intersect3]
    o3 <- order(score3,decreasing = TRUE)
    score3 <- score3[o3]
    signal3 <- signal_intersect3[o3]
    ref3 <- ref3[o3]
    
    scores_list3[[i]] <- score3
    ref_list3[[i]] <- ref3
    signals_list3[[i]] <- signal3
    
    signal_matrix[i,1:194] <- signal3
    
  }
  
  score_claims_LB95RR <- score_claims_LB95RR[data$id %in% signal_intersect3]
  o_LB95RR <- order(score_claims_LB95RR,decreasing = TRUE)
  score_claims_LB95RR <- score_claims_LB95RR[o_LB95RR]
  signal_LB95RR <- signal_intersect3[o_LB95RR]
  
  signal_matrix[(length(path_list)+1),1:194] <- signal_LB95RR
  
  score_claims_LB95RR2 <- score_claims_LB95RR2[data2$id %in% signal_intersect3]
  o_LB95RR2 <- order(score_claims_LB95RR2,decreasing = TRUE)
  score_claims_LB95RR2 <- score_claims_LB95RR2[o_LB95RR2]
  signal_LB95RR2 <- signal_intersect3[o_LB95RR2]
  
  signal_matrix[(length(path_list)+2),1:194] <- signal_LB95RR2
  return(signal_matrix)
}	


### Rank aggregation part can be ignored
signal_matrix <- rank_aggre_function(path_list)
save(signal_matrix, file = "signal_matrix.RData")

aggre <- RankAggreg(score_matrix, 246)

setwd("/home/cxiao/workspace/Signal_Combination")
load("aggre_result.RData")
rank_list <- aggre$top.list

match_index <- match(signal1,rank_list)
ranked_score <- score1_prev[match_index]
ranked_ref <- ref1_prev[match_index]

roc(ranked_ref, ranked_score, percent=TRUE, plot=TRUE)
