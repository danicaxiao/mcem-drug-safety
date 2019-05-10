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
  
  path3 <- paste(path, "/","solo.RData",sep="")
  load(path3)
  next.mgps <- solo.mgps
  
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
  ref = read.table("/Users/inci.baytasibm.com/Box Sync/SummerInternProject2017/data/mcem_mtl/omop_adr_drug_concept_cut.tsv", sep="\t",comment="",header=TRUE)
  #ref = read.table("/home/baytasin/IBM/omop_adr_drug_concept_cut.tsv", sep="\t",comment="",header=TRUE)
  
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
  return_list <- get_score(path,baseline_file)
  refGPS <- return_list$refGPS1
  refGPS_base <- return_list$refGPS
  score <- refGPS$`Q_0.05(lambda)`
  var <- refGPS$Se
  id <- refGPS$id
  ref <- refGPS$ref
  score_base <- refGPS_base$`Q_0.05(lambda)`
  var_base <- refGPS_base$Se
  id_base <- refGPS_base$id
  ref_base <- refGPS_base$ref
  score_base <- score_base[id_base %in% id]
  var_base <- var_base[id_base %in% id]
  ref_base <- ref_base[id_base %in% id]
  
  return_list <- list("ref"=ref, "score"=score, "var"=var, "id"=id, "ref_base"=ref_base, "score_base"=score_base, "var_base"=var_base)
  return(return_list)
}

get_auc_individual = function(path_list){
  path_base <- "/Users/inci.baytasibm.com/Documents/cloud_"
  for(i in 1:length(path_list)){
    print(i)
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    ref <- return_list$ref
    score <- return_list$score
    print(roc(ref, score, percent=TRUE, plot=TRUE))
    ref_base <- return_list$ref_base
    score_base <- return_list$score_base
    print(roc(ref_base, score_base, percent=TRUE, plot=TRUE))
  }
}

get_claims_signal = function(claims_path){
  setwd(claims_path)
  claims_data <- "claims_403"
  claims_path1 <- paste(claims_path,"/",claims_data,sep="")
  data <- read.table(claims_path1,header=TRUE, sep="|", comment="")
  data <- as.data.frame(data)
  
  claims_data2 <- "claims_404" 
  claims_path2 <- paste(claims_path,"/",claims_data2,sep="")
  data2 <- read.table(claims_data2,header=TRUE, sep="|", comment="")
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
  var_403 <- data$SELOGRR
  score_403 <- data$LB95RR
  ref_403 <- data$GROUND_TRUTH
  id_403 <- data$id
  
  data2$id = paste(data2$DRUG_CONCEPT_ID,data2$CONDITION_CONCEPT_NAME,sep="_")
  o2 = order(data2$id)
  data2 = data2[o2,]
  var_404 <- data2$SELOGRR
  score_404 <- data2$LB95RR
  ref_404 <- data2$GROUND_TRUTH
  id_404 <- data2$id
  
  return_list <- list("score_403"=score_403, "var_403"=var_403, "ref_403"=ref_403, "id_403"=id_403, "score_404"=score_404, "var_404"=var_404, "ref_404"=ref_404, "id_404"=id_404)
  return(return_list)
}

EM_Alg = function(theta_ini, tau_ini, sj, yj_bar, J){
  max_iter = 1000000
  theta <- theta_ini
  tau <- tau_ini
  for(i in 1:max_iter){
    theta_prev <- theta
    tau_prev <- tau 
    # E-Step
    T1 <- 0
    T2 <- 0
    Bj <- tau / (tau + sj)
    for(j in 1:J){
      T1 <- T1 + Bj*yj_bar[j] + (1-Bj)*theta
      T2 <- T2 + (Bj * yj_bar[j] + (1-Bj)*theta)^2 + Bj*sj
    }
    # M-Step
    theta <- T1 / J
    tau <- T2 / J - theta^2
    if((abs(theta - theta_prev) + (abs(tau - tau_prev))) < 1e-10){
      print("Algorithm Converged!")
      break
    }
  }
  return(list("theta"=theta, "tau"=tau))
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
  
  # return_list <- EM_Alg(theta_ini, tau_ini, sj, yj_bar, J)
  # theta <- return_list$theta
  # tau <- return_list$tau
  
  Bj <- tau / (tau + sj)
  
  mu <- Bj * yj_bar + (1-Bj) * theta
  
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
  #yj_bar <- yj_bar1 + yj_bar2
  yj_bar <- yj_bar1 / median(median1) + yj_bar2 / median(median2)
  
  
  sj <- var(yj_bar)
  J <- length(yj_bar)
  
  theta <- sum(yj_bar) / J
  diff <- (yj_bar - theta)^2
  tau <- sum(diff) / J
  
  # return_list <- EM_Alg(theta_ini, tau_ini, sj, yj_bar, J)
  # theta <- return_list$theta
  # tau <- return_list$tau
  
  Bj <- tau / (tau + sj)
  
  mu <- Bj * yj_bar + (1-Bj) * theta
  
}

get_auc_combined = function(path_list, claims_path){
  path_base <- "/Users/inci.baytasibm.com/Documents/cloud_"
  combine_list <- list()
  for (i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="")
    return_list <- get_signal_info(path)
    combine_list[[i]] = return_list$id
  }
  
  # Claims + FAERS
  claims_return_list <- get_claims_signal(claims_path)
  
  ref_403 <- claims_return_list$ref_403
  ref_404 <- claims_return_list$ref_404
  var_403 <- claims_return_list$var_403
  var_404 <- claims_return_list$var_404
  score_403 <- claims_return_list$score_403
  score_404 <- claims_return_list$score_404
  id_403 <- claims_return_list$id_403
  id_404 <- claims_return_list$id_404
  
  #combine_list1 <- combine_list
  # combine_list2 <- combine_list
  combine_list3 <- combine_list
  
  #combine_list1[[length(path_list)+1]] <- id_403
  # combine_list2[[length(path_list)+1]] <- id_404
  combine_list3[[length(path_list)+1]] <- id_403
  combine_list3[[length(path_list)+2]] <- id_404
  
  #FAERS only
  signal_intersect <- Reduce(intersect,combine_list)
  
  #signal_intersect1 <- Reduce(intersect,combine_list1)
  # signal_intersect2 <- Reduce(intersect,combine_list2)
  signal_intersect3 <- Reduce(intersect,combine_list3)
  
  #write.csv(as.data.frame(signal_intersect3),"signal_intersect.csv",row.names=FALSE)
  
  signals_list <- list()
  vars_list <- list()
  
  #signals_list1 <- list()
  #vars_list1 <- list()
  # signals_list2 <- list()
  # vars_list2 <- list()
  signals_list3 <- list()
  vars_list3 <- list()
  
  for(i in 1:length(path_list)){
    path <- paste(path_base,path_list[i], sep ="", collapse = " ")
    return_list <- get_signal_info(path)
    ref <- return_list$ref_base
    score <- return_list$score_base
    var <- return_list$var_base
    id <- return_list$id
    # ref <- return_list$ref
    # score <- return_list$score
    # var <- return_list$var
    # id <- return_list$id
    
    score0 <- score[id %in% signal_intersect]
    var0 <- var[id %in% signal_intersect]
    signals_list[[i]] <- score0
    vars_list[[i]] <- var0
    
    # score1 <- score[id %in% signal_intersect1]
    # var1 <- var[id %in% signal_intersect1]
    # signals_list1[[i]] <- score1
    # vars_list1[[i]] <- var1
    # 
    # score2 <- score[id %in% signal_intersect2]
    # var2 <- var[id %in% signal_intersect2]
    # signals_list2[[i]] <- score2
    # vars_list2[[i]] <- var2
    
    score3 <- score[id %in% signal_intersect3]
    var3 <- var[id %in% signal_intersect3]
    signals_list3[[i]] <- score3
    vars_list3[[i]] <- var3
    
    #name <- paste("Score_",i,".csv",sep="")
    #write.csv(as.data.frame(score1),name,row.names=FALSE)
    
  }
  ref0 <- ref[id %in% signal_intersect]
  #ref1 <- ref[id %in% signal_intersect1]
  # ref2 <- ref[id %in% signal_intersect2]
  ref3 <- ref[id %in% signal_intersect3]
  
  
  #score_403_1 <- score_403[id_403 %in% signal_intersect1]
  #var_403_1 <- var_403[id_403 %in% signal_intersect1]
  
  score_403_3 <- score_403[id_403 %in% signal_intersect3]
  var_403_3 <- var_403[id_403 %in% signal_intersect3]
  
  # score_404_2 <- score_404[id_404 %in% signal_intersect2]
  # var_404_2 <- var_404[id_404 %in% signal_intersect2]
  
  score_404_3 <- score_404[id_404 %in% signal_intersect3]
  var_404_3 <- var_404[id_404 %in% signal_intersect3]
  
  #write.csv(as.data.frame(score_403_1),"Claimes_403_score.csv",row.names=FALSE)
  
  
  
  #Signal combination Claims + FAERS
  #claims_signal_list1 <- list(score_403_1)
  #claims_var_list1 <- list(var_403_1)
  # claims_signal_list2 <- list(score_404_2)
  # claims_var_list2 <- list(var_404_2)
  claims_signal_list3 <- list(score_403_3, score_404_3)
  claims_var_list3 <- list(var_403_3, var_404_3)
  
  print("CLAIMS + FAERS")
  #combined_score_ours1 <- combine_signals_ours(signals_list1, vars_list1, claims_signal_list1, claims_var_list1)
  print("ours")
  #print("Faers + 403")
  #print(roc(ref1, combined_score_ours1, percent=TRUE, plot=TRUE))

  # combined_score_ours2 <- combine_signals_ours(signals_list2, vars_list2, claims_signal_list2, claims_var_list2)
  # print("Faers + 404")
  # print(roc(ref2, combined_score_ours2, percent=TRUE, plot=TRUE))
  
  combined_score_ours3 <- combine_signals_ours(signals_list3, vars_list3, claims_signal_list3, claims_var_list3)
  #print("Faers + 404 + 403")
  #print(roc(ref3, combined_score_ours3, percent=TRUE, plot=TRUE))
  #write.csv(as.data.frame(combined_score_ours3),"Combine_score3.csv",row.names=FALSE)
  #signals_list1[[length(path_list)+1]] <- score_403_1
  #vars_list1[[length(path_list)+1]] <- var_403_1
  # 
  # signals_list2[[length(path_list)+1]] <- score_404_2
  # vars_list2[[length(path_list)+1]] <- var_404_2
  
  combined_score_kdd13_avg <- combine_signals_kdd13(signals_list3, vars_list3,1)
  write.csv(as.data.frame(combined_score_kdd13_avg),"Combine_score_faersonly_base.csv",row.names=FALSE)
  
  combined_score_kdd13_avg <- combine_signals_kdd13(claims_signal_list3, claims_var_list3,1)
  write.csv(as.data.frame(combined_score_kdd13_avg),"Combine_score_claimsonly_base.csv",row.names=FALSE)
  
  #combined_score_kdd13_avg <- combine_signals_kdd13(signals_list3, vars_list3,1)
  # combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list3, vars_list3,0)
  #print("Faers + 403 + 404")
  # print("kdd13_avg")
 # print(roc(ref3, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  #write.csv(as.data.frame(combined_score_kdd13_avg),"Combine_score_faersonly_base.csv",row.names=FALSE)
  
  #claims_score_kdd13_avg <- combine_signals_kdd13(claims_signal_list3, claims_var_list3,1)
  #print(roc(ref3, claims_score_kdd13_avg, percent=TRUE, plot=TRUE))
  #write.csv(as.data.frame(claims_score_kdd13_avg),"Combine_score_claimsonly_base.csv",row.names=FALSE)

  # signals_list3[[length(path_list)+1]] <- score_403_3
  # signals_list3[[length(path_list)+2]] <- score_404_3
  # vars_list3[[length(path_list)+1]] <- var_403_3
  # vars_list3[[length(path_list)+2]] <- var_404_3


  #combined_score_kdd13_avg <- combine_signals_kdd13(signals_list1, vars_list1,1)
  # combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list1, vars_list1,0)
  #print("Faers + 403")
  #print("kdd13_avg")
  #print(roc(ref1, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  # print("kdd13_wavg")
  # print(roc(ref1, combined_score_kdd13_wavg, percent=TRUE, plot=TRUE))
  # 
  # combined_score_kdd13_avg <- combine_signals_kdd13(signals_list2, vars_list2,1)
  # combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list2, vars_list2,0)
  # print("Faers + 404")
  # print("kdd13_avg")
  # print(roc(ref2, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  # print("kdd13_wavg")
  # print(roc(ref2, combined_score_kdd13_wavg, percent=TRUE, plot=TRUE))

 #combined_score_kdd13_avg <- combine_signals_kdd13(signals_list3, vars_list3,1)
  #write.csv(as.data.frame(combined_score_kdd13_avg),"kdd13.csv",row.names=FALSE)
  # combined_score_kdd13_wavg <- combine_signals_kdd13(signals_list3, vars_list3,0)
 # print("Faers + 403 + 404")
  # print("kdd13_avg")
 # print(roc(ref3, combined_score_kdd13_avg, percent=TRUE, plot=TRUE))
  # print("kdd13_wavg")
  # print(roc(ref3, combined_score_kdd13_wavg, percent=TRUE, plot=TRUE))
  # 
  # print("ONLY FAERS")
  # only_kdd13_avg <- combine_signals_kdd13(signals_list, vars_list,1)
  # only_kdd13_wavg <- combine_signals_kdd13(signals_list, vars_list,0)
  # print("kdd13_avg")
  # print(roc(ref0, only_kdd13_avg, percent=TRUE, plot=TRUE))
  # print("kdd13_wavg")
  # print(roc(ref0, only_kdd13_wavg, percent=TRUE, plot=TRUE))
  
}



path_list <- c(2007,2008,2009,2010,2011,2012,2013,2014)
claims_path <- "/Users/inci.baytasibm.com/Documents"
get_auc_combined(path_list, claims_path)

