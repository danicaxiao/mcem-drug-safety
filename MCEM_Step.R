setwd("/home/cxiao/workspace")
rm(list = ls())

# function as.PhViD
as.PhViD = function (DATA.FRAME, MARGIN.THRES = 1) 
{
  data <- DATA.FRAME
  data[, 1] <- as.factor(DATA.FRAME[, 1])
  data[, 2] <- as.factor(DATA.FRAME[, 2])
  data[, 3] <- as.double(DATA.FRAME[, 3])
  coln <- names(data)
  names(data)[3] <- "n11"
  data_cont <- xtabs(n11 ~ ., data = data)
  n1._mat <- apply(data_cont, 1, sum)
  n.1_mat <- apply(data_cont, 2, sum)
  if (MARGIN.THRES > 1) {
    while (sum(n1._mat < MARGIN.THRES) > 0 | sum(n.1_mat < 
                                                 MARGIN.THRES) > 0) {
      data_cont <- data_cont[n1._mat >= MARGIN.THRES, ]
      data_cont <- data_cont[, n.1_mat >= MARGIN.THRES]
      n1._mat <- apply(data_cont, 1, sum)
      n.1_mat <- apply(data_cont, 2, sum)
    }
  }
  coord <- which(data_cont != 0, arr.ind = TRUE)
  coord <- coord[order(coord[, 1]), ]
  Nb_n1. <- length(n1._mat)
  Nb_n.1 <- length(n.1_mat)
  libel.medoc <- rownames(data_cont)[coord[, 1]]
  libel.effet <- colnames(data_cont)[coord[, 2]]
  n11 <- data_cont[coord]
  N <- sum(n11)
  n1. <- n1._mat[coord[, 1]]
  n.1 <- n.1_mat[coord[, 2]]
  RES <- vector(mode = "list")
  RES$L <- data.frame(libel.medoc, libel.effet)
  colnames(RES$L) <- coln[1:2]
  RES$data <- cbind(n11, n1., n.1)
  rownames(RES$data) <- paste(libel.medoc, libel.effet)
  RES$N <- N
  return(RES)
}

##formula 12 in the table
.lik2NB <- function(p,n11,E){
  sum(-log((p[5] * dnbinom(n11, size=p[1], prob=p[2]/(p[2]+E)) + (1-p[5]) * dnbinom(n11, size=p[3], prob=p[4]/(p[4]+E)))))
}

##likelihood function##
.likTronc2NB <- function(p, n11, E, tronc){
  nb1 <- dnbinom(n11, size = p[1], prob = p[2] / (p[2] + E)) /
    (1 - pnbinom(tronc, size = p[1], prob = p[2] / (p[2] + E)))
  
  nb2 <- dnbinom(n11, size = p[3], prob = p[4] / (p[4] + E)) /
    (1 - pnbinom(tronc, size = p[3], prob = p[4] / (p[4] + E)))
  L <- (p[5] * nb1) + ((1 - p[5]) * nb2)
  sum(-log(L))
}


# LB <- .QuantileDuMouchel(0.05, Q, PRIOR.PARAM[1] + n11, PRIOR.PARAM[2] + 
#                             E, PRIOR.PARAM[3] + n11, PRIOR.PARAM[4] + E)
#Seuil <- 0.05
#a1 <- PRIOR.PARAM[1] + n11
#b1 <- PRIOR.PARAM[2] + E
#a2 <- PRIOR.PARAM[3] + n11
#b2 <- PRIOR.PARAM[4] + E
.QuantileDuMouchel<-function(Seuil,Q,a1,b1,a2,b2) {
  m<-rep(-100000,length(Q))
  M<-rep(100000,length(Q))
  x<-rep(1,length(Q))
  Cout<- .FCoutQuantileDuMouchel(x,Seuil,Q,a1,b1,a2,b2)
  while (max(round(Cout*1e4))!=0)	{
    S<-sign(Cout)
    xnew<-(1+S)/2*((x+m)/2)+(1-S)/2*((M+x)/2) ##??? what does it do
    M<-(1+S)/2*x+(1-S)/2*M
    m<-(1+S)/2*m+(1-S)/2*x
    x<-xnew
    Cout<-.FCoutQuantileDuMouchel(x,Seuil,Q,a1,b1,a2,b2)
  }
  x
}
.FCoutQuantileDuMouchel<-function(p,Seuil,Q,a1,b1,a2,b2) {
  Q*pgamma(p,shape=a1,rate=b1)+(1-Q)*pgamma(p,shape=a2,rate=b2)-Seuil
}


# function GPS
GPS = function (DATABASE, RR0 = 1, MIN.n11 = 1, DECISION = 1, DECISION.THRES = 0.05, 
                RANKSTAT = 1, TRONC = FALSE, TRONC.THRES = 1, PRIOR.INIT = c(alpha1 = 0.2, beta1 = 0.06, alpha2 = 1.4, beta2 = 1.8, w = 0.1), PRIOR.PARAM = NULL)
{
  DATA <- DATABASE$data
  N <- DATABASE$N
  L <- DATABASE$L
  n11 <- DATA[, 1]
  n1. <- DATA[, 2]
  n.1 <- DATA[, 3]
  E <- DATA[, 2] * DATA[, 3]/N
  P_OUT <- TRUE
  if (is.null(PRIOR.PARAM)) {
    P_OUT <- FALSE
    if (TRONC == FALSE) {
      data_cont <- xtabs(DATA[, 1] ~ L[, 1] + L[, 2])
      n1._mat <- apply(data_cont, 1, sum)
      n.1_mat <- apply(data_cont, 2, sum)
      n1._c <- rep(n1._mat, times = length(n.1_mat))
      n.1_c <- rep(n.1_mat, each = length(n1._mat))
      E_c <- n1._c * n.1_c/N
      n11_c <- as.vector(data_cont)
      p_out <- suppressWarnings(nlm(.lik2NB, p = PRIOR.INIT, 
                                    n11 = n11_c, E = E_c, iterlim = 500))
    }
    if (TRONC == TRUE) {
      tronc <- TRONC.THRES - 1
      p_out <- suppressWarnings(nlm(.likTronc2NB, p = PRIOR.INIT, 
                                    n11 = n11[n11 >= TRONC.THRES], E = E[n11 >= TRONC.THRES], 
                                    tronc, iterlim = 500))
    }
    PRIOR.PARAM <- p_out$estimate 
    code.convergence <- p_out$code
  }
  if (MIN.n11 > 1) {
    E <- E[n11 >= MIN.n11]
    n1. <- n1.[n11 >= MIN.n11]
    n.1 <- n.1[n11 >= MIN.n11]
    LL <- data.frame(drugs = L[, 1], events = L[, 2], n11)
    LL1 <- LL[, 1][n11 >= MIN.n11]
    LL2 <- LL[, 2][n11 >= MIN.n11]
    rm(list = "L")
    L <- data.frame(LL1, LL2)
    n11 <- n11[n11 >= MIN.n11]
  }
  Nb.Cell <- length(n11)
  post.H0 <- vector(length = Nb.Cell)
  Q <- PRIOR.PARAM[5] * dnbinom(n11, size = PRIOR.PARAM[1], 
                                prob = PRIOR.PARAM[2]/(PRIOR.PARAM[2] + E))/(PRIOR.PARAM[5] * 
                                                                               dnbinom(n11, size = PRIOR.PARAM[1], prob = PRIOR.PARAM[2]/(PRIOR.PARAM[2] + 
                                                                                                                                            E)) + (1 - PRIOR.PARAM[5]) * dnbinom(n11, size = PRIOR.PARAM[3], 
                                                                                                                                                                                 prob = PRIOR.PARAM[4]/(PRIOR.PARAM[4] + E)))
  post.H0 <- Q * pgamma(RR0, PRIOR.PARAM[1] + n11, PRIOR.PARAM[2] + 
                          E) + (1 - Q) * pgamma(RR0, PRIOR.PARAM[3] + n11, PRIOR.PARAM[4] + 
                                                  E)
  # Posterior Expectation of log2(lambda) formula 10 in the paper
  postE <- log(2)^(-1) * (Q * (digamma(PRIOR.PARAM[1] + n11) - 
                                 log(PRIOR.PARAM[2] + E)) + (1 - Q) * (digamma(PRIOR.PARAM[3] + 
                                                                                 n11) - log(PRIOR.PARAM[4] + E)))
  LB <- .QuantileDuMouchel(0.05, Q, PRIOR.PARAM[1] + n11, PRIOR.PARAM[2] + 
                             E, PRIOR.PARAM[3] + n11, PRIOR.PARAM[4] + E)
  
  alpha1 <- PRIOR.PARAM[1] + n11
  beta1 <- PRIOR.PARAM[2] + E
  alpha2 <- PRIOR.PARAM[3] + n11
  beta2 <- PRIOR.PARAM[4] + E
  var <- Q*(1-Q)*(alpha1/beta1-alpha2/beta2)^2+Q*alpha1/(beta1^2)+(1-Q)*alpha2/(beta2^2)
  var <- var*1.0/N # s.e.
  
  if (RANKSTAT == 1) 
    RankStat <- post.H0
  if (RANKSTAT == 2) 
    RankStat <- LB
  if (RANKSTAT == 3) 
    RankStat <- postE
  if (RANKSTAT == 1) {
    FDR <- (cumsum(post.H0[order(RankStat)])/(1:length(post.H0)))
    FNR <- rev(cumsum((1 - post.H0)[order(1 - RankStat)]))/(Nb.Cell - 
                                                              1:length(post.H0))
    Se <- cumsum((1 - post.H0)[order(RankStat)])/(sum(1 - 
                                                        post.H0))
    Sp <- rev(cumsum(post.H0[order(1 - RankStat)]))/(Nb.Cell - 
                                                       sum(1 - post.H0))
  }
  if (RANKSTAT == 2 | RANKSTAT == 3) {
    FDR <- (cumsum(post.H0[order(RankStat, decreasing = TRUE)])/(1:length(post.H0)))
    FNR <- rev(cumsum((1 - post.H0)[order(1 - RankStat, decreasing = TRUE)]))/(Nb.Cell - 
                                                                                 1:length(post.H0))
    Se <- cumsum((1 - post.H0)[order(RankStat, decreasing = TRUE)])/(sum(1 - 
                                                                           post.H0))
    Sp <- rev(cumsum(post.H0[order(1 - RankStat, decreasing = TRUE)]))/(Nb.Cell - 
                                                                          sum(1 - post.H0))
  }
  if (DECISION == 1) 
    Nb.signaux <- sum(FDR <= DECISION.THRES)
  if (DECISION == 2) 
    Nb.signaux <- min(DECISION.THRES, Nb.Cell)
  if (DECISION == 3) {
    if (RANKSTAT == 1) 
      Nb.signaux <- sum(RankStat <= DECISION.THRES, na.rm = TRUE)
    if (RANKSTAT == 2 | RANKSTAT == 3) 
      Nb.signaux <- sum(RankStat >= DECISION.THRES, na.rm = TRUE)
  }
  Q_func <- Q
  Expect_Q <- postE
  RES <- vector(mode = "list")
  RES$INPUT.PARAM <- data.frame(RR0, MIN.n11, DECISION, DECISION.THRES, 
                                RANKSTAT, TRONC, TRONC.THRES)
  RES$STOPPING <- data.frame(Q_func, Expect_Q)
  RES$PARAM <- vector(mode = "list")
  if (P_OUT == TRUE) 
    RES$PARAM$PRIOR.PARAM <- data.frame(PRIOR.PARAM)
  if (P_OUT == FALSE) {
    RES$PARAM$PRIOR.INIT <- data.frame(PRIOR.INIT)
    RES$PARAM$PRIOR.PARAM <- PRIOR.PARAM
    RES$PARAM$CONVERGENCE <- code.convergence
  }
  if (RANKSTAT == 1) {
    RES$ALLSIGNALS <- data.frame(L[, 1][order(RankStat)], 
                                 L[, 2][order(RankStat)], n11[order(RankStat)], E[order(RankStat)], 
                                 RankStat[order(RankStat)], (n11/E)[order(RankStat)], 
                                 n1.[order(RankStat)], n.1[order(RankStat)], FDR, 
                                 FNR, Se, Sp,var)
    colnames(RES$ALLSIGNALS) <- c("drug", "event", "count", 
                                  "expected count", "postH0", "n11/E", "drug margin", 
                                  "event margin", "FDR", "FNR", "Se", "Sp","var")
  }
  if (RANKSTAT == 2 | RANKSTAT == 3) {
    RES$ALLSIGNALS <- data.frame(L[, 1][order(RankStat, decreasing = TRUE)], 
                                 L[, 2][order(RankStat, decreasing = TRUE)], n11[order(RankStat, 
                                                                                       decreasing = TRUE)], E[order(RankStat, decreasing = TRUE)], 
                                 RankStat[order(RankStat, decreasing = TRUE)], (n11/E)[order(RankStat, 
                                                                                             decreasing = TRUE)], n1.[order(RankStat, decreasing = TRUE)], 
                                 n.1[order(RankStat, decreasing = TRUE)], FDR, FNR, 
                                 Se, Sp, post.H0[order(RankStat, decreasing = TRUE)], var)
    if (RANKSTAT == 2) 
      colnames(RES$ALLSIGNALS) <- c("drug", "event", "count", 
                                    "expected count", "Q_0.05(lambda)", "n11/E", 
                                    "drug margin", "event margin", "FDR", "FNR", 
                                    "Se", "Sp", "postH0", "var")
    if (RANKSTAT == 3) 
      colnames(RES$ALLSIGNALS) <- c("drug", "event", "count", 
                                    "expected count", "post E(Lambda)", "n11/E", 
                                    "drug margin", "event margin", "FDR", "FNR", 
                                    "Se", "Sp", "postH0", "var")
  }
  RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux, ]
  RES$NB.SIGNALS <- Nb.signaux
  return(RES)
}


InitMgps <- function(sample, baseline.dump) {
  if (file_test("-f", baseline.dump)) {
    return(FALSE)
  }
  
  # calculate contengency table
  testTwoWayTable = base::table(sample$drug_id,sample$adr_id) # side effect of spark
  
  testDataFrame = as.data.frame(testTwoWayTable)
  colnames(testDataFrame) = c("drug_id","adr_id","freq")
  testDataFrame1 = testDataFrame[testDataFrame$freq>0,]
  
  
  # write.table(testDataFrame1,"omop_adr_related_random_contigency_table.tsv",row.names = FALSE,sep="\t")
  # load contengency table 
  # testDataFrame1 = read.table("omop_adr_related_random_contigency_table.tsv",sep="\t",comment="",header=TRUE)
  
  
  # for calculating disproportionality score
  PhViDdata <- as.PhViD(testDataFrame1)
  mgps <- GPS(PhViDdata, DECISION = 3, DECISION.THRES = 2, RANKSTAT = 2)
  # = 5% quantile of the posterior distribution of lambda
  mgps$ALLSIGNALS[1:5, c('drug', 'event', 'Q_0.05(lambda)', 'Se')]
  
  save(sample, mgps, file = baseline.dump)
  
  return(TRUE) 
}

faers.omop.random.report <- "faers2014_1on1.tsv"
#"faers2011.tsv"
# load raw data
#raw.data = read.table(faers.omop.random.report,header=TRUE, sep="\t", comment="")
#colnames(raw.data) = c("report_id","drug_id","role_cod","adr_id","report_year","age","gender","occu_cod")
#unique.data = unique(raw.data[,c("report_id","drug_id","adr_id")])
unique.data = read.table(faers.omop.random.report,header=TRUE, sep="\t", comment="")
unique.data = as.data.frame(unique.data)

baseline.faers.omop.dump <- "faers.RData"
# initialize risk score
InitMgps(unique.data, baseline.faers.omop.dump) 

load(baseline.faers.omop.dump)

#save.image("./unittestdump.RData")

#######################################################################
# the MCEM procedure 

library(dplyr)

# subset sample and prepare dict
# sample <- unique.data[1:800,]
sample <- unique.data
dict <- mgps$ALLSIGNALS[, c('Q_0.05(lambda)', 'Se')]  # row key: paste(1367268, 35104756)

# TODO try what stat to use
#rank.stat <- 'post E(Lambda)'
rank.stat <- 'Q_0.05(lambda)'
colnames(dict) = c(rank.stat, 'Se')


na.zero <- function (x) {
  x[is.na(x) | (x < 1e-2)] <- 1e-2
  return(x)
}


McStep <- function(sample, dict) {
  # lookup by dict, signal goes to 'lambda'
  sample['index'] <- paste(sample$drug_id, sample$adr_id)
  result <- sample %>% mutate(lambda = na.zero(dict[index, rank.stat]), Se = dict[index, 'Se'])
  
  # determine primary cause by multinom
  # selected set is a subset of overall sample so that drug ADR pairs in the omop reference
  # is always decreased
  selected <- result %>% 
    group_by(report_id, adr_id) %>% 
    summarise(drug_id = drug_id[which(rmultinom(1, 1, lambda)[,1]>0, arr.ind=TRUE)])
  
  return(as.data.frame(selected)[, c('report_id', 'drug_id', 'adr_id')])
}

MaxStep <- function(sample, prior) {
  # calculate contengency table
  testTwoWayTable = base::table(sample$drug_id,sample$adr_id)
  
  testDataFrame = as.data.frame(testTwoWayTable)
  colnames(testDataFrame) = c("drug_id","adr_id","freq")
  testDataFrame1 = testDataFrame[testDataFrame$freq>0,]
  
  PhViDdata <- as.PhViD(testDataFrame1)
  mgps <- GPS(PhViDdata, DECISION = 3, DECISION.THRES = 2, RANKSTAT = 2, PRIOR.INIT = prior)
  # RANKSTAT = 3: the posterior distribution of lambda
  #prior <- mgps$PARAM$PRIOR.PARAM
  #var <- get.lambda.Variance(PhViDdata, prior)
  
  print(mgps$ALLSIGNALS[1:5, c('drug', 'event', rank.stat, 'Se')])
  
  #return_list <- list("score"=mgps, "variance"=var)
  
  return(mgps)
}


# iteration 1
next.sample <- McStep(sample, dict)
# save sampling results 
write.csv(as.data.frame(next.sample),"sample1.csv",row.names=FALSE)
next.sample <- rbind(sample, next.sample)
next.mgps <- MaxStep(next.sample, mgps$PARAM$PRIOR.PARAM)
#next.mgps <- return_list$score
#next.var <- return_list$variance


print("First iteration is done")

save(next.sample,next.mgps, file = paste(1,".RData",sep=""))

# continue to iterate
max.iters <- 500
for (i in 2:max.iters) {
  message(sprintf("Iteration: %d\n", i))
  prev.mgps = next.mgps
  dict <- next.mgps$ALLSIGNALS[, c(rank.stat, 'Se')]
  #Running time
  ptm <- proc.time()
  prev.sample <- next.sample
  next.sample <- McStep(sample, dict)
  write.csv(as.data.frame(next.sample),paste("sample",i,".csv",sep=""),row.names=FALSE)
  next.sample <- rbind(prev.sample, next.sample)
  print(proc.time() - ptm)
  
  #next.mgps <- MaxStep(next.sample, mgps$PARAM$PRIOR.PARAM)  # TODO can also be next.mgps$PARAM$PRIOR.PARAM
  
  #error handle
  ptm <- proc.time()
  tryCatch({
    next.mgps <- MaxStep(next.sample, mgps$PARAM$PRIOR.PARAM) # TODO can also be next.mgps$PARAM$PRIOR.PARAM
    #next.mgps <- return_list$score
    #next.var <- return_list$variance
  }, warning = function(war) {
    # warning handler picks up where error was generated
    print('warning')
  }, error = function(err) {
    # error handler picks up where error was generated
    print('error')
  }, finally = {
  }) # END tryCatch
  print(proc.time() - ptm)
  
  if (i %% 1 == 0) {
    # TODO evaluate accuracy
    save(next.sample,next.mgps, file = paste(i,".RData",sep=""))
    Q.current <- next.mgps$STOPPING$Q_func
    Q.prev <- prev.mgps$STOPPING$Q_func
    diff <- mapply('-', Q.current, Q.prev)
    diff <- mapply('abs', diff)
    val <- Reduce("+", diff)/length(diff)
    message(sprintf("Difference: %f\n", val))
    if (val < 1e-3){
      print('Algorithm converged!')
      break
    }
  }
}
