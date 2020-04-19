#Packages
library(caret)
library(xgboost)

##data##
o_indx = c(rep(c(rep(T,4),rep(T,4)),3),
           rep(c(rep(T,4),rep(T,4)),3),
           rep(c(rep(T,4),rep(T,4)),3),
           rep(c(rep(T,4),rep(T,4)),3))


nn = sum(o_indx)

chk = apply(alg_dat!=0,2,sum)    
slt_dat = alg_dat[,chk>1]          
dim(slt_dat)
sum(chk>1)                         
std_dat = slt_dat[o_indx,]         
dim(std_dat)
resp = resp[o_indx]

#----------------------------------------------------------------------------
#Initializing the vecs & mats

dic.acc <- c()
idx_mat <- matrix(0,24,10)
m <- matrix(0,24,9)
prob <- rep(list(m),10)
est_cl_ftn  <- function(x) return(names(which.max(table(x))))
iter<-1

performance <- function(tb){
  sensitivity <- tb[2,2]/sum(tb[2,])
  specificity <- tb[1,1]/sum(tb[1,])
  FDR <- tb[1,2]/(tb[2,2]+tb[1,2])
  Accuracy <- (tb[1,1]+tb[2,2])/(tb[1,1]+tb[1,2]+tb[2,1]+tb[2,2])
  F1.score <- tb[2,2]/(2*tb[2,2]+tb[1,2]+tb[2,1]) 
  
  return(list(sensitivity=sensitivity, specificity=specificity, FDR=FDR, Accuracy=Accuracy, F1.score=F1.score))
}


#xgboost


while(iter<11){
  
  set.seed(iter)
  
  lab_idx1 <- sample(which(resp=="1"), 6)
  lab_idx2 <- sample(which(resp=="2"), 6)
  lab_idx3 <- sample(which(resp=="3"), 6)
  lab_idx4 <- sample(which(resp=="4"), 6)
  lab_idx  <- c(lab_idx1,lab_idx2,lab_idx3,lab_idx4)
  
  idx_mat[,iter]<- lab_idx
  
  
  ##build XGB
  for(i in 1:9){
    
    set.seed(iter*1000+i)
    
    #3?? class???? 6???? tissue???? 6?? ?? xdata & label
    std_dat2 <- std_dat[-lab_idx,]
    resp2    <- as.numeric(resp[-lab_idx])
    
    #36x18 binary class data
    idx123 <- sample(which(resp2!=4), 36, replace = T)
    idx4   <- which(resp2==4)
    idx    <- append(idx123,idx4)
    
    std_dat3 <- std_dat2[idx,]
    resp3    <- as.numeric(resp2[idx])
    resp3[which(resp3!=4)] <- 0
    resp3[which(resp3==4)] <- 1
    
    train_dt <- xgb.DMatrix(data = std_dat3, label=resp3)

    #test data
    test_xdt <- as.data.frame(std_dat[lab_idx,])
    
    test_ydt <- as.numeric(resp[lab_idx])
    test_ydt[which(test_ydt!=4)] <- 0
    test_ydt[which(test_ydt==4)] <- 1
    
    test_dt <- xgb.DMatrix(data=std_dat[lab_idx,], label=test_ydt)
    
    ##model fitting
    params1 <- list(booster = "gbtree", 
                   objective = "binary:logistic", 
                   eta=0.1, 
                   gamma=0, 
                   max_depth=3, 
                   min_child_weight=1, 
                   subsample=1,
                   colsample_bytree=1,
                   lambda = 0,
                   eval_metric = "error",
                   silent = 1,
                   alpha=0.1
                   )
    
    params2 <- list(booster = "dart", 
                    objective = "binary:logistic", 
                    eta=0.01, 
                    gamma=0, 
                    max_depth=4, 
                    sample_type='uniform',
                    min_child_weight = 1,
                    normalize_type = 'tree',
                    rate_drop=0.7,
                    skip_drop=0.5,
                    colsample_bytree=1,
                    eval_metric="error"
                   )

   xgb.cv <- xgb.cv(params = params1,
                    data = train_dt,
                    nrounds = 300,
                    nfold = 30,
                    showsd = TRUE,
                    stratified = TRUE,
                    print_every_n = 50,
                    early_stopping_rounds = 350,
                    maximize = F # error:F, auc:T
                    )
   print(xgb.cv)   
   xgb.fit <- xgboost(params = params1,
                      data = train_dt,
                      nrounds = xgb.cv$best_iteration
                      #watchlist = list(val=test_dt, train=train_dt)
                      )
  
    predict(xgb.fit, data.matrix(test_xdt), type='response')
    
    ##calculate probabilties
    prob[[iter]][,i] <- predict(xgb.fit, data.matrix(test_xdt), type='response')
    
    
  }
  
  prd_result     <- as.numeric(apply(ifelse(prob[[iter]]>=0.75,1,0), 1, est_cl_ftn))
  dic.acc[iter]  <- mean(test_ydt==prd_result)
  
  
  cat("***** End of ",iter,"th iteration*****\n")
  iter <- iter+1
  
}


#----------------------------------------------------------------------------
#Test result
tb.list <- rep(list(0),10)
ans <- c(rep(0,18), rep(1,6))
for(i in 1:10){
  p <- as.numeric(apply(ifelse(prob[[i]]>=0.5,1,0), 1, est_cl_ftn))
  tb.list[[i]] <- table(ans,p)
}

tb.list

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#Test result

result <- unlist(sapply(tb.list, performance))
mean(dic.acc)

res <- matrix(result, 5, 10, 
              dimnames = list(c("sensitivity","specificity", "FDR", "Accuracy", "F1.score"), 
                              c(1:10)))

apply(res,1,mean)

