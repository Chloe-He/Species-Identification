##Preparing data
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

grid_lam = c(seq(0.001,0.01,by=0.001), seq(0.01,2,by=0.01)) #glmnet lambda
grid_al = seq(0,1,by=0.1) #glmnet alpha

dic.lam <-matrix(0,11,10)
dic.al  <-matrix(0,11,10)
dic.acc <- c()
idx_mat <- matrix(0,24,10)
m <- matrix(0,24,11)
prob <- rep(list(m),10)
iter<-1

est_cl_ftn  <- function(x) return(names(which.max(table(x))))
performance <- function(tb){
  sensitivity <- tb[2,2]/sum(tb[2,])
  specificity <- tb[1,1]/sum(tb[1,])
  FDR <- tb[1,2]/(tb[2,2]+tb[1,2])
  Accuracy <- (tb[1,1]+tb[2,2])/(tb[1,1]+tb[1,2]+tb[2,1]+tb[2,2])
  F1.score <- tb[2,2]/(2*tb[2,2]+tb[1,2]+tb[2,1]) 
  
  return(list(sensitivity=sensitivity, specificity=specificity, FDR=FDR, Accuracy=Accuracy, F1.score=F1.score))
}



#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#glmnet + bagging

while(iter<11){   #no. of weak learner : 11
  
  set.seed(iter)
  
  lab_idx1 <- sample(which(resp=="1"), 6)
  lab_idx2 <- sample(which(resp=="2"), 6)
  lab_idx3 <- sample(which(resp=="3"), 6)
  lab_idx4 <- sample(which(resp=="4"), 6)
  lab_idx  <- c(lab_idx1,lab_idx2,lab_idx3,lab_idx4)
  
  idx_mat[,iter]<- lab_idx
  
  
  ##build 11 models for bagging
  for(i in 1:11){
    
    set.seed(1000*iter+i)
    #3?? class???? 6???? tissue???? 6?? ?? xdata & label
    std_dat2 <- std_dat[-lab_idx,]
    resp2    <- as.numeric(resp[-lab_idx])
    
    #36x18 binary class data
    idx123 <- sample(which(resp2!=4), 36, replace=T) 
    idx4   <- which(resp2==4)
    idx    <- append(idx123,idx4)
  
    std_dat3 <- std_dat2[idx,]
    resp3    <- as.numeric(resp2[idx])
    resp3[which(resp3!=4)] <- 0
    resp3[which(resp3==4)] <- 1
  
  
    cv_mat = matrix(0,length(grid_lam),length(grid_al))
    colnames(cv_mat) = grid_al
    rownames(cv_mat) = grid_lam[length(grid_lam):1] 
  
    
    #estimate best paramter for 11 models
    for(k in 1:length(grid_al))
    {
      alpha = grid_al[k]
      cvfit=cv.glmnet(std_dat3,resp3,family='binomial',alpha=alpha,lambda=grid_lam,
                      type.measure="deviance",nfold=nn,standardize=F)               #nfold=nn => loocv
      cv_mat[,k] = cvfit$cvm                                                        #cv-error
      cat("(",i,"-", k,")th parameter searching iteration \n")
    }
  
    
    #select lambda & alpha 
    lam <- as.numeric(names(which.min(apply(cv_mat,1,min))))
    al  <- grid_al[which.min(apply(cv_mat,2,min))]
    dic.lam[i,iter] <- lam
    dic.al[i,iter]  <- al
    
    
    ##model fitting
    test <- glmnet(std_dat3,resp3,family='binomial', alpha=al, lambda=lam, standardize=F)
summary(test)
    b0 <- test$a0
    bs <- test$beta
    bs[which(bs==0)] <- 0
    
    #test data
    new_dat <- std_dat[lab_idx,]
    new_lab <- as.numeric(resp[lab_idx])
    new_lab[which(new_lab!=4)] <- 0
    new_lab[which(new_lab==4)] <- 1
  
    
    #calculate probabilties
    for(kk in 1:24)
    {
      new_data <- t(as.matrix(new_dat[kk,]))
      prob[[iter]][kk,i] <- as.numeric(exp(b0+new_data%*%bs)/(1+exp(b0+new_data%*%bs)))
    }

    
    
  }

  #voting
  prd_result     <- as.numeric(apply(ifelse(prob[[iter]]>=0.75,1,0), 1, est_cl_ftn)) #threshold : 0.75
  dic.acc[iter]  <- mean(new_lab==prd_result)
  
  cat("********** End of ",iter,"th iteration*********\n")
  iter <- iter+1
  
}


#----------------------------------------------------------------------------


dic.acc
dic.lam
dic.al
r1 <- as.numeric(apply(ifelse(prob[[1]]>=0.75,1,0), 1, est_cl_ftn))
r2 <- as.numeric(apply(ifelse(prob[[2]]>=0.75,1,0), 1, est_cl_ftn))
r10 <- as.numeric(apply(ifelse(prob[[10]]>=0.75,1,0), 1, est_cl_ftn))


#----------------------------------------------------------------------------
#Test result

tb.list <- rep(list(0),10)
ans <- c(rep(0,18), rep(1,6)) ##true species##
for(i in 1:10){
  p <- as.numeric(apply(ifelse(prob[[i]]>=0.75,1,0), 1, est_cl_ftn))
  tb.list[[i]] <- table(ans,p)
}
p <- as.numeric(apply(round(prob[[5]],0), 1, est_cl_ftn))
tb.list


result <- unlist(sapply(tb.list, performance))###from here not running###
mean(dic.acc)

res <- matrix(result, 5, 10, 
              dimnames = list(c("sensitivity","specificity", "FDR", "Accuracy", "F1.score"), 
                              c(1:10)))

apply(res,1,mean)
