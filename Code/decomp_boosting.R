library(glmnet)

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


#glmnet parameter
grid_lam = c(seq(0.001,0.01,by=0.001),seq(0.01,2,by=0.01)) #glmnet lambda
grid_al = seq(0,1,by=0.1) 
dic.lam <-matrix(0,10,10)
dic.al  <-matrix(0,10,10)
dic.acc <- c()

idx_mat <- matrix(0,12,4)
m <- matrix(0,12,11)
prob <- list(m,m,m,m)
est_cl_ftn  <- function(x) return(names(which.max(table(x))))


#beta estimator
b0 <-c()
bs <- list()


#validation value
a <- matrix(0,10,10)
error <- matrix(0,10,10)
w0 <- rep(1/9,9)
weight <- matrix(0,9,11)
weight[,1] <- w0
weight <- rep(list(weight), 10)
mis_prd <- rep(list(matrix(0,9,10)),10)

prd <- matrix(0,9,10)
prob_te <- matrix(0,24,10)
pp <- matrix(0,24,10)
tb <- list()



#performance measure
accuracy <- rep(0,10)
dff <- matrix(0, 4, 10, 
              dimnames = list(c("sensitivity","specificity", "FDR",  "F1.score"), 
                              c(1:10)))

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
#glmnet + boosting

iter <- 1
while(iter<11){
  
  set.seed(iter)
  
  lab_idx1 <- sample(which(resp=="1"), 6)
  lab_idx2 <- sample(which(resp=="2"), 6)
  lab_idx3 <- sample(which(resp=="3"), 6)
  lab_idx4 <- sample(which(resp=="4"), 6)
  lab_idx  <- c(lab_idx1,lab_idx2,lab_idx3,lab_idx4)
  
  
  for(i in 1:10){
    
    set.seed(1000*iter+i)
    #?? class???? 6?? ?? xdata & label
    std_dat2 <- std_dat[-lab_idx,]
    resp2    <- as.numeric(resp[-lab_idx])
  
    #(30+15) train data, (3+3) validataion data
    idx123 <- sample(which(resp2!=4), 39, replace=T) #3class???? ?????ϰ? 39 ????
    idx123_val <- idx123[1:6] #6?? val
    idx123_tr  <- idx123[7:length(idx123)] #36?? train
  
    idx4   <- which(resp2==4)
    idx4_val <- idx4[1:3]
    idx4_tr  <- idx4[4:length(idx4)]
  
    idx_val <- append(idx123_val, idx4_val)
    idx <- append(idx123_tr,idx4_tr)
  
    std_dat3 <- std_dat2[idx,]
    resp3    <- as.numeric(resp2[idx])
    resp3[which(resp3!=4)] <- 0
    resp3[which(resp3==4)] <- 1
  
  
    cv_mat = matrix(0,length(grid_lam),length(grid_al))
    colnames(cv_mat) = grid_al
    rownames(cv_mat) = grid_lam[length(grid_lam):1] 
  
  
    #estimate best paramter for 10 models
    for(k in 1:length(grid_al))
    {
      alpha = grid_al[k]
      cvfit=cv.glmnet(std_dat3, 
                      resp3, 
                      family='binomial', 
                      alpha=alpha, 
                      lambda=grid_lam,
                      type.measure="deviance", 
                      nfold=nn, 
                      standardize=F)     
      
      cv_mat[,k] = cvfit$cvm             
      cat(k," - ")
    }
    cat("\n")
  
    #select lambda & alpha 
    lam <- as.numeric(names(which.min(apply(cv_mat,1,min))))
    al  <- grid_al[which.min(apply(cv_mat,2,min))]
    dic.lam[i,iter] <- lam
    dic.al[i,iter]  <- al
  
  
    ##model fitting
    test <- glmnet(std_dat3, resp3, family='binomial', alpha=al, lambda=lam, standardize=F)
    b0[i] <- test$a0
    bs[[i]] <- test$beta
  
  
    #validation data
    val_dat <- std_dat2[idx_val,]
    val_lab <- as.numeric(resp2[idx_val])
    val_lab[which(val_lab!=4)] <- 0
    val_lab[which(val_lab==4)] <- 1
  
    prd[,i] <- predict(test, newx=val_dat, type='response', s=lam, standardize=F)
    pr <- prd[,i]
    pr2 <- ifelse(pr>=0.75,1,0)
    pr2[pr2==0] <- 0
  
    mis_prd[[iter]][,i] <- as.numeric(val_lab!=pr2)
  
    ##error, a(evaluation of weak learner), z(normalize constant)
    error[i,iter] <- sum(weight[[iter]][,i] *mis_prd[[iter]][,i])
    a[i,iter] <- 0.5*log((1-error[i,iter])/(error[i,iter]+0.01))
  
  
    ##weight updating
    for(jj in 1:length(val_lab)){
      if(val_lab[jj]==pr2[jj]) 
        weight[[iter]][jj,i+1] <- weight[[iter]][jj,i]*exp(-a[i,iter])
      else 
        weight[[iter]][jj,i+1] <- weight[[iter]][jj,i]*exp(a[i,iter])
    }
  
    ##normalize
    weight[[iter]][,i+1] <- weight[[iter]][,i+1]/sum(weight[[iter]][,i+1])
  
  }


  #Strong classifier 
  out <-rep(0,324) 
  for(kk in 1:10){
    beta <- bs[[kk]]
    aa <- a[kk,iter]
    tmp <- as.vector(beta*aa)
    out <- out+tmp
  }  

  betas <- out
  beta0 <- sum(b0*a[,iter])



  new_dat <- std_dat[lab_idx,]
  resp_te    <- as.numeric(resp[lab_idx])
  resp_te[which(resp_te!=4)] <- 0
  resp_te[which(resp_te==4)] <- 1


  ##test probability
  for(ss in 1:24){
    new_data <- t(as.matrix(new_dat[ss,]))
    prob_te[ss,iter] <- as.numeric(exp(sign(beta0+new_data%*%betas))/(1+exp(sign(beta0+new_data%*%betas))))
  }

  ##test table
  pp[,iter] <- ifelse(prob_te[,iter]>=0.70,1,0)
  tb[[iter]] <- table(resp_te, pp[,iter])

  
  cat(iter, "th loop finished \n")
  iter <- iter+1
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------


tb 


pp
mis_prd
error
weight
a

dic.al
dic.lam



#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#Test result

result <- unlist(sapply(tb, performance))

res <- matrix(result, 5, 10, 
              dimnames = list(c("sensitivity","specificity", "FDR", "Accuracy", "F1.score"), 
                              c(1:10)))
 

aaaaa <- apply(res,1,mean)

aaaaa
