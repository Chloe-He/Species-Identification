#Packages
library(caret)
library(tree)
library(rpart)

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
m <- matrix(0,24,300) ##300 trees will be built
prob <- rep(list(m),10)
np <- sqrt(ncol(std_dat))
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

#----------------------------------------------------------------------------
#random forest


while(iter<11)
{
  
  set.seed(iter)
  #test set randomly selected
  lab_idx1 <- sample(which(resp=="1"), 6)
  lab_idx2 <- sample(which(resp=="2"), 6)
  lab_idx3 <- sample(which(resp=="3"), 6)
  lab_idx4 <- sample(which(resp=="4"), 6)
  ## Extract the index number for the test observations selected
  lab_idx  <- c(lab_idx1,lab_idx2,lab_idx3,lab_idx4) 
  ## The test set
  idx_mat[,iter]<- lab_idx
  
  
  ##build 300 tree models for randomforest
  for(i in 1:300)
  {
    
    set.seed(1000*iter+i)
    
    ## The training set that will be used to build 300 trees 
    std_dat2 <- std_dat[-lab_idx,]
    resp2    <- as.numeric(resp[-lab_idx])
    
    #36x18 binary class data
    idx123 <- sample(which(resp2!=4), 36, replace=T) 
    idx4   <- which(resp2==4)
    idx    <- append(idx123,idx4)
    npp <- sample(1:ncol(std_dat), np, replace=F)
    
    std_dat3 <- std_dat2[idx,npp]
    resp3    <- as.numeric(resp2[idx])
    resp3[which(resp3!=4)] <- 0
    resp3[which(resp3==4)] <- 1
    
    levels <- unique(resp3) 
    resp3 <- factor(resp3, labels=make.names(levels))
    

    #train data
    train_boots_data <- as.data.frame(cbind(resp3, std_dat3))
    levels <- unique(train_boots_data$resp3) 
    train_boots_data$resp3 <- factor(train_boots_data$resp3, labels=make.names(levels))
    
    #test data
    test_xdt <- as.data.frame(std_dat[lab_idx,])
    
    test_ydt <- as.numeric(resp[lab_idx])
    test_ydt[which(test_ydt!=4)] <- 1
    test_ydt[which(test_ydt==4)] <- 2
    levels <- unique(test_ydt) 
    test_ydt <- factor(test_ydt, labels=make.names(levels))
    

    ##model fitting
    my.con <- list(mincut=1, minsize=2, nmax=100, mindev=0.01)
    tree <- tree(resp3~., data=train_boots_data, control=my.con)
    
    summary(tree)
    plot(tree); text(tree)
    

    #pruning using cv(parameter tuning)
    cv.trees <- cv.tree(tree, FUN=prune.misclass, K=76)
    plot(cv.trees)
  
    min_idx <- which.min(cv.trees$dev)
    best_pram <- cv.trees$size[min_idx]
    
    #checking the single node
    if(best_pram==1)
    {
      dev <- cv.trees$dev[-min_idx]
      size <- cv.trees$size[-min_idx]
      
      min_idx2 <- which.min(dev)
      best_pram2 <- size[min_idx2]
      prune.trees <- prune.misclass(tree, best=best_pram2)
    }
    else{
      prune.trees <- prune.misclass(tree, best=best_pram)
    }
 
  
    ##calculate probabilties
    prob[[iter]][,i] <- as.numeric(predict(prune.trees, test_xdt,type='class'))
  
  }
   
    prd_result     <- as.numeric(apply(ifelse(prob[[iter]]>=0.75,0,1), 1, est_cl_ftn))
    dic.acc[iter]  <- mean(test_ydt==prd_result)
  
    cat("***** End of ",iter,"th iteration*****\n")
    iter <- iter+1
}


#Test result
tb.list <- rep(list(0),10)
ans <- c(rep(0,18), rep(1,6))
for(i in 1:10)
{
  p <- as.numeric(apply(ifelse(prob[[i]]>=2,1,0), 1, est_cl_ftn))
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

