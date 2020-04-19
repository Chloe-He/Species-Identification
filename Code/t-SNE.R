##Package##
library(Rtsne)
library(randomForest)

#Initialization of matrices#
m <- matrix(0,24,10)
prob <- rep(list(m),10)
tb.list <- list()
iter <- 1
est_cl_ftn  <- function(x) return(names(which.max(table(x))))
performance <- function(tb){
  sensitivity <- tb[2,2]/sum(tb[2,])
  specificity <- tb[1,1]/sum(tb[1,])
  FDR <- tb[1,2]/(tb[2,2]+tb[1,2])
  Accuracy <- (tb[1,1]+tb[2,2])/(tb[1,1]+tb[1,2]+tb[2,1]+tb[2,2])
  F1.score <- tb[2,2]/(2*tb[2,2]+tb[1,2]+tb[2,1]) 
  
  return(list(sensitivity=sensitivity, specificity=specificity, FDR=FDR, Accuracy=Accuracy, F1.score=F1.score))
}

##Data Partition##
while(iter<11){
  
  set.seed(123)
  
  idx1 <- sample(which(resp=="1"),6, replace=FALSE)
  idx2 <- sample(which(resp=="2"),6, replace=FALSE)
  idx3 <- sample(which(resp=="3"),6, replace=FALSE)
  idx4 <- sample(which(resp=="4"),6, replace=FALSE)
  idx <- c(idx1, idx2, idx3, idx4)
  
  resp2 <- resp
  
  resp2[which(resp %in% c("1","2","3"))] <- "0"
  resp2[which(resp=="4")] <- "1"
  for(i in 1:10)
  {
    train_xdt <- alg_dat[-idx,]
    train_ydt <- as.numeric(resp2[-idx])
    
    train <- as.data.frame(cbind(train_ydt, train_xdt))
    train$train_ydt <- as.factor(train$train_ydt)
    
    test_xdt <- alg_dat[idx,]
    test_ydt <- resp2[idx]
    test <- as.data.frame(cbind(test_ydt, test_xdt))
    
    trn <- data.matrix(train)
    tst <- data.matrix(test)
    all <- rbind(trn, tst)
    
    tsne <- Rtsne(as.matrix(all[,2:419]), check_duplicates = FALSE, pca = FALSE, 
                  perplexity=30, theta=0.5, dims=2)
    cols <- rainbow(10)
    plot(tsne$Y, t='n')
    text(tsne$Y, labels=all[,1], col=cols[all[,1] +1])
    all <- cbind(all[,2:419], tsne$Y, all[,1])
    trn <- all[1:72,]
    tst <- all[73:96,]
    
    x1 <- trn[,1:420]
    y1 <- trn[,421]
    x2 <- tst[,1:420]
    y2 <- tst[,421]
    colnames(x1)[419:420] <- c("TSNEx", "TSNEy")
    colnames(x2)[419:420] <- c("TSNEx", "TSNEy")
    rf <- randomForest(x=x1,y=factor(y1),ytest=factor(y2), ntree=500, proximity=TRUE, 
                       importance=TRUE, keep.forest=TRUE, do.trace=TRUE)
    for (k in 1:24)
    {
      prob[[iter]][,i] <- predict(object=rf, newdata=x2)
    }
  }
    prd_result     <- as.numeric(apply(ifelse(prob[[iter]]>=2,1,0), 1, est_cl_ftn))#threshold : 0.75
    dic.acc[iter]  <- mean(test_ydt==prd_result)
    
    cat("********** End of ",iter,"th iteration*********\n")
    iter <- iter+1
}

tb.list <- rep(list(0),10)
ans <- c(rep(0,18), rep(1,6)) ##true species##
  for(i in 1:10)
  {
    p <- as.numeric(apply(ifelse(prob[[i]]>=2,1,0), 1, est_cl_ftn))
    tb.list[[i]] <- table(ans,p)
  }

tb.list

result <- unlist(sapply(tb.list, performance))
mean(dic.acc)

res <- matrix(result, 5, 10, 
              dimnames = list(c("sensitivity","specificity", "FDR", "Accuracy", "F1.score"), 
                              c(1:10)))

apply(res,1,mean)


