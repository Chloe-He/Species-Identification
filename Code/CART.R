library(rpart)
library(caret)
library(rpart.plot)
library(randomForest)
library(mlr)
library(tree)
library(e1071)
alg_dat

idx_mat <- matrix(0,24,10)
iter<-1
m <- matrix(0,24,10)
prob <- rep(list(m),10)
dic.acc <- c()

est_cl_ftn  <- function(x) return(names(which.max(table(x))))
performance <- function(tb){
  sensitivity <- tb[2,2]/sum(tb[2,])
  specificity <- tb[1,1]/sum(tb[1,])
  FDR <- tb[1,2]/(tb[2,2]+tb[1,2])
  Accuracy <- (tb[1,1]+tb[2,2])/(tb[1,1]+tb[1,2]+tb[2,1]+tb[2,2])
  F1.score <- tb[2,2]/(2*tb[2,2]+tb[1,2]+tb[2,1]) 
  
  return(list(sensitivity=sensitivity, specificity=specificity, FDR=FDR, Accuracy=Accuracy, F1.score=F1.score))
}

##CART##
while(iter < 11)
{ 
  set.seed(iter)
  lab_idx1 <- sample(which(resp=="1"), 6)
  lab_idx2 <- sample(which(resp=="2"), 6)
  lab_idx3 <- sample(which(resp=="3"), 6)
  lab_idx4 <- sample(which(resp=="4"), 6)
  lab_idx  <- c(lab_idx1,lab_idx2,lab_idx3,lab_idx4)
  
  idx_mat[,iter]<- lab_idx
  
  resp2 <- resp
  resp2[which(resp %in% c("1","2","3"))] <- "0"
  resp2[which(resp=="4")] <- "1"
  for(i in 1:10)
  {
    train_xdt <- as.data.frame(alg_dat[-lab_idx,])
    train_ydt <- as.vector(resp2[-lab_idx])
    training <- cbind.data.frame(train_xdt,train_ydt)
    
    test_xdt <- as.data.frame(alg_dat[lab_idx,])
    test_ydt <- as.vector(resp2[lab_idx])
    testing <- cbind.data.frame(test_xdt,test_ydt)
    ##find the best tuning params##
    
    #fitControl <- trainControl(method = "repeatedcv",number = 10,repeats = 10)
    #?train
    #model <- train(train_xdt, train_ydt, method = "rpart", trControl = fitControl, tuneLength =10)
    #plot(model)
    #model$bestTune

    rpart.control <- rpart.control(minsplit = 20, cp = 0.03395, maxcompete = 4, 
                                   maxsurrogate = 5, usesurrogate = 2, xval = 10,
                                   surrogatestyle = 1, maxdepth = 30)
                             
    cart <- rpart(train_ydt~.,data=training, control = rpart.control)
    cart

    pred <- predict(cart,testing,type="class")
    pred
    
    #printcp(cart) # display the results 
    #plotcp(cart) # visualize cross-validation results 
    #summary(cart) # detailed summary of splits
    
    #calculate probabilties
    for (k in 1:24)
    {
      prob[[iter]][,i] <- predict(cart, testing,type = "class",na.action = na.pass)
    }
  }
  #voting
  prd_result     <- as.numeric(apply(ifelse(prob[[iter]]>=2,1,0), 1, est_cl_ftn))#threshold : 0.75
  dic.acc[iter]  <- mean(test_ydt==prd_result)
  
  cat("********** End of ",iter,"th iteration*********\n")
  iter <- iter+1
}

tb.list <- rep(list(0),10)
ans <- c(rep(0,18), rep(1,6)) ##true species##
for(i in 1:10){
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

