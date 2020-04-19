##packages##
library(e1071)
library(rpart)

#Initialization#
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
#make training and testing#
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

  train_xdt <- std_dat[-idx,]
  train_ydt <- as.numeric(resp2[-idx])
  train <- as.data.frame(cbind(train_ydt, train_xdt))
  train$train_ydt <- as.factor(train$train_ydt)

  test_xdt <- std_dat[idx,]
  test_ydt <- resp2[idx]
  test <- as.data.frame(cbind(test_ydt, test_xdt))
  test$test_ydt <- as.factor(test$test_ydt)
  
  for(i in 1:10)
  {
  #find the best tuning#
  tune.out <- tune(svm, train_xdt,train_ydt, data = train, kernel="polynomial",
                 ranges = list(cost = c(0.1,1,10,100,1000),
                               gamma = c(0.5,1,2,3,4)))
  tune.out$best.model

  #build the model#
  library(e1071)
  svmfit <- svm(train_ydt~.,data=train,gamma=3,cost=0.1,epsilon=0.1,kernel="polynomial")
  summary(svmfit)

  #prediction on testing#
  prob[[iter]][,i] <- as.numeric(predict(svmfit, test_xdt, type='class'))
  }
  
  cat("***** End of ",iter,"th iteration*****\n")
  iter <- iter+1
  
}

prob
tb.list <- rep(list(0),10)
ans <- c(rep(0,18), rep(1,6))
for(i in 1:10){
  p <- as.numeric(apply(ifelse(prob[[i]]>=2,1,0), 1, est_cl_ftn))
  tb.list[[i]] <- table(ans,p)
}

tb.list

#get results#
result <- unlist(sapply(tb.list, performance))

res <- matrix(result, 5, 10, 
              dimnames = list(c("sensitivity","specificity", "FDR", "Accuracy", "F1.score"), 
                              c(1:10)))

apply(res,1,mean)

# visualize (classes by color, SV by crosses):
plot(cmdscale(dist(train[,-419])),
     col = as.integer(train[,419]),
     pch = c("o","+")[1:150 %in% svmfit$index + 1])

