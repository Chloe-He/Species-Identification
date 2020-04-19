#######################
##Deep Neural Network##
#######################

#Packages
library(caret)
library(h2o)
search()

#data#
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

h2o.init(nthreads = -1)

performance <- function(tb){
  sensitivity <- tb[2,2]/sum(tb[2,])
  specificity <- tb[1,1]/sum(tb[1,])
  FDR <- tb[1,2]/(tb[2,2]+tb[1,2])
  Accuracy <- (tb[1,1]+tb[2,2])/(tb[1,1]+tb[1,2]+tb[2,1]+tb[2,2])
  F1.score <- tb[2,2]/(2*tb[2,2]+tb[1,2]+tb[2,1]) 
  
  return(list(sensitivity=sensitivity, specificity=specificity, FDR=FDR, Accuracy=Accuracy, F1.score=F1.score))
}

dnn.fit <- list()
tb.list <- list()
iter <- 1


#DNN

while(iter<11){

  set.seed(iter)
  
  ##data partitioning
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
  train <- as.h2o(train) # h2o???? ??????????!

 
  test_xdt <- as.h2o(std_dat[idx,])
  test_ydt <- resp2[idx]



  ##model fitting
  dnn.fit[[iter]] <- h2o.deeplearning(x = 2:ncol(train), #column index
                                      y = 1,
                                      training_frame = train,
                                      distribution = "bernoulli",
                                      activation = "RectifierWithDropout",
                                      hidden = c(ncol(train), ncol(train)*1.6, ncol(train)*1.2, ncol(train)*0.8),
                                      epochs = 300,
                                      loss = "CrossEntropy",
                                      mini_batch_size = 12,
                                      hidden_dropout_ratios = rep(0.5,4), #default:0.5
                                      seed = iter*10
                                      )

  ##prediction
  pred <- h2o.predict(dnn.fit[[iter]], newdata=test_xdt, type="response")
  pred <- as.matrix(pred)
  pp <- c()
  for(i in 1:length(test_ydt)){
    if(as.numeric(pred[i,2]) > as.numeric(pred[i,3])) pp[i] <- "0"
    else pp[i] <- "1"
  }

  tb.list[[iter]] <- table(test_ydt,pp)
  
  cat('End of ',iter,' th iteration')
  
  iter <- iter+1
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------




#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#Test result
tb.list
dnn.res <- unlist(sapply(tb.list, performance))

res <- matrix(dnn.res, 5, 10, 
              dimnames = list(c("sensitivity","specificity", "FDR", "Accuracy", "F1.score"), 
                              c(1:10)))

apply(res,1,mean)

dnn.fit[[1]]


