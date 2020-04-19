setwd('C:/Users/dyu/Dropbox/18_NFS_odor/raw_data_csv')
load("Decomposition_Odor_Raw_data (1).Rdata")

ndata = length(raw_dat) 
#dat1d_list #this is a dataset containing RT in minutes and Intensity as summation#
dat1d_list = list()
for(i in 1:ndata)
{
  dat1d_list[[i]] = data.frame(rt=raw_dat[[i]][,2],ints = apply(raw_dat[[i]][,-(1:3)],1,sum))
}
head(dat1d_list[[1]])

range_mat = matrix(0,ndata,2)
for(i in 1:ndata)
{
  range_mat[i,] = range(dat1d_list[[i]][,1])
}
rt_range = range(range_mat)
rt_range


# Binning Procedure

bin = 0.02
le = 0.2 #smallest RT
ri = 41.9 #largest RT
bret = seq(le,ri,by=bin)
nb_rt = length(bret)
nb_rt ## total number of bins

nn=96
bre_dat = matrix(0,nn,nb_rt)  #create a matrix with nn rows and nb_rt columns#
bcount = matrix(0,nn,nb_rt) #create a matrix with nn rows and nb_rt columns#
for(i in 1:nn)  # for every matrix#
{
  temp  = dat1d_list[[i]] #get any one of the 96 matrices#
  for(k in 1:nb_rt)  #for the kth bin#
  {
    #in that particular set, extract the index in which RT are in the range#
    indx = which(temp[,1]>= bret[k]-bin/2 & temp[,1] < bret[k]+bin/2) 
    #cat(indx)
    if(length(indx)==0) next
    if(length(indx)==1)
    {
      bcount[i,k] = 1 ## count=1 if there's ONE observation that has RT falls into the bin
      bre_dat[i,k] = temp[indx,2] ## Fill in that ONE observation's intensity
    }
    if(length(indx)>1)
    {
      ## count=length if there're more than ONE observation that have RT falls into the bin
      bcount[i,k] = length(indx) 
      bre_dat[i,k] = sum(temp[indx,2]) ## Fill in the SUM of those intensity  
    }
    
  }
}

dim(bre_dat)
tt_area = range(bre_dat)
tt_area

rat = 0.1
rat/bin
ncol(bre_dat)/(rat/bin) #round up=418

source('opt_binning.R')
fit = opt_binning(bre_dat,nbins=418,slack=0.9)
alg_dat = log(fit$dat+1)
bd_t = bret[fit$index[-1]]
nbd_t = length(bd_t)
ints_range = range(alg_dat)
dim(alg_dat)
nbd_t

resp = spec_class #1,2,3,4 corresponding to species#

set.seed(1111) ##set 75% as training, the rest as testing set##
tr_indx = list()
ts_indx = list()
for(i in 1:4)
{
  temp = which(spec_class==i)
  tr_indx[[i]] = sample(temp,18)
  ts_indx[[i]] = setdiff(temp,tr_indx[[i]])
}
tr_indx = unlist(tr_indx)
ts_indx = unlist(ts_indx)


chk = apply(alg_dat!=0,2,sum)
slt_dat = alg_dat[,chk>1]
sum(chk>1)

tr_dat = slt_dat[tr_indx,]
tr_resp = resp[tr_indx]

ts_dat = slt_dat[ts_indx,]
ts_resp = resp[ts_indx]

n_tr = length(tr_indx)
n_ts = length(ts_indx)

library(glmnet)

cand_lam = seq(0.01,5,by=0.05)
cand_al = seq(0,1,by=0.1)
##trying different lambda##
test = glmnet(tr_dat,tr_resp,family='multinomial',alpha=0,lambda=cand_lam,standardize=F)
prd = predict(test,newx=tr_dat,type='response',s=2.96,standardize=F)
est = coef(test,s=2.96)

est_class = apply(prd,1,function(x) which(x==max(x)))
cbind(tr_resp,est_class)
tb = table(tr_resp,est_class)
tb ## perfect classification ##

cand_lam = seq(0.01,5,by=0.05)
cand_al = seq(0,1,by=0.05)


cv_mat = matrix(0,length(cand_lam),length(cand_al))
colnames(cv_mat) = cand_al
rownames(cv_mat) = cand_lam[length(cand_lam):1]

for(i in 1:length(cand_al))
{
  alpha = cand_al[i]
  cvfit=cv.glmnet(tr_dat,tr_resp,family='multinomial',alpha=alpha,lambda=cand_lam,
                  type.measure="deviance",nfold=5,standardize=F)
  cv_mat[,i] = cvfit$cvm
  cat(i,"-th iteration \n")
}
apply(cv_mat,2,min)
min(apply(cv_mat,2,min))
which.min(cv_mat[,1])
cv_mat[60:100,] #choose tuning parameter alpha and lambda#

lam_fin = 0.21 
al_fin = 0
test = glmnet(tr_dat,tr_resp,family='multinomial',alpha=al_fin,lambda=cand_lam,standardize=F)
est_coef = predict(test,type='coefficients',s=lam_fin,standardize=F)

prd_res = predict(test,newx=tr_dat,type='response',s=lam_fin,alpha=al_fin,standardize=F)
est_class = apply(prd_res,1,function(x) which(x==max(x)))
cbind(tr_resp,est_class)
tb = table(tr_resp,est_class)
sum(diag(tb))/n_tr*100
#[1] 100
# overfitting issue

prd_ts = predict(test,newx=ts_dat,type='response',s=lam_fin,alpha=al_fin,standardize=F)
est_ts = apply(prd_ts,1,function(x) which(x==max(x)))
cbind(ts_resp,est_ts)
tb_ts = table(ts_resp,est_ts)
sum(diag(tb_ts))/n_ts*100
#[1] 75
library(caret)
auc_res <- multiclass.roc(ts_resp,est_ts) 
auc_res


