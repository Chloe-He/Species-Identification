
####################################################################################
###Decomposition###
# B : Beef
# F : Fish
# P : Pork
# T : Tissue
# Blank
####################################################################################

#-----------------------------------------------------------------------------------
#Set working Directory
setwd("D:/dyu2017/decomposition_csv")
#-----------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
#Loading data 
data <- list()
name <- dir()[1:114]

blank <-substr(name,3,5)=="BLA"
sum(blank)
name2 <- cbind(name, blank)
name <- name2[name2[,2]=='FALSE',1]
nn <- length(name)


symb = substr(name,3,3)
resp = rep(0,nn)
resp[symb=='B'] = 1 # Beef
resp[symb=='F'] = 2 # Fish
resp[symb=='P'] = 3 # Pork
resp[symb=='T'] = 4 # Tissue (Human)


rdata <-list()
for (i in 1:nn){
  rdata[[i]] <- read.table(name[i], sep=";", head=T, stringsAsFactors = T)  
}


save(list=ls(),file='get_data_decomp.rdata')

####################################################################################
rm(list=ls())

load(file='get_data_decomp.rdata')
####################################################################################


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#Checking the data
len_x <-c()
len_y <-c()
for(i in 1:length(rdata)){
  len_x[i] <- ncol(rdata[[i]])-3
  len_y[i] <- nrow(rdata[[i]])
}
len_x   

par(mfrow=c(1,2))
hist(len_x, xlab = "The number of variables(M/Z)", main = "Mass to charge ratio")
legend("center", legend=name[which(len_x<350)], col = 'red', cex=0.85)
hist(len_y, xlab = "The number of observed RT", main = "Retention Time")
legend("center", legend=name[which(len_y!=9643)], col = 'red', cex=0.85)

name[which(len_x<350)] 
name[which(len_y!=9643)]




#Defining Intensities (rowsum of x-variables)
data <-list()
intensities <- c()
for(i in 1:length(rdata)){
  intensities <- apply(rdata[[i]][,c(-1,-2,-3)], 1, sum)  # exclude time variable
  data[[i]] <- cbind(rdata[[i]][,2], intensities)
  colnames(data[[i]]) <- c("r_time", "intensities")
}


#Range of the ret_time 
range_mat = matrix(0,nn,2)
for(i in 1:nn){
  range_mat[i,] = range(data[[i]][,1]) 
}
range_mat
tt_range <- range(range_mat) 
tt_range   #(0.1909667, 41.9930333)



#Range of Intensity
intens_mat <- matrix(0,nn,2)
for(i in 1:nn){
  intens_mat[i,] = range(data[[i]][,2]) 
}
intens_mat
tt_intens <- range(intens_mat)
tt_intens   #(1370, 30637313)

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#Plots
seq(1,nn,by=4)
for(i in seq(1,nn,by=4))
{
  st = 3
  if(substr(name[i],3,3)=="T") {
    ed = 10
  } else {
    st = 3
    ed = 5
  } 
  png(file=paste0('fig/max_',substr(name[i],st,ed),'_peaks_norm_avg.png'),
      width = 1024, height = 1280,type="cairo-png",pointsize=25)
  par(mfrow=c(4,1))
  plot(data[[i]][,1],data[[i]][,2],main=substr(name[i],st,ed),xlab='Retention time',
       ylab='Intensity (1st)',col=1,type='h'
       ,xlim=tt_range,ylim=tt_intens)
  plot(data[[i+1]][,1],data[[i+1]][,2],main=substr(name[i+1],st,ed),xlab='Retention time',
       ylab='Intensity (2nd)',col=1,type='h'
       ,xlim=tt_range,ylim=tt_intens)
  plot(data[[i+2]][,1],data[[i+2]][,2],main=substr(name[i+2],st,ed),xlab='Retention time',
       ylab='Intensity (3rd)',col=1,type='h'
       ,xlim=tt_range,ylim=tt_intens)
  plot(data[[i+3]][,1],data[[i+3]][,2],main=substr(name[i+3],st,ed),xlab='Retention time',
       ylab='Intensity (4th)',col=1,type='h'
       ,xlim=tt_range,ylim=tt_intens)
  dev.off()
}


save(list=ls(),file='before_binning')




#######################################################################################

rm(list=ls())

load(file='before_binning')

#######################################################################################

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
## Binning Procedure ##

bin = 0.01       #initial bin size
le = 0.19        
ri = 40.54

length(seq(le, ri,by=bin))   

bret <- seq(le,ri,by=bin)   
nb_rt <- length(bret)        



#binning process
nn=96 
bre_dat <- matrix(0,nn,nb_rt)
bcount <- matrix(0,nn,nb_rt)
for(i in 1:nn)
{
  temp  = data[[i]] 
  for(k in 1:nb_rt)    
  {
    indx = which(temp[,1]>= bret[k]-bin/2 & temp[,1] < bret[k]+bin/2)  
    #cat(indx)
    if(length(indx)==0) next
    if(length(indx)==1)
      bre_dat[i,k] = temp[indx,2]         
    if(length(indx)>1)
    {
      bcount[i,k] = bcount[i,k] + 1
      bre_dat[i,k] = sum(temp[indx,2])     
    }
    
  }
}

dim(bre_dat)


#plotting with initial binning
bin_forplot <- apply(bre_dat, 2, function(x) sum(x)/length(x))
length(bin_forplot)
plot(x=bret,
     y=bin_forplot, 
     main = "Initial binning",
     xlab = 'Retention Time', 
     ylab = 'Intensitiy',
     xlim = c(0.19,5),
     ylim = c(0,8000000),  
     type='l',
     lwd=2,
     col="darkblue",
     cex.lab=1.5,
     cex.axis=1.5)

bret2 <- seq(le,ri,length=324)
abline(v =bret2, col = 'black', lty=2, lwd=1)

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
## Optimized Bucketing process ##

rat <- 0.05  #optimize binning 위해 bin의 이동 비율 설정
rat/bin
ncol(bre_dat)/(rat/bin) 


source("D://dyu2017//EH_4class_data//opt_binning.R")
fit <- opt_binning(bre_dat,
                   nbins=18*18, 
                   slack=0.95)
optbin_dat <- fit$dat
alg_dat <- log(fit$dat+1) 
bd_t <- bret[fit$index]
nbd_t <- length(bd_t)
tt_area <- range(alg_dat)
tt_area

dim(alg_dat)


#plotting with opt_binning 
bret2 <-seq(le, ri, length.out = 18*18)
optbin_forplot <- apply(fit$dat, 2, function(x) sum(x)/length(x))
length(optbin_forplot)
length(bret2)

plot(x=bret2, 
     y=optbin_forplot,  
     main = "Average spectrum of Opt_binning", 
     xlab = 'Retention Time', 
     ylab = 'Intensitiy',
     xlim = c(0,5),
     ylim=c(0,80000000), 
     type='l',
     lwd=2,
     col="darkblue",
     cex.lab=1.5,
     cex.axis=1.5)
abline(v =bd_t, col = 'black', lty=2, lwd=1)

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
##Spectrum Plots

#plot with no binning
init_rt <- data[[1]]
init_rt <- init_rt[,1]

for(i in 1:nrow(bre_dat)){
  plot(x=init_rt, 
       df[i,], 
       type='l', col=i, 
       xlim=c(0,5), 
       ylim =c(0,55000000),
       main="All spectrums without binning",
       xlab = 'Retention time', 
       ylab = 'Intensitiy',
       cex.lab=1.5,
       cex.axis=1.5)
  
  par(new=TRUE)
}


#init binning plot 
for(i in 1:nrow(bre_dat)){
  plot(x=bret, 
       bre_dat[i,], 
       type='l', col=i, 
       xlim=c(0,5), 
       ylim =c(0,55000000),
       main="All spectrums after Init_binning",
       xlab = 'Retention Time', 
       ylab = 'Intensitiy',
       cex.lab=1.5,
       cex.axis=1.5)
  
  par(new=TRUE)
}



#opt binning plot
for(i in 1:nrow(fit$dat)){
  plot(x=bret2, 
       fit$dat[i,],
       type='l', 
       col=i, 
       xlim=c(0,5), 
       ylim =c(0,550000000),
       main="All spectrums after Opt_binning",
       xlab = 'Retention Time', 
       ylab = 'Intensitiy',
       cex.lab=1.5,
       cex.axis=1.5)
  
  par(new=TRUE)
}



save(list=ls(),file='after_binning')

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
