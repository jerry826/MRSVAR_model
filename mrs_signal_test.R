# author:  Jiayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: mrs_signal.R
# time: 2016/9/12
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ========  markov regime switch model signal generator  ===============================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#初始化
# source("Default/Index.R")
setwd("~/Code/R/")
source('MRSVAR_model/mrs_func.R',encoding='utf-8')
source('QuantFunction/basic_func.R',encoding='utf-8')
source('RiskParityModel/risk_parity_func.R',encoding='utf-8')

library(MSBVAR)
library(WindR)
library(xts)
library(TTR)
library(expm)
library(Rdonlp2)
library(mvtnorm)
library(zoo)
library(vars)
library(fBasics)
library(TTR)
library(parallel)




w.start()

# 设置数据路径
path <- "MRSVAR_model/input/data_c.csv"
# # 保存数据:只需运行一次
# dataset <- get_data2(start='2002-03-01',end='2016-09-05',type='discrete')
# write.csv(as.data.frame(dataset),file="MRSVAR_model/input/data_d.csv")

# 1.  H11137.CSI      中国互联网(美元)             Equity    2007-07-01
# 2.  HSCEI.HI        恒生国企指数(港币)           Equity    2005-04-14
# 3.  000016.SH       上证50                       Equity    2004-01-01
# 4.  399006.SZ       创业板指                     Equity    2010-06-02
# 5.  000905.SH       中证500                      Equity    2005-01-04
# 6.  000300.SH       沪深300                      Equity    2002

# 7.  H11017.CSI      中期国债                     Bond      2008-01-01
# 8.  H11078.CSI      中证中高信用                 Bond      2008-01-01
# 9.  0401E.CS        中债固定利率全价(7-10)指数   Bond      2006-11-20
# 14  070025.OF       嘉实信用A                    Bond      2011-09-14
# 10. 037.CS          中债总财富指数               Bond      2002-01-01

# 15. 000198.OF       天弘余额宝                   Cash      2013-06-03
# 16. 482002.OF       工银货币                     Cash      2006-03-26

# 11. CL.NYM NYMEX    原油(美元)                   Commodity 2011-12-21
# 12. AU9999.SGE      黄金                         Commodity 2004-01-01
# 13. CCFI.WI         WIND商品                     Commodity 2002-01-01

path1 <- "MRSVAR_model/input/data_c.csv"
path2 <- "MRSVAR_model/input/data_d.csv"

dt <- '2016-08-26'
id <- c(5,6,10)
h <- 5
num <- length(id)



end <- '2016-08-30'   # 最后一期
start <- '2013-08-01'   # 第一期
tdays <- (w.tdays(start,end,'Period=m'))$Data$DATETIME
t <- length(tdays)

weight <- matrix(0,t,9)

# weight <- matrix(0,t,9)
input <- as.data.frame(tdays)


####################################################
#####  Part 1.0: mrsvar mean=variance signal   #####
####################################################

system.time({
  func <- function(dt){
    setwd("~/Code/R/")
    source('MRSVAR_model/mrs_func.R',encoding='utf-8')
    library(MSBVAR)
    library(WindR)
    library(xts)
    library(Rdonlp2)
    w.start()
    rr <- mrs_signal(dt=dt)
    return(rr)
  }
  x <- 1:t
  cl <- makeCluster(4) # 初始化四核心集群
  results <- parLapply(cl,tdays,func) # lapply的并行版本
  res.df <- do.call('rbind',results) # 整合结果
  stopCluster(cl) # 关闭集群
})

cc2 <- matrix(0,37,14)
for (i in 1:t){
  cc <- func(tdays[i])
  cc2[i,] <- cc
}




# 处理信号
re2 <- xts(cc2,tdays)
ww1 <- re2[,1:3]  # 分开
ww2 <- re2[,4:6]  # 加总优化
rr1 <- re2[,7:9] # 预测收益率
pp <- re2[,10:14]
# 37 signals

weight <- xts(adjust_weight(ww2)$weight_adj,tdays)
# 产生回测
dataset <- read.csv(path)
ret <- xts(dataset[,2:dim(dataset)[2]],as.Date(dataset[,1]))
ret <- ret[paste(tdays[1],'/',sep=''),id]
backtest <- na.locf(merge(ret,weight))
perf <-perf_analysis(backtest[,4:6],backtest[,1:3],freq='d')




######################################
#####  Part 1.1: mrsvar signal   #####
######################################
# 产生权重
signal1 <- apply(rr1,1,which.max)
# weight1 <- matrix(0.2,t,num)
# for (i in 1:t){
#   weight1[i,signal1[i]] <- 0.4
# }
# weight <- xts(adjust_weight(weight1)$weight_adj,tdays)
# # 产生回测
# dataset <- read.csv(path)
# ret <- xts(dataset[,2:dim(dataset)[2]],as.Date(dataset[,1]))
# ret <- ret[paste(tdays[1],'/',sep=''),id]
# backtest <- na.locf(merge(ret,weight))
# perf <-perf_analysis(backtest[,5:8],backtest[,1:4],freq='d')


########################################
#####    Part 1.2: trend signal    #####
########################################


dataset1 <- read.csv(path1)
ret <- xts(dataset1[,2:dim(dataset1)[2]],as.Date(dataset1[,1]))
ret <- ret[paste(start,'/',sep=''),id]

ret_actual_5lag <- apply.monthly(lag(ret,5)[6:dim(ret)[1]],colSums)  # 5 days lag trend
signal2<- apply(ret_actual_5lag[1:t,],1,which.max)
# weight2 <- matrix(0,t,num)
# for (i in 1:(t)){
#   weight2[i,signal2[i]] <- 1
# }
# weight_adj <- xts(adjust_weight(weight2)$weight_adj,tdays[1:37])

# 产生回测
# ret_raw <- read.csv(path2)
# backtest_ret <- xts(ret_raw[,2:dim(ret_raw)[2]],as.Date(ret_raw[,1]))
# backtest_ret <- backtest_ret[paste(tdays[1],'/',sep=''),id]
# backtest <- na.locf(merge(backtest_ret,weight_adj))
# perf <-perf_analysis(backtest[,5:8],backtest[,1:4],freq='d')




############################################
#####       Part 1.3: best signal      #####
############################################

dataset1 <- read.csv(path1)
ret <- xts(dataset1[,2:dim(dataset1)[2]],as.Date(dataset1[,1]))
ret <- ret[paste(start,'/',sep=''),id]
ret_actual <- apply.monthly(ret,colSums)
signal3 <- apply((lag(ret_actual,-1))[1:37,],1,which.max)


##################################################
#####    Part 1.4: signal efficient test     #####
##################################################

# choose a signal: sigbal1, sigbnal2, signal3
signal <- signal1
# adjust weight
w0 <- 0
weight <- matrix(w0,t,num)
for (i in 1:(t)){
  weight[i,signal[i]] <- 1-3*w0
}
weight_adj <- xts(adjust_weight(weight)$weight_adj,tdays)

# backtesting
ret_raw <- read.csv(path2)
backtest_ret <- xts(ret_raw[,2:dim(ret_raw)[2]],as.Date(ret_raw[,1]))
backtest_ret <- backtest_ret[paste(tdays[1],'/',sep=''),id]
backtest <- na.locf(merge(backtest_ret,weight_adj))
perf <-perf_analysis(backtest[,(num+1):(2*num)],backtest[,1:num],freq='d')


####################################################
#####         Part 2.1: risk parity            #####
####################################################


weight<- matrix(0,t,num)
risk_weight <- matrix(0,t,num)
u_limit <- rep(1,num)  # 资产权重上限 
l_limit <- rep(0,num)                  # 资产权重下限
u_sigma <- +Inf                      # 年化波动率上限
l_sigma <- -Inf                      # 年化波动率下限
currency <- FALSE                    # 是否进行货币转换
signal <- as.data.frame(signal2)

for (i in 1:t){
  risk_weight0 <- rep(0.3,num)
  # risk_weight0 <- c(0.2)
  
  risk_weight0[signal[i,]] <- 0.4
  # risk_weight[floor(runif(1)*4)+1] <- 0.7
  result_temp <- risk_parity3(id,risk_weight0,date=tdays[i],u_limit=u_limit,l_limit=l_limit,u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,path=path1,window=60,init=c(0))
  k <- 1
  while (result_temp$error > 0.001 & k <= 20){
    print('===========      re     ============')
    result_temp <- risk_parity3(id,risk_weight0,date=tdays[i],u_limit=u_limit,l_limit=l_limit,u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,path=path1,window=60,init=c(0))
    k <- k+1
  }
  risk_weight[i,] <- result_temp$risk_weight
  weight[i,] <- result_temp$weight
} 
weight_adj <- xts(adjust_weight(weight)$weight_adj,tdays)
# backtesting
ret_raw <- read.csv(path2)
backtest_ret <- xts(ret_raw[,2:dim(ret_raw)[2]],as.Date(ret_raw[,1]))
backtest_ret <- backtest_ret[paste(tdays[1],'/',sep=''),id]
backtest <- na.locf(merge(backtest_ret,weight_adj))
perf <-perf_analysis(backtest[,(1+num):(2*num)],backtest[,1:num],freq='d')



####################################################
#####         Part 2.2: risk parity            #####
####################################################
q1 <- quantile(re2[,1])
q2 <- quantile(re2[,2])
q4 <- quantile(re2[,4])



weight<- matrix(0,t,num)
risk_weight <- matrix(0,t,num)
u_limit <- rep(1,num)  # 资产权重上限 
l_limit <- rep(0,num)                  # 资产权重下限
u_sigma <- +Inf                      # 年化波动率上限
l_sigma <- -Inf                      # 年化波动率下限
currency <- FALSE                    # 是否进行货币转换
signal <- as.data.frame(signal3)

for (i in 1:t){
  risk_weight0 <- rep(0.,4)
  if (re2[i,1]>0.0){
    risk_weight0[1] <- 15
  }else if(re2[i,1]<q1[2]){
    risk_weight0[1] <- 1
  }else{
    risk_weight0[1] <- 3
  }

  if (re2[i,2]>0.02){
    risk_weight0[2] <- 3
  }else if(re2[i,2]<q2[2]){
    risk_weight0[2] <- 1
  }else{
    risk_weight0[2] <- 3
  }
  
  if (re2[i,4]>q4[4]){
    risk_weight0[4] <- 3
  }else if(re2[i,4]<q4[2]){
    risk_weight0[4] <- 1
  }else{
    risk_weight0[4] <- 2
  }  
  
  risk_weight0[3] <- 2
  risk_weight0 <- risk_weight0/sum(risk_weight0)
  
  # risk_weight[floor(runif(1)*4)+1] <- 0.7
  result_temp <- risk_parity3(id,risk_weight0,date=tdays[i],u_limit=u_limit,l_limit=l_limit,u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,path=path1,window=60,init=c(0))
  k <- 1
  while (result_temp$error > 0.001 & k <= 10){
    print('===========      re     ============')
    result_temp <- risk_parity3(id,risk_weight0,date=tdays[i],u_limit=u_limit,l_limit=l_limit,u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,path=path1,window=60,init=c(0))
    k <- k+1
  }
  risk_weight[i,] <- result_temp$risk_weight
  weight[i,] <- result_temp$weight
} 
weight_adj <- xts(adjust_weight(weight)$weight_adj,tdays)
# backtesting
ret_raw <- read.csv(path2)
backtest_ret <- xts(ret_raw[,2:dim(ret_raw)[2]],as.Date(ret_raw[,1]))
backtest_ret <- backtest_ret[paste(tdays[1],'/',sep=''),id]
backtest <- na.locf(merge(backtest_ret,weight_adj))
perf <-perf_analysis(backtest[,5:8],backtest[,1:4],freq='d')

plot.zoo(risk_weight)
