# author: Jiayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: risk_budget_backtest2.R
# time: 2016/9/6

rm(list = ls())

setwd("~/Code/R/MRSVAR_model")
source('function\\risk_parity_func.R')
source('function\\basic_func.R')
library(WindR)
library(xts)
library(TTR)
library(Rdonlp2)

w.start()


# risk parity


id <- c(2,1,12,17,13)

# 初始化
num <- length(id)
start <- '2004-05-01'
end <- '2016-09-01'
tdays <- (w.tdays(start,end,'Period=m'))$Data$DATETIME
t <- length(tdays)
weight <- matrix(0,t,num)
risk_weight <- matrix(0,t,num)
u_limit <- rep(1,num)  # 资产权重上限 
l_limit <- rep(0,num)                  # 资产权重下限
u_sigma <- +5.1                   # 年化波动率上限
l_sigma <- +4.9                     #  年化波动率下限
currency <- FALSE                    # 是否进行货币转换

# # 保存数据:只需运行一次
# dataset <- get_data2(start='2002-03-01',end='2016-09-05',type='discrete')
# write.csv(as.data.frame(dataset),"C:/Users/lvjia/Documents/Code/R/MRSVAR_model/data/data0906.csv")

# 读取本地数据

dataset <- read.csv("C:/Users/lvjia/Documents/Code/R/MRSVAR_model/data/data0906.csv")
dataset1 <- xts(dataset[,id+1],as.Date(dataset[,1]))
dataset1 <- dataset1[paste(start,"/",sep='')]




risk_parity3(id,risk_weight,date='2016-08-26',u_limit=rep(1,9),l_limit=rep(0,9),u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,window=60,init=c(0))



x <- xts(dataset[,2:17],as.Date(dataset[,1]))
# 获取日期
# 设置风险权重
risk_weight0 <- c(1/4,1/4,1/4,1/4)
# 初始化上一个调仓日权重
last <- rep(0,num)
# risk parity模型
for (i in 1:t){
  date <- tdays[i]
  print('---------------------------------------------')
  print(paste('START : ',as.character(date),sep=''))
  k <- 1
  if (i>1){
    # result <- risk_parity2(id,date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=last)
    result <- risk_parity3(id,risk_weight0,date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=last)
    # 记录最优结果
    best_weight <- result$weight
    best_risk_weight <- result$risk_weight
    best_error <- result$error
    while (sum(abs(result$weight-last))>0.1*num & k < 10){
    if (k>3){
      init <- runif(num,min=0,max=1)
    }else{
      init <- last
    } 
    # result <- risk_parity2(id,date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=init)
    result <- risk_parity3(id,risk_weight0,date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=last)
    # 打压过程
    print('############## Re optimizing ################')
    print(paste('Trial: ',as.character(k),sep=''))
    print(round(result$weight,2))
    print(round(last,2))
    print(sum(abs(result$weight-last)))
    # 是否好于历史最优值
    if (best_error > result$error){
      print('UPDATE')
      best_error <- result$error
      best_risk_weight <- result$risk_weight
      best_weight <- result$weight
    }
    k <- k+1
    }
    
  }else{
    # result <- risk_parity2(id,date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=last)
    result <- risk_parity3(id,risk_weight0,date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=last)
    # 初始化最优值
    best_error <- result$error
    best_risk_weight <- result$risk_weight
    best_weight <- result$weight    
    
  }
  # 记录调仓期的结果
  weight[i,] <- best_weight
  risk_weight[i,] <- best_risk_weight
  last <- result$weight 
  # print(result$risk_weight)
}

# 优化调仓
weight_adj <- adjust_weight(weight,qu=0.05,type=1)
w1 <- xts(weight_adj$weight_adj,tdays)
# 数据整合
xx <- merge.xts(dataset1,w1)
w2 <- xx[,(num+1):(2*num)]
w2 <- na.locf(w2)
# 模型表现输出
result2 <- perf_analysis(w2[35:(dim(w2)[1]-1),],na.fill(dataset1[36:dim(w2)[1],],0),cost=0.001,freq='d')

dev.new()

plot.zoo(risk_weight)
print('风险权重')
print(summary(risk_weight))

