# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ========  fof   ===============================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#初始化
# source("Default/Index.R")
setwd("~/Code/R/MRSVAR_model")
source('function\\risk_parity_func.R')
source('function\\basic_func.R')
library(Rdonlp2)
library(timeSeries)
library(WindR)
library(xts)
library(TTR)
w.start()


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


path <- "C:/Users/lvjia/Documents/Code/R/MRSVAR_model/data/data0906.csv"
# # 保存数据:只需运行一次
# dataset <- get_data2(start='2002-03-01',end='2016-09-05',type='discrete')
# write.csv(as.data.frame(dataset),"output/data0906.csv")


# 模型1
#资产 恒生+300+商品+黄金
id1 <- c(2,6,13,12)
risk_weight1 <- c(1/4,1/4,1/4,1/4)
# 初始化
num <- length(id1)
u_limit <- rep(1,num)  # 资产权重上限 
l_limit <- rep(0,num)                  # 资产权重下限
u_sigma <- +Inf                  # 年化波动率上限
l_sigma <- -Inf                     #  年化波动率下限
currency <- FALSE                    # 是否进行货币转换
# 设置风险权重配比
# 输出某个时间点往前t天的波动率
asset_vol(id1,c(20,60,90,250),'2015-08-26',path)

risk_parity_result1 <- risk_parity3(id1,risk_weight1,'2016-08-26',u_limit,l_limit,u_sigma,l_sigma,currency,path)
ret1 <- risk_parity_result1$ret
w1 <- risk_parity_result1$weight

# 模型2
#资产 中期国债+中期高信用
id2 <- c(7,8)
risk_weight2 <- c(1/2,1/2)
# 初始化
num <- length(id2)
u_limit <- rep(1,num)  # 资产权重上限 
l_limit <- rep(0,num)                  # 资产权重下限
u_sigma <- +Inf                  # 年化波动率上限
l_sigma <- -Inf                     #  年化波动率下限
currency <- FALSE                    # 是否进行货币转换
# 设置风险权重配比
# 输出某个时间点往前t天的波动率
asset_vol(id2,c(20,60,90,250),'2015-08-26',path)
risk_parity_result2 <- risk_parity3(id2,risk_weight2,'2016-08-26',u_limit,l_limit,u_sigma,l_sigma,currency,path)
ret2 <- risk_parity_result2$ret
w2 <- risk_parity_result2$weight

# 模型3：结合1,2
ret3 <- as.xts(cbind(ret1,ret2))
# 设置风险权重配比
risk_weight3 <- c(3/4,1/4)
# 初始化
num <- 2
u_limit <- rep(1,num)  # 资产权重上限 
l_limit <- rep(0,num)                  # 资产权重下限
u_sigma <- +Inf                  # 年化波动率上限
l_sigma <- -Inf                     #  年化波动率下限
currency <- FALSE                    # 是否进行货币转换
risk_parity_result3 <- risk_parity4(ret3,risk_weight3,u_limit,l_limit,u_sigma,l_sigma,init=c(0.2,0.2))
w3 <- risk_parity_result3$weight

## 计算总权重
w <- c(w1*w3[1],w2*w3[2])
print(round(w,2))















# 回测模型
# 读取本地数据

dataset <- read.csv(path)
dataset1 <- xts(dataset[,id+1],as.Date(dataset[,1]))
dataset1 <- dataset1[paste(start,"/",sep='')]

# x <- xts(dataset[,2:17],as.Date(dataset[,1]))
start <- '2004-05-01'
end <- '2016-09-01'
tdays <- (w.tdays(start,end,'Period=m'))$Data$DATETIME
t <- length(tdays)
weight <- matrix(0,t,num)
risk_weight <- matrix(0,t,num)


# 模型开始
# 初始化上一个调仓日权重
last <- rep(0,num)
# risk parity模型
for (i in 1:t){
  date <- tdays[i]
  print('---------------------------------------------')
  print(paste('Trading date ',as.character(i),'/',as.character(t), ': ',as.character(date),sep=''))
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
print('资产权重')
result2 <- perf_analysis(w2[35:(dim(w2)[1]-1),],na.fill(dataset1[36:dim(w2)[1],],0),cost=0.001,freq='d')

dev.new()
plot.zoo(risk_weight)
print('风险权重')
print(summary(risk_weight))

