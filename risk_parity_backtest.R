# author: Jiayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: risk_parity_backtest.R
# time: 2016/9/2

rm(list = ls())
source('risk_parity_func.R')
library(timeSeries)
library(zoo)
library(parallel)

adjust_weight <- function(weight,qu=0.05,type=1){
  # 调仓方式优化
  # weight:初始仓位
  # qu    :阈值
  # type  :调整方法:1根据最大调仓量;2根据总调仓量
  n <- dim(weight)[1]
  m <- dim(weight)[2]
  
  weight_adj <- matrix(rep(0,m*n),n)
  weight_adj[1,] <- as.numeric(weight[1,])
  
  
  weight_change <- matrix(rep(0,n),n)
  weight_change[1] <- sum(abs(weight_adj[1,]))
  
  change_pos <- matrix(rep(0,n),n)
  change_pos[1] <- TRUE
  
  
  for (i in 2:n){
    if (type==1){
      delta_w <- (max(abs(weight[i,]-weight_adj[i-1,])))
    }else{
      delta_w <- (sum(abs(weight[i,]-weight_adj[i-1,])))
    }
    
    if (delta_w > qu){
      # 超过阈值,调仓
      change_pos[i] <- TRUE
      weight_adj[i,] <- as.numeric(weight[i,])
    }else{
      # 未超过阈值,不调仓
      change_pos[i] <- FALSE
      weight_adj[i,] <- as.numeric(weight_adj[i-1,])
    }
    # 调整后的每日资产持仓变化量
    delta_w1 <-  (sum(abs(weight_adj[i,]-weight_adj[i-1,])))
    weight_change[i] <- delta_w1
  }
  return(list(weight=weight,weight_adj=weight_adj,weight_change=weight_change,change_pos=change_pos)) 
}

perf_analysis <- function(weight,asset_ret,cost=0.001,freq='m'){
  # 策略表现分析
  # weight   :权重list
  # asset_ret:对应资产收益率
  # cost     :单边交易手续费费
  # 权重序列长度
  n <- dim(weight)[1]
  
  if (freq=='w'){
    m <- 52
  } else if(freq=='m'){
    m <- 12
  } else{m <- 252}
  
  ret_before_fee <- apply(weight*as.matrix(asset_ret),1,sum)
  ret_after_fee <- matrix(rep(0,n-1))
  ret_after_fee[1] <- ret_before_fee[1]-cost*sum(abs(weight[1,]))
  for (i in (2:n)){
    ret_after_fee[i] <- ret_before_fee[i]-cost*sum(abs(weight[i,]-weight[i-1,]))
    
  }
  
  cum_ret_bf <- cumprod(ret_before_fee+1)
  cum_ret_af <- cumprod(ret_after_fee+1)
  
  # annual return
  mean_ret_bf <- (cum_ret_bf[n])^(m/n)-1
  mean_ret_af <- (cum_ret_af[n])^(m/n)-1
  
  # mean_ret_af <- (mean(ret_after_fee)+1)^m 
  sharpe_bf <- sharpe(ret_before_fee,freq)
  sharpe_af <- sharpe(ret_after_fee,freq)
  
  # annual volatility
  vol_bf <- sqrt(m)*sd(ret_after_fee)
  vol_af <- sqrt(m)*sd(ret_after_fee)
  
  # draw_down
  drawdown_bf <-  min(drawdowns(timeSeries(ret_before_fee)))
  drawdown_af <-  min(drawdowns(timeSeries(ret_after_fee)))
  
  # before fee
  print('----------Before fee result----------')
  print(paste('Annual return: ',as.character(round(mean_ret_bf,4)) ,sep=''))
  print(paste('Annual volatility: ',as.character(round(vol_bf,4)) ,sep=''))
  print(paste('Sharpe ratio: ',as.character(round(sharpe_bf,3)) ,sep=''))
  print(paste('Maximum drawdown: ',as.character(round(drawdown_bf,4)) ,sep=''))
  
  
  # after fee
  print('----------After fee result----------')
  print(paste('Annual return: ',as.character(round(mean_ret_af,4)) ,sep=''))
  print(paste('Annual volatility: ',as.character(round(vol_af,4)) ,sep=''))
  print(paste('Sharpe ratio: ',as.character(round(sharpe_af,3)) ,sep=''))
  print(paste('Maximum drawdown: ',as.character(round(drawdown_af,4)) ,sep='')) 
  
  print(summary(weight))
  dev.new()
  plot.zoo(xts(cbind(cum_ret_bf,cum_ret_af),index(asset_ret)),col=c('red','blue'),plot.type="single")
  dev.new()
  plot.zoo(xts(weight,index(asset_ret)))
  return(list(cum_ret_af=cum_ret_af,cum_ret_bf=cum_ret_bf,dd=drawdowns(timeSeries(ret_after_fee)),weight=weight))
}

sharpe <- function(ret,freq){
  # 计算年化夏普值
  # ret:收益率序列
  # freq:收益率序列频率
  print(freq)
  if (freq=='w'){
    sharpe <- sqrt(52)*mean(ret)/sqrt(var(ret))
  }else if (freq=='m'){
    sharpe <- sqrt(12)*mean(ret)/sqrt(var(ret))
  }else if (freq=='d'){
    sharpe <- sqrt(252)*mean(ret)/sqrt(var(ret))
  }else{
    print('ERROR')
  }
  return(sharpe)
}

func1 <- function(t){
  source('risk_parity_func.R')
  source('func.R')
  library(timeSeries)
  library(zoo)
  data1 <- (get_data(start='2010-08-10',end=t,type='continuous'))[,c(3,6,7)]
  data1 <- apply.weekly(data1,colSums)
  mrs_m <- msvar(data1,p=1,h=5,niterblkopt=30)
  last_p <- mrs_m$fp[dim(mrs_m$fp)[1],]
  next_p <- t(last_p) %*% mrs_m$Q
  last_ret <- as.matrix(data1[dim(data1)[1],])
  next_ret_est <- matrix(0,3,1)
  for (i in 1:5){
    beta <- as.matrix(mrs_m$hreg$Bk[,,i])
    next_ret_est <- next_ret_est+t(beta)%*%(as.matrix(c(as.numeric(last_ret),1)))*next_p[i]
  }
  # w0 <- rep(0,3)
  # w0[which.max(next_ret_est)] <- 1
  return(t(next_ret_est))
}


# 1. H11137.CSI      中国互联网(美元)      Equity    2007-07-01
# 2. HSCEI.HI        恒生国企指数(港币)    Equity    2005-04-14
# 3. H11017.CSI      中期国债              Bond      2008-01-01
# 4. CL.NYM NYMEX    原油(美元)            Commodity 2011-12-21
# 5. AU9999.SGE      黄金                  Commodity 2004-01-01
# 6. 000016.SH       上证50                Equity    2004-01-01
# 7. 399006.SZ       创业板指              Equity    2010-06-02

# 8. 070025.OF       嘉实信用A             Bond      2011-09-14
# 9. 000198.OF       天弘余额宝            Bond      2013-06-03
#    482002.OF       工银货币                        2006-03-26
# 设置参数
start <- '2013-10-01'
end <- '2016-09-01'

u_limit <- c(0,1,1,1,1,1,1,1,0)  # 资产权重上限 
l_limit <- rep(0,9)                  # 资产权重下限
u_sigma <- +Inf                      # 年化波动率上限
l_sigma <- -Inf                      # 年化波动率下限
currency <- FALSE                    # 是否进行货币转换

# 获取日期并初始化

tdays <- (w.tdays(start,end,'Period=m'))$Data$DATETIME
t <- length(tdays)
weight <- matrix(0,t,9)


for (i in 1:t){
  date <- tdays[i]
  print(date)
  # l_limit1 <-index1[i,]*0.2
    
  # w <- risk_parity(date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency)
  if (i > 20){
    opt_weight <- risk_parity(date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=last_weight)
    k <-1
    # 优化方法
    ## 如果优化结果与前一天差距太大，则重新优化直到满足
    ## 1:10次  用前一次结果
    ## 11:50次 用随机数
    
    while (sum(abs(opt_weight-last_weight))>0.50 & k < 2){
      if (k>3){
        init <- init<-runif(9,min=0,max=1)
      }else{
        init <- last_weight
      }
      opt_weight <- risk_parity(date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,init=init)
      print('############## Re optimizing ################')
      print(paste('Trial: ',as.character(k),sep=''))
      print(round(opt_weight,2))
      print(round(last_weight,2))
      print(sum(abs(opt_weight-last_weight)))
      k <- k+1
    }
    weight[i,] <- opt_weight
    last_weight <- opt_weight
    
  }else{
    opt_weight <- risk_parity(date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency)
    weight[i,] <- opt_weight
    last_weight <- opt_weight
}
}

# 权重优化
weight_adj <- adjust_weight(weight,qu=0.05,type=1)

ww <- xts(weight_adj$weight_adj,tdays)

# 获取收益率
ret <- get_data(start=start,end=end,type='discrete')
xx <- merge.xts(ret,ww)
ww <- xx[,10:18]
ww <- na.locf(ww)

# 计算表现
result <- perf_analysis(ww[21:713,],na.fill(ret[22:714,],0),cost=0.001,freq='d')

