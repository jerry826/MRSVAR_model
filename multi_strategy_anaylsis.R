source('MRS.R')

# strat1: 等权重组：每周再平衡
strategy_1 <- function(dataset,freq){
  m <- dim(dataset)[2]
  n <- dim(dataset)[1]
  weight <- matrix(0.2,n,m)
  perf <- list(weight=weight)
  weight2 <- adjust_weight(perf)
  perf_analysis(weight2$weight_adj,dataset,freq=freq)
  print(dim(weight2$weight_adj))
  
}

# strat2: 等权重组：无再平衡
strategy_2 <- function(dataset,freq){
  m <- dim(dataset)[2]
  n <- dim(dataset)[1]
  weight <- matrix(0,n,m)
  cum_ret <- cumprod(dataset+1)
  weight[1,] <- rep(0.2,m)
  weight[2:n,] <- cum_ret[1:n-1,]/apply(cum_ret[1:n-1,],1,sum)
  
  # weight <- matrix(rep(0.2,m*n),n)
  perf <- list(weight=weight)
  weight2 <- adjust_weight(perf,qu=0.05,type=2)
  perf_analysis(weight2$weight_adj,dataset,freq=freq)
  print(dim(weight2$weight_adj))
}
# strat3: 等权

# 趋势追踪策略
## 每一期购买上一期表现最好的组合
strategy_3 <- function(dataset,freq){
  m <- dim(dataset)[2]
  n <- dim(dataset)[1]
  id <- apply(dataset,1,which.max)
  weight <- matrix(0.0,n,m)
  weight[1,] <- rep(0.0,m)
  for (i in 1:(n-1)){
    weight[i+1,id[i]] <- weight[i+1,id[i]]+1
  }
  
  perf <- list(weight=weight)
  weight2 <- adjust_weight(perf,qu=0.05,type=2)
  perf_analysis(weight2$weight_adj,dataset,freq=freq)  
  print(dim(weight2$weight_adj))
  
}

# 风险平价模型
strategy_4 <- function(dataset,freq){
  window <- 60
  h <- weight_cal(covar2(dataset,window))
  perf <- list(weight=h)
  # 仓位再优化
  weight2 <- adjust_weight(perf)
  
  asset_ret <- dataset
  n <- dim(dataset)[1]
  asset_ret <- dataset[(window+1):n]
  perf_analysis(weight2$weight_adj[1:(n-window),],asset_ret ,cost=0.001)
  print(dim(asset_ret))
  cum_p <- cbind(cumprod(apply(weight2$weight_adj[1:(n-window),]*asset_ret,1,sum)+1),index(asset_ret))
  write.csv(cum_p,'result.csv')
}



assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE','GF0002.GF')

w.start()
start <- '2005-02-01'
end <- '2016-08-16'
mid <- '2013-08-01'
freq <- 'w'
assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE')

# assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE','GF0002')
ret_all <- collect_data(assets,start,end,mid,freq)
# write.csv(cumprod(1+ret_all$test_d),file='data.csv')
plot.zoo(cumprod(ret_all$test+1))
dataset <- ret_all$test_d

# 单个分析
for (i in 1:length(assets)){
  ret_analysis(dataset[,i],freq=freq)
  
} 
source('MRS.R')


strategy_1(dataset,freq)
strategy_2(dataset,freq)
strategy_3(dataset,freq)  


start <- '2005-02-01'
end <- '2016-08-16'
mid <- '2012-05-20'
freq <- 'w'
assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE')

# assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE','GF0002.GF')
ret_all <- collect_data(assets,start,end,mid,freq)
dataset <- ret_all$test_d

strategy_4(dataset,freq)
