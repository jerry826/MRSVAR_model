# author: Jiayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: basic_func.R
# time: 2016/9/6


sharpe <- function(ret,freq){
  # 计算年化夏普值
  # ret:收益率序列
  # freq:收益率序列频率
  if (freq=='w'){
    sharp <- sqrt(52)*mean(ret)/sqrt(var(ret))
  }else if (freq=='m'){
    sharp <- sqrt(12)*mean(ret)/sqrt(var(ret))
  }else if (freq=='d'){
    sharp <- sqrt(252)*mean(ret)/sqrt(var(ret))
  }
  return(sharp)
}

adjust_weight <- function(weight,qu=0.05,type=1){
  # 调仓方式优化
  # perf:mrs mv模型优化结果
  # qu  :阈值
  # type:调整方法:1根据最大调仓量;2根据总调仓量
  n <- dim(weight)[1]
  m <- dim(weight)[2]
  
  weight <- as.data.frame(weight)
  
  weight_adj <- matrix(rep(0,m*n),n)
  weight_adj[1,] <- as.numeric(weight[1,])
  
  
  weight_change <- matrix(rep(0,n),n)
  weight_change[1] <- sum(abs(weight_adj[1,]))
  
  change_pos <- matrix(rep(0,n),n)
  change_pos[1] <- TRUE
  
  
  for (i in 2:length(weight$V1)){
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
  # weight$DATETIME <- index(ret$test)[2:dim(ret$test)[1]]
  # weight[,12:15] <- as.data.frame(ret$test[2:dim(ret$test)[1],1:4])
  # names(weight)<-c('w1','w2','w3','w4','delta_v','change','w1_adj','w2_adj','w3_adj','w4_adj','delta_v_adj','ret1','ret2','ret3','ret4')
  
  return(list(weight=weight,weight_adj=weight_adj,weight_change=weight_change,change_pos=change_pos)) 
}

ret_analysis <- function(ret,freq='w'){
  ret <- as.xts(ret)
  n <- dim(ret)[1]
  if (freq=='w'){
    m <- 52
  } else if(freq=='m'){
    m <- 12
  } else{m <- 252}
  cum_ret <- cumprod(1+ret)
  vol <- sqrt(m)*sd(ret)
  mean_ret <- (cum_ret[n])^(m/n)-1
  drawdown <-  min(drawdowns(timeSeries(ret)))
  sharpe2 <- (mean_ret-0.02)/vol
  sharpe1 <- (mean_ret)/vol
  
  
  print('----------Return result----------')
  print(paste('Annual return: ',as.character(round(mean_ret,4)) ,sep=''))
  print(paste('Annual volatility: ',as.character(round(vol,4)) ,sep=''))
  print(paste('Sharpe ratio: ',as.character(round(sharpe1,3)) ,sep=''))
  print(paste('Maximum drawdown: ',as.character(round(drawdown,4)) ,sep=''))  
  print(paste('Sharpe ratio2: ',as.character(round(sharpe2,3)) ,sep=''))
  
  plot.zoo(xts(cum_ret,index(ret)),screens=c(1),col=c('blue'))
  
}

perf_analysis <- function(weight,asset_ret,cost=0.001,freq='w'){
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
    ret_after_fee[i] <- ret_before_fee[i]-cost*sum(abs(as.numeric(weight[i,])-as.numeric(weight[i-1,])))
    # print(cost*sum(abs(weight[i,]-weight[i-1,])))
    # print(abs(as.numeric(weight[i,])-as.numeric(weight[i-1,])))
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
  print(paste('Sharpe ratio1: ',as.character(round(sharpe_bf,3)) ,sep=''))
  print(paste('Sharpe ratio2: ',as.character(round(mean_ret_bf/vol_bf,3)),sep=''))
  
  print(paste('Maximum drawdown: ',as.character(round(drawdown_bf,4)) ,sep=''))
  
  # after fee
  print('----------After fee result----------')
  print(paste('Annual return: ',as.character(round(mean_ret_af,4)) ,sep=''))
  print(paste('Annual volatility: ',as.character(round(vol_af,4)) ,sep=''))
  print(paste('Sharpe ratio1: ',as.character(round(sharpe_af,3)) ,sep=''))
  print(paste('Sharpe ratio2: ',as.character(round(mean_ret_af/vol_af,3)),sep=''))
  print(paste('Sharpe ratio3: ',as.character(round((mean_ret_af-0.02)/vol_af,3)),sep=''))
  print(paste('Maximum drawdown: ',as.character(round(drawdown_af,4)) ,sep='')) 
  print(mean_ret_af/vol_af)
  
  
  print(summary(weight))
  dev.new()
  plot.zoo(xts(cbind(cum_ret_bf,cum_ret_af),index(asset_ret)),col=c('red','blue'),screens=c(1,1))
  dev.new()
  plot.zoo(xts(weight,index(asset_ret)))
  return(list(cum_ret_af=cum_ret_af,cum_ret_bf=cum_ret_bf,dd=drawdowns(timeSeries(ret_after_fee)),weight=weight))
}
