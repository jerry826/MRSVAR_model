# 

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

collect_data <- function(assets,start='2010-01-01',end='2016-08-16',mid='2015-01-01',freq='w'){
  # 获取资产收益率序列
  # assets: 资产代码序列
  # start:  训练集起始日期
  # end:    测试集结束日期
  # mid:    训练集、测试集间隔点
  # freq:   数据频率，w代表周，m代表月，d代表日
  for (i in 1:length(assets)){
    raw_data <- (w.wsd(assets[i],"pre_close,open,high,low,close",start,end))$Data
    if (i == 1){
      s <- raw_data$CLOSE
      time <- as.Date(raw_data$DATETIME)
    }else{
      s <- cbind(s,raw_data$CLOSE)
    }
  }
  dataset <- xts(s,time)
  # 数据重整
  if (freq=='w'){
    dataset <- apply.weekly(dataset,last)
  }else if (freq=='m'){ 
    dataset <- apply.monthly(dataset,last)
  }
  # 计算收益率
  ret <- ROC(dataset)
  ret <- ret[2:dim(ret)[1]]
  # 训练、测试集划分
  train <- ret[paste(start,mid,sep='/')]
  test <- ret[paste(mid,end,sep='/')]
  return(list(train=train,test=test) )
}

mrs_model <- function(ret,rlag=1,h=5){
  # 输入多元时间序列，训练模型
  # p: VAR模型滞后阶数
  # h: Markov隐状态数量
  
  model <- msvar(ret,p=rlag,h=h)
  var_beta <- model$hreg$Bk      # beta 5*4                                 VAR参数估计
  cov_error <- model$hreg$Sigmak # the covariance matrix of residuals 4*4   残差协方差矩阵
  error <- model$hreg$e          # errors 1607*4*5                          残差
  p <- model$fp                  # 1607*5                                   状态概率  
  q <- model$Q                   # 5*5                                      状态转移概率矩阵
  return(list(var_beta=var_beta,cov_error=cov_error,error=error,p=p,q=q,h=h,lag=rlag))
}

insample_mv_validation <- function(result,ret,h){
  # 样本内均值表现验证
  # result:
  #
  ret <- ret$train
  h <- result$h
  # 初始化
  weight <- c()
  ret_s <- c()
  # 
  n <- dim(ret)[1]
  
  
  for (i in (2:(n-1)))
  {
    ret_model <- 0
    weight_model <-c(0,0,0,0)
    
    p_next <- t(result$p[i,])%*%result$q
    for (j in (1:h))
    {
      ret_today <- ret[i,]  #1*4
      ret_next_est <- t(result$var_beta[,,j])%*%c(as.numeric(ret_today),1)
      S_est <- result$cov_error[,,j] # 4*4
      # print(ret_next_est)
      
      w <- weight_optimizer(S_est,ret_next_est)
      # w <- t(ret_next_est)%*%solve(S_est)
      # w <- w/sum(w)
      ret_model_state <- w%*%t(ret[i+1,])
      weight_model <- weight_model + p_next[j]*w
      
      ret_model <- ret_model + ret_model_state[1]*p_next[j]
    }
  
  weight <- rbind(weight,as.numeric(weight_model))  
  ret_s[i] <- ret_model
  }
  return(list(ret=ret_s,weight=weight))
}



cal_func_1 <- function(w,S,r){
  # 计算风险收益比
  # w:资产权重
  # S:资产收益率协方差
  # r:资产收益率均值
  p_r <- w%*%r
  p_std <- sqrt(t(w)%*%S%*%w)
  s <- p_r/p_std
  return(s)
}

cal_func_2 <- function(w,S,r){
  # 计算风险调整后的收益
  # w:资产权重
  # S:资产收益率协方差
  # r:资产收益率均值  
  p_r <- w%*%r
  p_var <- (t(w)%*%S%*%w)
  adj_ret <- p_r-0.5*p_var
  return(adj_ret)
}

cal_func_3 <-function(w,S){
  # 计算各资产对组合的风险贡献度
  # w:资产权重
  # S:资产收益率协方差
  n <- length(w)
  sigma_ij <- diag(w)%*%(S*10000)%*%w
  sigma_delta <- matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      sigma_delta[i,j] <- (sigma_ij[i] - sigma_ij[j])^2
    }
  }
  return(sum(sigma_delta)) 
}
  
f1 <- function(w){
  # 包装函数
  ans <- cal_func_1(w,S,r)
  return(-ans)
}

f2 <- function(w){
  # 包装函数
  ans <- cal_func_2(w,S,r)
  return(-ans)
}

f3 <- function(w){
  # 包装函数
  ans <- cal_func_3(w,S)
  return(ans)
}


weight_optimizer <- function(S0,r0,f){
  # 根据目标函数r做资产组合最优化
  # w:资产权重
  # S:资产收益率协方差
  # r:资产收益率均值 
  n <- length(r0)
  # 限制空头和杠杆
  par.l = rep(0,n); par.u = rep(1,n)
  # 限制总资产权重和为1
  A = matrix(rep(1,n),1)
  lin.l = c(1)
  lin.u = c(1)
  # 组合波动率限制
  nlcon1 = function(w){
    return(sqrt(w%*%S%*%w)*sqrt(252))
  }
  nlin.l = c(-Inf)  ; nlin.u = c(+Inf)  #目标年化波动率
  
  # 更改变量属性
  S <<- S0
  r <<- r0
  # 随机初始化参数
  init <- runif(n,min=0,max=1)
  # 模型优化
  m <- donlp2(init,f,par.u=par.u,par.l=par.l,A,lin.l=lin.l,lin.u=lin.u,
              nlin=  list(nlcon1) ,        #list(nlcon1,nlcon2), 
              nlin.u = nlin.u, nlin.l=nlin.l,control = donlp2Control())
  # 返回结果
  # print(nlcon1(m$par))
  return(m$par)
}

# weight_optimizer_3 <- function(S0){
#   n = dim(S0)[1]
#   par.l = rep(0,n)
#   par.u = rep(1,n)
#   A = matrix(rep(1,n),1)
#   lin.l = c(1)
#   lin.u = c(1)
#   
#   S <<- S0*10000
# 
#   
# 
#   
#   m = donlp2(runif(n,min=0,max=1), f3, par.u=par.u, par.l=par.l,
#                A,　lin.l=lin.l,lin.u=lin.u,  
#                nlin=  list(nlcon1) ,        #list(nlcon1,nlcon2), 
#                nlin.u = nlin.u, nlin.l=nlin.l,control = donlp2Control())
#   print(nlcon1(m$par))
#   
#   return(m$par)
# } 

sharpe <- function(ret,freq){
  # 计算年化夏普值
  # ret:收益率序列
  # freq:收益率序列频率
  if (freq=='w'){
    sharpe <- sqrt(52)*mean(ret)/sqrt(var(ret))
  }else if (freq=='m'){
    sharpe <- sqrt(12)*mean(ret)/sqrt(var(ret))
  }
  return(sharpe)
}

gx <- function(r,r_1,p,beta,sigma,h){
  # 贝叶斯概率更新函数
  # beta[i]代表第i个状态所有回归的系数，sigma[i]i个状态的协方差， r是今天的收益率 , r_1代表昨天的收益率
  # p是原来预测今天的概率
  dens = rep(0,h)
  for (i in 1:h) {
    dens[i]= dmvnorm( r , t(beta[,,i])%*%c(r_1,1), sigma[,,i])
  }
  return(p*dens/sum(p*dens))
}

mrs_predict <- function(result,test,train,h,type=1){
  # 样本外均值协方差矩阵预测，并进行最优权重优化 
  # result:mrs模型输出结果 
  # test  :测试集数据  
  # train :训练集数据
  # type  :权重最优化函数：f1分状态最优化风险收益比 f2 f3 
  h <- result$h
  m <- dim(test)[2]
  n <- dim(test)[1]
  p_list <- t(matrix(result$p[dim(result$p)[1],])) # 
  ret <- c() # 组合收益率
  weight <- c()
  beta <- result$var_beta # VAR系数
  sigma <- result$cov_error 
  # 在1到n-1期每期进行均值方差预测，
  for (i in 1:(n-1))
  { 
    # 当天资产收益率
    print(i)
    r <- test[i,]
    # 获取前一期的资产收益率，如果是训练集第一期，则前一期的数据来自于测试集的最后一期
    if (i==1){
      r_1 <- train[dim(train)[1],] # 训练集最后一期数据
    }else{
      r_1 <- test[i-1,]
    }
    
    # 计算当期状态概率
    p <- p_list[i,]%*%result$q  # q是转移概率矩阵
    # 贝叶斯更新当前状态概率
    p_update <- gx(as.numeric(r),as.numeric(r_1),p,beta,sigma,h) # 更新概率
    # 记录
    p_list <- rbind(p_list,p_update) # 记录当天的p
    # 通过状态转移矩阵预测下一期概率
    p_next <- (p_update)%*%result$q # 预测下一期的p
    
    # 组合权重最优化
    ret_model <- 0
    weight_model <- rep(0,m)
    # 模型1，2，3
    if (type==1){
    for (j in (1:h))
    { 
      # 在每个状态下进行均值方差模型预测
      ret_next_est <- t(beta[,,j])%*%c(as.numeric(r),1) # 根据VAR预测下一期收益率
      S_est <- result$cov_error[,,j]                    # 4*4 根据状态获取协方差矩阵
      w <- weight_optimizer(S_est,ret_next_est,f1)
      weight_model <- weight_model + p_next[j]*w
      
    }
    }else if (type==2){
      # 对每个状态的均值和方差进行概率加权
      ret_next_est <- t(beta[,,1])%*%c(as.numeric(r),1)*p_next[1]
      S_est <- result$cov_error[,,1]*p_next[1]                    
      
      for (j in (2:h))
      {
        ret_next_est <- ret_next_est + t(beta[,,j])%*%c(as.numeric(r),1)*p_next[j] 
        S_est <- S_est + result$cov_error[,,j]*p_next[j]                    
      }     
      # 均值方差预测
      weight_model <- weight_optimizer(S_est,ret_next_est,f1)
    }else if (type==3){
      # 对每个状态的协方差均值进行加权
      S_est <- result$cov_error[,,1]*p_next[1]
      ret_next_est <- rep(0,m)
      for (j in (2:h))
      {
        S_est <- S_est + result$cov_error[,,j]*p_next[j]                    
      }
      # 风险平价方法求最优权重
      weight_model <- weight_optimizer(S_est,ret_next_est,f3)
    
    }
    # 计算收益率
    ret_model <- (test[i+1,])%*%weight_model
    # 记录收益率和权重
    ret[i] <- ret_model
    weight <- rbind(weight,as.numeric(weight_model))
    
  }
  
  return(list(ret=ret,p_list=p_list,weight=weight))
}

adjust_weight <- function(perf,qu=0.05,type=1){
  # 调仓方式优化
  # perf:mrs mv模型优化结果
  # qu  :阈值
  # type:调整方法:1根据最大调仓量;2根据总调仓量
  n <- dim(perf$weight)[1]
  m <- dim(perf$weight)[2]
  
  weight <- as.data.frame(perf$weight)
  
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
    ret_after_fee[i] <- ret_before_fee[i]-cost*sum(abs(weight[i,]-weight[i-1,]))
    
  }
  
  cum_ret_bf <- cumsum(ret_before_fee)
  cum_ret_af <- cumsum(ret_after_fee)
  
  
  mean_ret_bf <- mean(ret_before_fee)*m
  mean_ret_af <- mean(ret_after_fee)*m 
  sharpe_bf <- sharpe(ret_before_fee,freq)
  sharpe_af <- sharpe(ret_after_fee,freq)
  print(paste('Sharpe before fee: ',as.character(sharpe_bf) ,sep=''))
  print(paste('Sharpe after fee: ',as.character(sharpe_af) ,sep=''))
  print(paste('Annual return before fee: ',as.character(mean_ret_bf) ,sep=''))
  print(paste('Annual return after fee: ',as.character(mean_ret_af) ,sep='')) 
  
  print(summary(weight))
  plot.zoo(xts(cbind(cum_ret_bf,cum_ret_af),index(asset_ret)),screens=c(1,1),col=c('red','blue'))
  plot.zoo(xts(weight2$weight_adj,index(asset_ret)))
}


covar <- function(A, h){  #A 收益率矩阵
  #h = 60  #计算协方差矩阵时用到收益的长度
  n = dim(A)[1]
  covar_2 = c()
  for (i in 1:(n-h+1)) {
    a = cov(A[i:(i+h-1),])
    covar_2 = c(covar_2, list(a))
  }
  return(covar_2)
}

weight_cal <- function(h){
  n=length(h)
  weight = matrix(0, n, dim(h[[1]])[1])
  last_weight <- rep(0, dim(h[[1]])[1])
  for (i in 1:n) {
    print('-----------------------------------------------------------')
    if (i > 1){
      opt_weight = weight_optimizer_3(h[[i]])
      while (sum(abs(opt_weight-last_weight))>0.30){
      opt_weight = weight_optimizer_3(h[[i]])
      print(i)
      print(opt_weight)
      print(last_weight)
      print(sum(abs(opt_weight-last_weight)))
      }
      weight[i,] <- opt_weight
      last_weight <- opt_weight
    }else{
      opt_weight <- weight_optimizer_3(h[[i]])
      weight[i,] <- opt_weight
      last_weight <- opt_weight
      
    }
    
  }
  return(weight)
}

# # h = weight_cal(covar(ret, 60))
# cov_long=c()
# for (j in c(1,4,5,6)) {
#   
#   
#   S = result$cov_error[,,j]
#   S_2 =  result$cov_error[,,j]
#   A = t(result$var_beta[,,j][1:4,])
#   for (i in 1:1000) {
#     S_1 = A%*%S%*%t(A)+S_2
#     S = S_1
#     #print(S)
#   }
#   cov_long = c(cov_long, list(S))
# }

