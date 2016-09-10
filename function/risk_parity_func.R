# author: Jiayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: risk_parity.R
# time: 2016/9/2


risk_parity2 <- function(id,date='2016-08-26',u_limit=rep(1,9),l_limit=rep(0,9),u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,window=60,init=c(0)){
  window_str <- paste('ED-',as.character(window+10),'TD',sep='')
  ret <- get_data2(start=window_str,end=date,currency=currency)
  ret <- ret[,id]
  cov <- ew_var(ret)
  if (sum(init)!=1){
    init<-runif(dim(cov)[1],min=0,max=1)
  }
  weight <- optimize_weight1(init,cov,u_limit,l_limit,u_sigma,l_sigma)
  return(weight)
}

risk_parity3 <- function(id,risk_weight,date='2016-08-26',u_limit=rep(1,9),l_limit=rep(0,9),u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,path=NA,window=60,init=c(0)){
  window_str <- paste('ED-',as.character(window+10),'TD',sep='')
  # 两种数据读取方法
  ## Wind获取
  if (is.na(path)){
  ret <- get_data2(start=window_str,end=date,currency=currency)
  }else{
  ## 本地读取
  dataset <- read.csv(path)
  ret <- xts(dataset[,2:dim(dataset)[2]],as.Date(dataset[,1]))
  end <- (w.tdaysoffset(-1,date))$Data[,1]
  start <- (w.tdaysoffset(-window,date))$Data[,1]  
  ret <- ret[paste(as.character(start),'/',as.character(end),sep=''),id]
  }
  # 计算协方差
  cov <- ew_var(ret)
  if (sum(init)!=1){
    init<-runif(dim(cov)[1],min=0,max=1)
  }
  opt_result <- optimize_weight2(init,risk_weight,cov,u_limit,l_limit,u_sigma,l_sigma)
  opt_result$ret <-  apply(ret*t((as.matrix(opt_result$weight)%*%matrix(1,1,dim(ret)[1]))),1,sum)
  return(opt_result)
}

risk_parity4 <- function(ret,risk_weight,u_limit=rep(1,9),l_limit=rep(0,9),u_sigma=+Inf,l_sigma=-Inf,init=c(0)){
  # 计算协方差
  cov <- ew_var(ret)
  if (sum(init)!=1){
    init<-runif(dim(cov)[1],min=0,max=1)
  }
  opt_result <- optimize_weight2(init,risk_weight,cov,u_limit,l_limit,u_sigma,l_sigma)
  opt_result$ret <-  apply(ret*t((as.matrix(opt_result$weight)%*%matrix(1,1,dim(ret)[1]))),1,sum)
  return(opt_result)
}
# 新数据
# 1.  H11137.CSI      中国互联网(美元)             Equity    2007-07-01
# 2.  HSCEI.HI        恒生国企指数(港币)           Equity    2005-04-14
# 3.  000016.SH       上证50                       Equity    2004-01-01
# 4.  399006.SZ       创业板指                     Equity    2010-06-02
# 5.  000905.SH       中证500                      Equity    2005-01-04
# 6.  000300.SH       沪深300                      Equity    2002

# 14.  SPX.GI          标普500                     Equity    1966
# 15.  IXIC.GI         纳斯达克指数                Equity    1971
# 16.  GDAXI.GI        德国DAX                     Equity    1965

# 7.  H11017.CSI      中期国债                     Bond      2008-01-01
# 8.  H11078.CSI      中证中高信用                 Bond      2008-01-01
# 9.  0401E.CS        中债固定利率全价(7-10)指数   Bond      2006-11-20
# 17  070025.OF       嘉实信用A                    Bond      2011-09-14
# 10. 037.CS          中债总财富指数               Bond      2002-01-01

# 18. 000198.OF       天弘余额宝                   Cash      2013-06-03
# 19. 482002.OF       工银货币                     Cash      2006-03-26

# 11. CL.NYM NYMEX    原油(美元)                   Commodity 2011-12-21
# 12. AU9999.SGE      黄金                         Commodity 2004-01-01
# 13. CCFI.WI         WIND商品                     Commodity 2002-01-01
get_data2 <- function(start='2016-08-10',end='2016-08-26',currency=FALSE,type='continuous'){
  w.start()
  # 采用收盘价部分
  sid1 <- "H11137.CSI,HSCEI.HI,000016.SH,399006.SZ,000905.SH,000300.SH,H11017.CSI,H11078.CSI,0401E.CS,037.CS,CL.NYM,AU9999.SGE,CCFI.WI,SPX.GI,IXIC.GI,GDAXI.GI"
  if (currency) {
    data1 <- (w.wsd(sid1,"close",start,end,"Fill=Previous;Currency=CNY"))$Data
  }else{
    data1 <- (w.wsd(sid1,"close",start,end,"Fill=Previous"))$Data
  }
  # 采用万份净值部分
  sid2 <- "000198.OF,482002.OF"
  data2 <- (w.wsd(sid2,"mmf_unityield",start,end,"Fill=Previous"))$Data
  names(data2) <- c('TIME',sid2)
  # 采用基金净值部分
  sid3 <- "070025.OF"
  data3 <- (w.wsd(sid3,"nav",start,end,"Fill=Previous"))$Data
  # names(data3) <- c('TIME',sid3)
  # 数据整合
  dataset <- cbind(data1,data3[,2:dim(data3)[2]])
  names(dataset) <- c('TIME',unlist(strsplit(sid1,split=',')),unlist(strsplit(sid3,split=',')))
  # 计算收益率
  ret <- ROC(xts(dataset[,2:dim(dataset)[2]],dataset[,1]),type=c(type))
  
  # names(data2) <- c('TIME',unlist(strsplit(sid2,split=',')))
  temp <- xts(data2[,2:dim(data2)[2]],data2[,1])*365/(10000*252)
  ret2 <- merge(ret,temp)
  ret2 <- ret2[2:dim(ret2)[1],]
  names(ret2) <- c(unlist(strsplit(sid1,split=',')),unlist(strsplit(sid3,split=',')),unlist(strsplit(sid2,split=',')))
  return(ret2)
}


risk_parity <- function(date='2016-08-26',u_limit=rep(1,9),l_limit=rep(0,9),u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,window=60,init=c(0)){
  window_str <- paste('ED-',as.character(window+10),'TD',sep='')
  ret <- get_data(start=window_str,end=date,currency=currency)
  ret <- ret[11:dim(ret)[1]]
  cov <- ew_var(ret)
  if (sum(init)!=1){
    init<-runif(dim(cov)[1],min=0,max=1)
  }
  weight <- optimize_weight1(init,cov,u_limit,l_limit,u_sigma,l_sigma)
  return(weight)
}

risk_budget <- function(date='2016-08-26',u_limit=rep(1,9),l_limit=rep(0,9),u_sigma=+Inf,l_sigma=-Inf,currency=FALSE,window=60,init=c(0),rw=c(0)){
  window_str <- paste('ED-',as.character(window+10),'TD',sep='')
  ret <- get_data(start=window_str,end=date,currency=currency)
  ret <- ret[11:dim(ret)[1]]
  cov <- ew_var(ret)
  init<- runif(dim(cov)[1],min=0,max=1)
  weight <- optimize_weight2(init,rw,cov,u_limit,l_limit,u_sigma,l_sigma)
  return(weight)
}

asset_vol <- function(id,window_list=c(60),date='2016-08-26',path=NA){
  for (i in 1:length(window_list)){
  window <- window_list[i]
  # 两种数据读取方法
  ## Wind获取
  if (is.na(path)){
    window_str <- paste('ED-',as.character(window+10),'TD',sep='')
    ret <- get_data2(start=window_str,end=date,currency=currency)
  }else{
    ## 本地读取
    dataset <- read.csv(path)
    ret <- xts(dataset[,2:dim(dataset)[2]],as.Date(dataset[,1]))
    end <- (w.tdaysoffset(-1,date))$Data[,1]
    start <- (w.tdaysoffset(-window,date))$Data[,1]  
    ret <- ret[paste(as.character(start),'/',as.character(end),sep=''),id]
  }
  
  # 计算协方差
  print(paste('Annualized asset volatility for ',as.character(window),' days (%):',sep=''))
  cov <- ew_var(ret)
  print(as.numeric(round(sqrt(diag(cov)*10000*252),2)))
  }
}


get_data <- function(start='2016-08-10',end='2016-08-26',currency=FALSE,type='continuous'){
  # H11137.CSI      中国互联网(美元)
  # HSCEI.HI        恒生国企指数(港币)
  # 000198.OF       天弘余额宝
  # H11017.CSI      中期国债
  # CL.NYM NYMEX    原油(美元)
  # AU9999.SGE      黄金
  # 070025.OF       嘉实信用A
  # 000016.SH       上证50
  # 399006.SZ       创业板指
  # 采用收盘价部分
  sid1 <- "H11137.CSI,HSCEI.HI,H11017.CSI,CL.NYM,AU9999.SGE,000016.SH,399006.SZ"
  if (currency) {
    data1 <- (w.wsd(sid1,"close",start,end,"Fill=Previous;Currency=CNY"))$Data
  }else{
    data1 <- (w.wsd(sid1,"close",start,end,"Fill=Previous"))$Data
  }
  # 采用万份净值部分
  sid2 <- "000198.OF"
  data2 <- (w.wsd(sid2,"mmf_unityield",start,end,"Fill=Previous"))$Data
  names(data2) <- c('TIME',sid2)
  # 采用基金净值部分
  sid3 <- "070025.OF"
  data3 <- (w.wsd(sid3,"nav",start,end,"Fill=Previous"))$Data
  names(data3) <- c('TIME',sid3)
  # 数据整合
  dataset <- cbind(data1,data3[,2:dim(data3)[2]])
  names(dataset) <- c('TIME',unlist(strsplit(sid1,split=',')),unlist(strsplit(sid3,split=',')))
  # 计算收益率
  ret <- ROC(xts(dataset[,2:dim(dataset)[2]],dataset[,1]),type=c(type))
  
  names(data2) <- c('TIME',unlist(strsplit(sid2,split=',')))
  temp <- xts(data2[,2],data2[,1])/10000
  ret2 <- merge(ret,temp)
  ret2 <- ret2[2:dim(ret2)[1],]
  return(ret2)
}

ew_var <- function(X){
  m <- dim(X)[1]
  n <- dim(X)[2]
  # deman
  X_demean <- X-matrix(1,m,1)%*%(apply(X,2,mean))
  V <- array(rep(0,m*n*n),c(m,n,n)) 
  ewvar <- array(rep(0,n*n),c(n,n))
  for (i in 1:m){
    x <- as.matrix(X_demean[i,]) 
    S <- t(x)%*%x  
    V[i,,] <- S
  } 
  for (j in (1:n)){
    for (k in (1:n)){
      v <- V[,j,k]
      ewvar[j,k] <- WMA(v,n=m)[m]
    }
  }
  return(ewvar)
}

optimize_weight1 <- function(init,cov,u_limit,l_limit,u_sigma=+Inf,l_sigma=(-Inf)){
  n <- dim(cov)[1]
  # 设置空头和杠杆
  par.l <- l_limit
  par.u <- u_limit
  # 限制总资产权重和为1
  A = matrix(rep(1,n),1)
  lin.l <- c(1)
  lin.u <- c(1)
  # 组合波动率限制
  nlcon1 = function(w){
    return(sqrt((w%*%S*10000)%*%w)*sqrt(252)) #年化波动率
  }
  nlin.l = c(l_sigma)  ; nlin.u = c(u_sigma)  #目标年化波动率
  
  # 更改变量属性
  S <<- cov
  # 随机初始化参数
  
  # 模型优化
  m <- donlp2(init,f1,par.u=par.u,par.l=par.l,A,lin.l=lin.l,lin.u=lin.u,
              nlin=  list(nlcon1) ,        #list(nlcon1,nlcon2), 
              nlin.u = nlin.u, nlin.l=nlin.l,control = donlp2Control())
  risk_weight <- rr(m$par,S)
  # 返回结果
  print(paste('Volatility returned : ', as.character(nlcon1(m$par))))
  print(paste('Target function value : ',as.character(f1(m$par))))
  return(list(weight=m$par,risk_weight=rr(m$par,S),error=f1(m$par)))
}

optimize_weight2 <- function(init,risk_weight,cov,u_limit,l_limit,u_sigma=+Inf,l_sigma=(-Inf)){
  n <- dim(cov)[1]
  # 设置空头和杠杆
  par.l <- l_limit
  par.u <- u_limit
  # 限制总资产权重和为1
  A = matrix(rep(1,n),1)
  lin.l <- c(1)
  lin.u <- c(1)
  # 组合波动率限制
  nlcon1 = function(w){
    return(sqrt((w%*%S*10000)%*%w)*sqrt(252)) #年化波动率
  }
  nlin.l = c(l_sigma)  ; nlin.u = c(u_sigma)  #目标年化波动率
  # 更改变量属性
  S <<- cov
  rw <<- risk_weight
  # 随机初始化参数
  asset_vol <- as.numeric(round(sqrt(diag(S)*10000*252),4))
  
  # 模型优化
  m <- donlp2(init,f2,par.u=par.u,par.l=par.l,A,lin.l=lin.l,lin.u=lin.u,
              nlin=  list(nlcon1) ,        #list(nlcon1,nlcon2), 
              nlin.u = nlin.u, nlin.l=nlin.l,control = donlp2Control())
  # 返回结果
  print(paste('Volatility returned : ', as.character(round(nlcon1(m$par),3))))
  # print(paste('Target function value : ',as.character(round(f2(m$par),3))))
  print('Asset weight allocation :')
  
  print(round(m$par,2))
  print('Risk weight returned : ' )
  print(round(as.numeric(rr(m$par,S)),2))
  print('Asset volatility (%): ')
  print(asset_vol)
  print(paste('ERROR : ',as.character(f2(m$par)),sep=''))
  return(list(weight=m$par,risk_weight=rr(m$par,S),error=f2(m$par),asset_vol=asset_vol))
}



risk_contribution <-function(w,S){
  # 计算各资产对组合的风险贡献度
  # w:资产权重
  # S:资产收益率协方差
  n <- length(w)
  sigma_ij <- diag(w)%*%(S*10000)%*%w
  sigma_delta <-  (sigma_ij%*%matrix(rep(1,n),1)-t(sigma_ij%*%matrix(rep(1,n),1)))^2
  return(sum(sigma_delta)) 
}

rr <- function(w,S){
  sigma_ij <- diag(w)%*%(S*10000)%*%w
  sigma <- w%*%(S*10000)%*%w
  # print(paste('Sigma : ',as.character(sigma),sep=''))
  weight <- sigma_ij/as.numeric(sigma)
  return(t(weight))
}

risk_budget_error <- function(w,S,rw){
  sigma_ij <- diag(w)%*%(S*10000)%*%w
  sigma <- w%*%(S*10000)%*%w
  # weight <- sigma_ij/as.numeric(sigma)
  error <- 10000*sum((sigma_ij-rw*sigma)^2)
  return(error)
}


f1 <- function(w){
  # 包装函数1
  ans <- risk_contribution(w,S)
  return(ans)
}

f2 <- function(w){
  # 包装函数2
  error <- risk_budget_error(w,S,rw)
  return(error)
}

test <- function(){
  # 1. H11137.CSI      中国互联网(美元)
  # 2. HSCEI.HI        恒生国企指数(港币)
  # 3. 000198.OF       天弘余额宝
  # 4. H11017.CSI      中期国债
  # 5. CL.NYM NYMEX    原油(美元)
  # 6. AU9999.SGE      黄金
  # 7. 070025.OF       嘉实信用A
  # 8. 000016.SH       上证50
  # 9. 399006.SZ       创业板指
  # 参数设置
  date <- '2013-10-08'
  u_limit <- c(1,1,1,1,1,1,1,1,1)  # 资产权重上限 
  l_limit <- rep(0,9)              # 资产权重下限
  u_sigma <- +Inf                  # 年化波动率上限
  l_sigma <- -Inf                     # 年化波动率下限
  currency <- FALSE                # 是否进行货币转换
  
  # 
  w <- risk_budget(date=date,u_limit=u_limit,l_limit=l_limit,u_sigma=u_sigma,l_sigma=l_sigma,currency=currency,rw=rep(1/9,9))
  return(w)
}





