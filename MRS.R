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
  # ��ȡ�ʲ�����������
  # assets: �ʲ���������
  # start:  ѵ������ʼ����
  # end:    ���Լ���������
  # mid:    ѵ���������Լ������
  # freq:   ����Ƶ�ʣ�w�����ܣ�m�����£�d������
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
  # ��������
  if (freq=='w'){
    dataset <- apply.weekly(dataset,last)
  }else if (freq=='m'){ 
    dataset <- apply.monthly(dataset,last)
  }
  # ����������
  ret <- ROC(dataset)
  ret <- ret[2:dim(ret)[1]]
  # ѵ�������Լ�����
  train <- ret[paste(start,mid,sep='/')]
  test <- ret[paste(mid,end,sep='/')]
  return(list(train=train,test=test) )
}

mrs_model <- function(ret,rlag=1,h=5){
  # �����Ԫʱ�����У�ѵ��ģ��
  # p: VARģ���ͺ����
  # h: Markov��״̬����
  
  model <- msvar(ret,p=rlag,h=h)
  var_beta <- model$hreg$Bk      # beta 5*4                                 VAR��������
  cov_error <- model$hreg$Sigmak # the covariance matrix of residuals 4*4   �в�Э�������
  error <- model$hreg$e          # errors 1607*4*5                          �в�
  p <- model$fp                  # 1607*5                                   ״̬����  
  q <- model$Q                   # 5*5                                      ״̬ת�Ƹ��ʾ���
  return(list(var_beta=var_beta,cov_error=cov_error,error=error,p=p,q=q,h=h,lag=rlag))
}

insample_mv_validation <- function(result,ret,h){
  # �����ھ�ֵ������֤
  # result:
  #
  ret <- ret$train
  h <- result$h
  # ��ʼ��
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
  # ������������
  # w:�ʲ�Ȩ��
  # S:�ʲ�������Э����
  # r:�ʲ������ʾ�ֵ
  p_r <- w%*%r
  p_std <- sqrt(t(w)%*%S%*%w)
  s <- p_r/p_std
  return(s)
}

cal_func_2 <- function(w,S,r){
  # ������յ����������
  # w:�ʲ�Ȩ��
  # S:�ʲ�������Э����
  # r:�ʲ������ʾ�ֵ  
  p_r <- w%*%r
  p_var <- (t(w)%*%S%*%w)
  adj_ret <- p_r-0.5*p_var
  return(adj_ret)
}

cal_func_3 <-function(w,S){
  # ������ʲ�����ϵķ��չ��׶�
  # w:�ʲ�Ȩ��
  # S:�ʲ�������Э����
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
  # ��װ����
  ans <- cal_func_1(w,S,r)
  return(-ans)
}

f2 <- function(w){
  # ��װ����
  ans <- cal_func_2(w,S,r)
  return(-ans)
}

f3 <- function(w){
  # ��װ����
  ans <- cal_func_3(w,S)
  return(ans)
}


weight_optimizer <- function(S0,r0,f){
  # ����Ŀ�꺯��r���ʲ�������Ż�
  # w:�ʲ�Ȩ��
  # S:�ʲ�������Э����
  # r:�ʲ������ʾ�ֵ 
  n <- length(r0)
  # ���ƿ�ͷ�͸ܸ�
  par.l = rep(0,n); par.u = rep(1,n)
  # �������ʲ�Ȩ�غ�Ϊ1
  A = matrix(rep(1,n),1)
  lin.l = c(1)
  lin.u = c(1)
  # ��ϲ���������
  nlcon1 = function(w){
    return(sqrt(w%*%S%*%w)*sqrt(252))
  }
  nlin.l = c(-Inf)  ; nlin.u = c(+Inf)  #Ŀ���껯������
  
  # ���ı�������
  S <<- S0
  r <<- r0
  # �����ʼ������
  init <- runif(n,min=0,max=1)
  # ģ���Ż�
  m <- donlp2(init,f,par.u=par.u,par.l=par.l,A,lin.l=lin.l,lin.u=lin.u,
              nlin=  list(nlcon1) ,        #list(nlcon1,nlcon2), 
              nlin.u = nlin.u, nlin.l=nlin.l,control = donlp2Control())
  # ���ؽ��
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
#                A,��lin.l=lin.l,lin.u=lin.u,  
#                nlin=  list(nlcon1) ,        #list(nlcon1,nlcon2), 
#                nlin.u = nlin.u, nlin.l=nlin.l,control = donlp2Control())
#   print(nlcon1(m$par))
#   
#   return(m$par)
# } 

sharpe <- function(ret,freq){
  # �����껯����ֵ
  # ret:����������
  # freq:����������Ƶ��
  if (freq=='w'){
    sharpe <- sqrt(52)*mean(ret)/sqrt(var(ret))
  }else if (freq=='m'){
    sharpe <- sqrt(12)*mean(ret)/sqrt(var(ret))
  }
  return(sharpe)
}

gx <- function(r,r_1,p,beta,sigma,h){
  # ��Ҷ˹���ʸ��º���
  # beta[i]������i��״̬���лع��ϵ����sigma[i]i��״̬��Э��� r�ǽ���������� , r_1���������������
  # p��ԭ��Ԥ�����ĸ���
  dens = rep(0,h)
  for (i in 1:h) {
    dens[i]= dmvnorm( r , t(beta[,,i])%*%c(r_1,1), sigma[,,i])
  }
  return(p*dens/sum(p*dens))
}

mrs_predict <- function(result,test,train,h,type=1){
  # �������ֵЭ�������Ԥ�⣬����������Ȩ���Ż� 
  # result:mrsģ�������� 
  # test  :���Լ�����  
  # train :ѵ��������
  # type  :Ȩ�����Ż�������f1��״̬���Ż���������� f2 f3 
  h <- result$h
  m <- dim(test)[2]
  n <- dim(test)[1]
  p_list <- t(matrix(result$p[dim(result$p)[1],])) # 
  ret <- c() # ���������
  weight <- c()
  beta <- result$var_beta # VARϵ��
  sigma <- result$cov_error 
  # ��1��n-1��ÿ�ڽ��о�ֵ����Ԥ�⣬
  for (i in 1:(n-1))
  { 
    # �����ʲ�������
    print(i)
    r <- test[i,]
    # ��ȡǰһ�ڵ��ʲ������ʣ������ѵ������һ�ڣ���ǰһ�ڵ����������ڲ��Լ������һ��
    if (i==1){
      r_1 <- train[dim(train)[1],] # ѵ�������һ������
    }else{
      r_1 <- test[i-1,]
    }
    
    # ���㵱��״̬����
    p <- p_list[i,]%*%result$q  # q��ת�Ƹ��ʾ���
    # ��Ҷ˹���µ�ǰ״̬����
    p_update <- gx(as.numeric(r),as.numeric(r_1),p,beta,sigma,h) # ���¸���
    # ��¼
    p_list <- rbind(p_list,p_update) # ��¼�����p
    # ͨ��״̬ת�ƾ���Ԥ����һ�ڸ���
    p_next <- (p_update)%*%result$q # Ԥ����һ�ڵ�p
    
    # ���Ȩ�����Ż�
    ret_model <- 0
    weight_model <- rep(0,m)
    # ģ��1��2��3
    if (type==1){
    for (j in (1:h))
    { 
      # ��ÿ��״̬�½��о�ֵ����ģ��Ԥ��
      ret_next_est <- t(beta[,,j])%*%c(as.numeric(r),1) # ����VARԤ����һ��������
      S_est <- result$cov_error[,,j]                    # 4*4 ����״̬��ȡЭ�������
      w <- weight_optimizer(S_est,ret_next_est,f1)
      weight_model <- weight_model + p_next[j]*w
      
    }
    }else if (type==2){
      # ��ÿ��״̬�ľ�ֵ�ͷ�����и��ʼ�Ȩ
      ret_next_est <- t(beta[,,1])%*%c(as.numeric(r),1)*p_next[1]
      S_est <- result$cov_error[,,1]*p_next[1]                    
      
      for (j in (2:h))
      {
        ret_next_est <- ret_next_est + t(beta[,,j])%*%c(as.numeric(r),1)*p_next[j] 
        S_est <- S_est + result$cov_error[,,j]*p_next[j]                    
      }     
      # ��ֵ����Ԥ��
      weight_model <- weight_optimizer(S_est,ret_next_est,f1)
    }else if (type==3){
      # ��ÿ��״̬��Э�����ֵ���м�Ȩ
      S_est <- result$cov_error[,,1]*p_next[1]
      ret_next_est <- rep(0,m)
      for (j in (2:h))
      {
        S_est <- S_est + result$cov_error[,,j]*p_next[j]                    
      }
      # ����ƽ�۷���������Ȩ��
      weight_model <- weight_optimizer(S_est,ret_next_est,f3)
    
    }
    # ����������
    ret_model <- (test[i+1,])%*%weight_model
    # ��¼�����ʺ�Ȩ��
    ret[i] <- ret_model
    weight <- rbind(weight,as.numeric(weight_model))
    
  }
  
  return(list(ret=ret,p_list=p_list,weight=weight))
}

adjust_weight <- function(perf,qu=0.05,type=1){
  # ���ַ�ʽ�Ż�
  # perf:mrs mvģ���Ż����
  # qu  :��ֵ
  # type:��������:1������������;2�����ܵ�����
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
      # ������ֵ,����
      change_pos[i] <- TRUE
      weight_adj[i,] <- as.numeric(weight[i,])
    }else{
      # δ������ֵ,������
      change_pos[i] <- FALSE
      weight_adj[i,] <- as.numeric(weight_adj[i-1,])
    }
    # �������ÿ���ʲ��ֱֲ仯��
    delta_w1 <-  (sum(abs(weight_adj[i,]-weight_adj[i-1,])))
    weight_change[i] <- delta_w1
  }
  # weight$DATETIME <- index(ret$test)[2:dim(ret$test)[1]]
  # weight[,12:15] <- as.data.frame(ret$test[2:dim(ret$test)[1],1:4])
  # names(weight)<-c('w1','w2','w3','w4','delta_v','change','w1_adj','w2_adj','w3_adj','w4_adj','delta_v_adj','ret1','ret2','ret3','ret4')

  return(list(weight=weight,weight_adj=weight_adj,weight_change=weight_change,change_pos=change_pos)) 
}





perf_analysis <- function(weight,asset_ret,cost=0.001,freq='w'){
  # ���Ա��ַ���
  # weight   :Ȩ��list
  # asset_ret:��Ӧ�ʲ�������
  # cost     :���߽��������ѷ�
  # Ȩ�����г���
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


covar <- function(A, h){  #A �����ʾ���
  #h = 60  #����Э�������ʱ�õ�����ĳ���
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
