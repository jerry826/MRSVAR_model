# author: iayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: MRS.R
# time: 2016/9/22
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ========  Markov Regime Switch Model  ===============================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#初始化

##### MRS #####
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(WindR)
library(MSBVAR)
library(xts)
library(TTR)
library(reshape2)

# 计算MRS模型
mrs_model <- function(ret,rlag=1,h=5,niterblkopt=20,disp=TRUE){
  # MRS 模型训练
  # 输入
  # rlag        : VAR模型滞后阶数
  # h           : Markov隐状态数量
  # niterblkopt :算法优化的迭代次数
  # disp        :是否打印主要结果
  # 返回
  # var_beta    :VAR模型参数估计结果
  # cov_error   :误差项协方差矩阵估计结果
  # error       :残差
  # regime_prob :各状态概率
  # q           :概率转移矩阵
  # mu          :各状态条件均值
  # model       :模型估计结果
  # ret_next_est:下一期收益率预测
  # S_next_est  :协方差矩阵估计
  model <- msvar(ret,p=rlag,h=h,niterblkopt=niterblkopt)
  n <- dim(ret)[1]
  m <- dim(ret)[2]
  var_beta <- model$hreg$Bk                                          # VAR参数估计
  cov_error <- model$hreg$Sigmak                                     # 残差协方差矩阵
  error <- model$hreg$e          #                                   # 残差
  regime_prob <- xts(model$fp,index(ret)[2:n])                       # 状态概率  
  q <- model$Q                                                       # 状态转移概率矩阵
  ##     输出各状态下资产的均值和方差    ####
  mu = matrix(0,m,h)
  for (i in 1:h) {
    aa =t(var_beta[,,i])
    mu[,i] = solve(diag(1,m)-aa[,1:m],aa[,m+1])
  }
  ## 预测下一期收益率
  p_last <- regime_prob[dim(regime_prob)[1],]
  p_next <- p_last%*%q
  ret_last <- ret[dim(ret)[1],]
  ret_next_est <- rep(0,m)
  S_next_est <- matrix(0,m,m)
  for (j in (1:h))
  {
    ret_next_est <- ret_next_est + (t(var_beta[,,j])%*%c(as.numeric(ret_last),1))*p_next[j]
    print( ret_next_est)
    print(t(var_beta[,,j])%*%c(as.numeric(ret_last),1)*p_next[j])
    S_next_est <- S_next_est+cov_error[,,j]
    
  }
  if (disp) {
    print('1: Beta estimators')
    print(var_beta)
    print('2: Conditional mean return')
    print(mu)
  }
  return(list(var_beta=var_beta,cov_error=cov_error,error=error,regime_prob=regime_prob,q=q,model=model,mu=mu,ret_next_est=ret_next_est,S_next_est=S_next_est))
}
# 分段计算
mrs_output <- function(ret,mrs_regime){
  n0 <- dim(mrs_regime)[1]  # 长度
  asset_num <- dim(ret)[2]  # 资产数量
  circle <- matrix(0,n0,1)  # 周期记录
  # 循环
  k <- 1
  for (i in 1:n0){
    if (i==1){
      circle[i] <- k
    }else{
      if (as.numeric(mrs_regime[i])==as.numeric(mrs_regime[i-1])){
        circle[i] <- k
      }else{
        k <- k+1   # 周期变更
        circle[i] <- k
      }
    }  
  }
  n1 <- circle[n0] # 周期总数
  # 输出结果
  output <- data.frame(matrix(0,n1,asset_num+5))
  for (i in 1:n1){
    id1 <- (circle==i)
    output[i,1:asset_num] <- round((apply(ret[id1],2,mean)),4)       # 周期内资产收益率均值
    output[i,asset_num+1] <- mean(mrs_regime[id1])                   # 周期所处状态
    output[i,asset_num+2] <- length(mrs_regime[id1])                 # 周期长度
    output[i,asset_num+3] <- (index(ret)[id1])[1]                    # 周期起始月份
    output[i,asset_num+4] <- (index(ret)[id1])[output[i,asset_num+2]]# 周期结束月份
    output[i,asset_num+5] <- i                                       # 周期结束月份
    
  }
  output[,asset_num+3] <- as.Date(output[,asset_num+3])
  output[,asset_num+4] <- as.Date(output[,asset_num+4])
  names(output)[(asset_num+1):(asset_num+5)] <- c('Regime','Length','Start','End','id')
  print('Regime summary')
  print(output)
  return(output)
}
# 多资产作图
multi_regime_plot <- function(regime,cum_ret,start='',end=''){
  m <- dim(cum_ret)[2]
  regime <- regime[paste(start,'/',end,sep='')]
  cum_ret = cum_ret[paste(start,'/',end,sep='')]
  regime <- factor(regime)
  time <- index(cum_ret)
  pl <- lapply(1:m,function(.x)  regime_plot(data.frame(time,cum_ret[,.x]),regime)     )
  ml <- marrangeGrob(pl, nrow=m, ncol=1)
  return(ml)
}
# 单资产作图
regime_plot <- function(data,regime){
  names(data)=c('time','cum_ret')
  fig <- ggplot(data)
  fig <- fig+theme_igray()
  fig <- fig+geom_line(aes(x=time,y=cum_ret),color='grey',size=1)
  fig <- fig+geom_point(mapping=aes(x=time,y=cum_ret,colour=regime),size=3,pch=16)
  return(fig)
}

w.start()

####################################
#####       MRS计算并作图      #####
####################################
# 获取收益率数据
start <- '2006-01-01'
end <- '2016-08-30'
# 输入资产数量
num <- 4 
raw_data1 <- (w.wsd('000300.SH,000905.SH,037.CS,AU9999.SGE',"close",start,end,"Period=M"))$Data
data1 <- xts(raw_data1[,2:(num+1)],raw_data1[,1])
# 计算收益率
ret1 <- ROC(data1,type=c('continuous'),na.pad=FALSE)
#  计算累计收益
cum_ret <- cumprod(ret1+1)[2:dim(ret1)[1]]
# 产生MRS regimes
## 训练模型
result <- mrs_model(ret1,rlag=1,h=4,niterblkopt=20)
## 各状态作图
#  作图
mrs_regime <- xts(apply(result$regime_prob,1,which.max),index(result$regime_prob))
multi_regime_plot(mrs_regime, cum_ret,start='2006-01-01',end='2011-01-01')
multi_regime_plot(mrs_regime, cum_ret)
# 获取下一期收益率估计
ret_next_est <- result$ret_next_est


####################################
#####     MRS分段结果作图      #####
####################################

# 输入的ret序列和regime序列长度必须相同
out <- mrs_output(ret1[2:dim(ret1)[1]],mrs_regime)
# 观察某些条件下的情况
out[out$Regime==2 & out$X3>0,]
out2 <- as.data.frame(out[c(c(1),c(5:8))])
out3 <- melt(out[,c(c(1:num),num+1,num+5)],id=c('id','Regime'))  
out3[,2] <- as.character(out3[,2])
out3[,2] <- paste('Regime',out3[,2],sep=' ')
out3[,2] <- ifelse(out3[,3]=='X1',out3[,2],'')
fig <- ggplot(out3)+geom_bar(position='dodge',stat='identity',aes(y=value,x=(id),fill=variable,width=0.5))+geom_text(aes(label=Regime,y=value,vjust=3,x=(id)),colour="black")
fig


####################################
#####    经济周期regime作图    #####
####################################

# 读取经济周期数据
raw_data <- read.csv("~/Code/R/MRSVAR_model/macro_regime.csv")
temp <- (xts(raw_data[,2:3],as.Date(raw_data[,1])))[paste(start,'/',end,sep='')]
# regime21: 1=宽信用;0=紧信用
regime21 <- temp[3:dim(temp)[1],1]
# regime22: 1=宽货币;0=紧货币
regime22 <- temp[3:dim(temp)[1],2]
# regime23: 1=紧货币紧信用;2=紧货币宽信用;3=宽货币紧信用;4=宽货币宽信用
regime23 <- regime21+2*regime22+1
# 共有三种regime类型
multi_regime_plot(regime21, cum_ret)
multi_regime_plot(regime22, cum_ret)
multi_regime_plot(regime23, cum_ret)
                     

####################################
#####    经济周期regime作图    #####
####################################

# 1.指定某一宏观经济状态后计算MRS每个状态的各资产收益
regime_id <- 1
temp1 <- ret1[2:dim(ret1)[1]]
temp2 <- merge(temp1[regime23==regime_id],mrs_regime[regime23==regime_id])  
names(temp2)[length(names(temp2))] <- 'id'
temp_agg <- aggregate(temp2,FUN=mean,by=list(temp2$id))
out3 <- melt(as.data.frame(temp_agg),id=c('id'))  
fig <- ggplot(out3)+geom_bar(position='dodge',stat='identity',aes(y=value,x=(id),fill=variable,width=0.5))
fig


# 2. 指定某一宏观经济阶段后计算MRS每个状态的各资产收益
k <- 2
out4 <- mrs_output(ret1[2:dim(ret1)[1]],regime23)
start_t <- out4[k,7]
end_t <- out4[k,8]

out5 <- ret1[paste(start_t,'/',end_t,sep='')]
out6 <- merge(out5,mrs_regime[paste(start_t,'/',end_t,sep='')])
names(out6)[length(names(out6))] <- 'id'
temp_agg <- aggregate(out6,FUN=mean,by=list(out6$id))
out7 <- melt(as.data.frame(temp_agg),id=c('id'))  

fig <- ggplot(out7)+geom_bar(position='dodge',stat='identity',aes(y=value,x=(id),fill=variable,width=0.5))
fig


    