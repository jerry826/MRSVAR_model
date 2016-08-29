source('MRS.R')
# Markow state VAR model
w.start()
# 设置时间变量
start <- '2005-02-01'
end <- '2016-08-16'
mid <- '2013-08-01'
freqs <- 'w'

# 设置资产池
assets <- c('000300.SH','000905.SH','037.CS')
assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE')

# 提取数据 
ret_all <- collect_data(assets,start,end,mid,freqs)


set.seed(123)

mu <- apply(ret_all$train,2,mean)
st <- apply(ret_all$train,2,sd)
for (i in 1:dim(ret_all$train)[1]){
  ret_all$train[i,] <- (ret_all$train[i,]-mu)/st
} 
for (i in 1:dim(ret_all$test)[1]){
  ret_all$test[i,] <- (ret_all$test[i,]-mu)/st
} 




names(dd) <- c('x1','x2','x3','x4')
h=4

xm <- msbvar(ts(ret_all$train), p=1, h=h,
             lambda0=0.8, lambda1=0.1,
             lambda3=2, lambda4=0.1, lambda5=0.01, mu5=5,
             mu6=5, qm=12,max.iter =15)

# Plot out the initial mode
plot(ts(xm$fp))
print(xm$Q)

# Now sample the posterior
N1 <- 4000
N2 <- 4000

# First, so this with random permutation sampling
# x1 <- gibbs.msbvar(xm, N1=N1, N2=N2, permute=TRUE)

# Identify the regimes using clustering in plotregimeid()

# plotregimeid(x2, type="all")

# Now re-estimate based on desired regime identification seen in the
# plots. Here we are using the intercept of the first equation, so

x2 <- gibbs.msbvar(xm, N1=N1, N2=N2, permute=FALSE, Beta.idx=c(3,2))
# Plot regimes
plot.SS(x2)
# Summary of transition matrix
summary(x2$Q.sample)
# Plot of the variance elements
plot(x2$Sigma.sample)
## End(Not run)

h <- 4
m <- 4

Q <- matrix(apply(x2$Q.sample,2,mean),h,h)
S <- matrix(apply(x2$Sigma.sample,2,mean),m*(m+1)/2,h)
beta <- matrix(apply(x2$Beta.sample,2,mean),m*(m+1),h)
# TT <- matrix(apply(x2$transition.sample,2,mean),h,h)

plotregimeid(x2, type="all")


beta1 <- array(as.vector(beta),c(1+m,m,h))
Q1 <- Q
S1 <- array(0,c(m,m,h))
for (i in (1:h)){
  ii <- 1
  for (j in (4:1)){
    for (k in (1:j)){
      S1[5-j,k+4-j,i] <- S[ii,i]
      S1[k+4-j,5-j,i] <- S[ii,i]
      
      ii <- ii +1
    }
  }
} 


p_list <- c()
TT <- dim(ret_all$test)[1]

for (t in 1:TT){
  if (t==1){
    r_yst <- as.numeric(ret_all$train[dim(ret_all$train)[1],])
    p0 <- result$p[100,]
    
  }else{
    r_yst <- as.numeric(ret_all$test[t-1,])
    p0 <- p_update
  }
  
  r_today <-as.numeric(ret_all$test[t,])
  
  p_next <- p0%*%Q1
  print(r_today)
  p_update<-gx(r_today,r_yst,p,beta1,S1,h)
  p_list <- rbind(p_list,p_update)
  
}
plot.zoo(ts(p_list))












model <- msvar(dd,p=1,h=5,niterblkopt=10)


init <- initialize.msbvar(ts(dd),p=1,h=5,lambda0=1, lambda1=1,
                  lambda3=0, lambda4=1.2, lambda5=0.1, mu5=0,
                  mu6=0,prior=0, qm=12,Q=model$Q)


szbvar



