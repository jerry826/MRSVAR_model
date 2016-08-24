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

# 设置模型变量
p <- 1 # 滞后阶数
h <- 6 # 状态个数 
# 训练模型
result <- mrs_model(ret_all$train,p,h)


# rr <- mrs_model(ret$train) # 样本内验证

# 模型权重最优化
perf <- mrs_predict(result,ret_all$test,ret_all$train,h,1)

# 仓位再优化
weight2 <- adjust_weight(perf)


# 表现分析
## 
n <- dim(ret_all$test)[1]
## 取2:n期的数据
asset_ret <- ret_all$test_d[2:n]
## 分析
perf_analysis(weight2$weight_adj,asset_ret ,cost=0.001,freq=freqs)

plot.zoo(result$p)

mu = matrix(0,4,h)
for (i in 1:h) {
  aa =t(result$var_beta[,,i])
  mu[,i] = solve(diag(1,4)-aa[,1:4],aa[,5])
}

print(mu)


p_all = rbind(result$p, perf$p_list) #把概率全部拼起来

#找到最大概率对应的状态
p_max = apply(p_all, 1, max)
p_max_loc = rep(0, length(p_max))
for (i  in 1:length(p_max)) {
  p_max_loc[i] = which(p_all[i,]==p_max[i])
}

windows()
plot(p_max_loc)

time_loc = c()
for (j in 1:h) {
  t_loc = which(p_max_loc==j)
  time_loc =c(time_loc, list(t_loc))
}
