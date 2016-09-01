# author: Min Song, Jiayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: mrs_mean_variance.R
# time: 2016/9/1

source('func.R')

# Markow state VAR model

w.start()

# 1. 基本参数设置

## 1.1 设置时间变量
start <- '2005-02-01' # 训练集第一期
end <- '2016-08-16'   # 测试集最后一期
mid <- '2013-08-01'   # 训练集最后一期
freq <- 'w' # 模型频率

## 1.2 设置资产池
assets <- c('000905.SH','037.CS')
assets <- c('000300.SH','000905.SH','037.CS')

## 提取数据 
ret_all <- collect_data(assets,start,end,mid,freqs)

## 1.3 设置模型变量
p <- 1 # 滞后阶数
h <- 5 # 状态个数 

# 2. 训练模型输出MRS-VAR模型结果
result <- mrs_model((ret_all$train),p,h,niterblkopt=20)


# 3. 模型资产配置最优化
weight <- mrs_predict(result,ret_all$test,ret_all$train,h,1)

# 4. 仓位再优化，减少调仓
weight_adj <- adjust_weight(weight)


# 5. 表现分析
## 
n <- dim(ret_all$test)[1]
## 取2:n期的数据
asset_ret <- ret_all$test_d[2:n]
## 分析
model_ <- perf_analysis(weight_adj$weight_adj,asset_ret,cost=0.001,freq=freq)

plot.zoo(result$p)


# 6. 其他结果输出


## 6.1 计算各状态均值
mu = matrix(0,length(assets),h)
for (i in 1:h) {
  aa =t(result$var_beta[,,i])
  mu[,i] = solve(diag(1,length(assets))-aa[,1:length(assets)],aa[,length(assets)+1])
}
print('各个状态下资产的均值')
print(mu)


## 6.2 计算每天各个状态的概率

p_all = rbind(result$p, perf$p_list) #把概率全部拼起来
plot.zoo(p_all)
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



