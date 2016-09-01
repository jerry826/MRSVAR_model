source('func.R')

# 1. 基本参数设置

## 1.1 设置时间变量
w.start()
start <- '2005-02-01' # 训练集第一期
end <- '2016-08-16'   # 测试集最后一期
mid <- '2013-08-01'   # 训练集最后一期
freq <- 'd' # 模型频率
assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE',"CUFI.WI")

# 
ret <- collect_data(assets,start,end,mid,freq)

h <- weight_cal(covar(ret$test, 60))

perf <- list(weight=h)
# 仓位再优化
weight2 <- adjust_weight(perf)

asset_ret <- ret$test
n <- dim(ret$test)[1]
asset_ret <- ret$test[60:n]
ana <- perf_analysis(weight2$weight_adj,asset_ret ,cost=0.001)

