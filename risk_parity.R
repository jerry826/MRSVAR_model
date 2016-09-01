# author: Min Song, Jiayang Lv
# contact: lvjy3.15@sem.tsinghua.edu.cn
# file: risk_parity.R
# time: 2016/9/1
source('func.R')

# 1. 基本参数设置

## 1.1 设置时间变量
w.start()
start <- '2005-02-01' # 训练集第一期
end <- '2016-08-16'   # 测试集最后一期
mid <- '2013-08-01'   # 训练集最后一期
freq <- 'd' # 模型频率
assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE',"CUFI.WI")

# 2. 获取数据
ret <- collect_data(assets,start,end,mid,freq)
# 3. 计算历史协方差矩阵
weight <- weight_cal(covar(ret$test, 60))

weight <- list(weight=weight)
# 仓位再优化
weight_adj <- adjust_weight(weight)

asset_ret <- ret$test
n <- dim(ret$test)[1]
asset_ret <- ret$test[60:n]
ana <- perf_analysis(weight2$weight_adj,asset_ret ,cost=0.001)

