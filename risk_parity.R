source('C:\\Users\\lvjia\\Google 云端硬盘\\研一暑假\\bd quant research\\大类资产配置\\MRS.R')

w.start()
start <- '2005-02-01'
end <- '2016-08-16'
mid <- '2013-08-01'
freq <- 'd'
assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE',"CUFI.WI")

# assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE','GF0002')
ret <- collect_data(assets,start,end,mid,freq)


h <- weight_cal(covar(ret$test, 60))
perf <- list(weight=h)
# 仓位再优化
weight2 <- adjust_weight(perf)

asset_ret <- ret$test
n <- dim(ret$test)[1]
asset_ret <- ret$test[60:n]
perf_analysis(weight2$weight_adj,asset_ret ,cost=0.001)
