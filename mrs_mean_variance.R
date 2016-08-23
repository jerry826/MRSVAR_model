source('C:\\Users\\lvjia\\Google �ƶ�Ӳ��\\��һ���\\bd quant research\\�����ʲ�����\\code\\MRS.R')



# Markow state VAR model
w.start()
# ����ʱ�����
start <- '2005-02-01'
end <- '2016-08-16'
mid <- '2013-08-01'
freqs <- 'w'

# �����ʲ���
assets <- c('000300.SH','000905.SH','037.CS')
assets <- c('000300.SH','000905.SH','037.CS','AU9999.SGE')



# ��ȡ���� 
ret_all <- collect_data(assets,start,end,mid,freqs)

# ����ģ�ͱ���
p <- 1 # �ͺ����
h <- 6 # ״̬���� 
# ѵ��ģ��
result <- mrs_model(ret_all$train,p,h)


# rr <- mrs_model(ret$train) # ��������֤

# ģ��Ȩ�����Ż�
perf <- mrs_predict(result,ret_all$test,ret_all$train,h,1)

# ��λ���Ż�
weight2 <- adjust_weight(perf)

# ���ַ���
n <- dim(ret_all$test)[1]
## ȡ2:n�ڵ�����
asset_ret <- ret_all$test[2:n]
## ����
perf_analysis(weight2$weight_adj,asset_ret ,cost=0.001,freq=freqs)

plot.zoo(result$p)

mu = matrix(0,4,h)
for (i in 1:h) {
  aa =t(result$var_beta[,,i])
  mu[,i] = solve(diag(1,4)-aa[,1:4],aa[,5])
}

print(mu)
