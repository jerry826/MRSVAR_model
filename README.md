# MRSVAR_model
Markov regime swtich and VAR model

## 基本函数
func.R:提供了基本的数据调取、整合、回测函数和模型训练函数等 \
mrs_mean_variance.R: 主要策略模型，基于MRS VAR模型进行资产配置输出结果\
multi_strategy_analysis.R:对比了几类常用的大类资产配置模型，例如简单趋势追踪、等权再平均法等\
msb_var.R :采用Gibbs抽样和贝叶斯回归方法计算MRS模型，结果不稳定不建议使用\ 
