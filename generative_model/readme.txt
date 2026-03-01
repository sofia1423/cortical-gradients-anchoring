①Z:\work\spower_all\figure_7\generative_2026_119
这里面是有了4bin的规则 先去预测8190的矩阵的权重

②Z:\work\spower_all\figure_7\generative_2026_119\asy_40area
有了8190的权重
2.1 下一步这个文件夹里面有一个asymmetry_code_final_65  里面最下面1153行的版本 是用40x40的实验矩阵去计算后验概率
两版：91x40和40x40的分了不同的bin 最后留了91x40分24个bin结果是data_nd_6182_40area2_91_40_6_30.mat里面的第10行
2.2 下一步使用cal_for_different_Previse_40_40weight_91_40P (9191P和9191Pparfor的都是错误的版本)
去随机留边 并计算和afferent的相关
结果是:
afferent_pre_432
matrix_save_40area_40_40weight_91_40P1670  
rangepre_40area_40_40weight_91_40P1670
第432是最后用的 为了保第三个cluster 选了0.4935的版本
以上就是40area的最终预测矩阵 以及预测的afferent range

validation-40area：
validation一共有两个版本：
1. 是40area得到的8190矩阵的权重的validation  舍弃
结果在：Z:\work\spower_all\figure_7\generative_2026_119\validation_40area
但是这个结果有问题，因为validation_4bin里面 我用的24脑区train 所以取边（2274(91x40脑区的边数)有问题 这块不对 不是2274  因为没用这个validation  所以不管了）
如果需要不加后验概率的validation需要和老师确认一下重新算

2. 是40area得到的8190矩阵的权重，然后进行贝叶斯随机留边得到的矩阵的validation 也是最终保留的validation版本
步骤①：Z:\work\spower_all\figure_7\generative_2026_119\validation_40area\validation_asy
使用asymmetry_train_test_code
得到1000个train矩阵（24脑区）的P值
步骤②：Z:\work\spower_all\figure_7\generative_2026_119\validation_40area\validation_asy\cluster
需要的数据和代码全在这 有了P值 然后在集群上算 数据保留在cluster_data文件夹了
这一步得到1000个train的贝叶斯删边之后的矩阵：结果在：load('Z:\work\spower_all\figure_7\generative_2026_119\validation_40area\validation_asy\cluster\results_992\matrix_1_1000.mat')
步骤③：返回Z:\work\spower_all\figure_7\generative_2026_119\validation_40area\validation_asy
代码里面 计算了sen 和roc(auc的值)
但是因为matrix_f 保留的边只有5~7k 所以roc画不全，会影响auc的值
所以需要继续在matrix_f的基础上 放边
代码在Z:\work\spower_all\figure_7\generative_2026_119\validation_40area\validation_asy\cluster\results_992
code_matrixf_to_8190
但是其实不用每个matrix_f都放到8190
因为sen分布已经有了 我只需要把auc最高的 接着放边（这一版是第808个trian的矩阵）

同样也可以把40脑区的接着放边 看看哪个效果好：也就是把上面2.2最终预测的数据 接着放边看一下ROC 和auc、sen
要么就是40脑区的  要么选20脑区train最高的

最后的结果：
①train最高的Z:\work\spower_all\figure_7\generative_2026_119\validation_40area\validation_asy
train_sen
test_sen
roc_fig_train_808_5270  这个是train最高的roc
②40脑区接着放边的




2026/2/7测试新的:40x40去预测权重 91x91(newp402)去算贝叶斯
首先预测权重规则pre_weight_325_40_40 通过40x40 的999条实验边
然后预测8190权重： pre_weight_2_40_40
