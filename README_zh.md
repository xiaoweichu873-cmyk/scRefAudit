# scRefAudit 平衡最优传输复现

本包对应 GigaScience Technical Note 的图 2B／表 1 和图 5。两个实验都能从指定的预处理嵌入启动。代码版本、上游来源、依赖、测试平台与命令均见 [README.md](README.md)；每个文件的校验值见 `SHA256SUMS`。

| 用途 | 命令 | 预期输出 |
|---|---|---|
| 从公开记录重画图、核算表 | `python scripts/public_tables.py --out output/tables` | 图 5A 的 11 个完整配对中有 10 个更快，中位 2.55 倍。图 5B 中二阶预测相对只复用曲率为 8/12、中位 1.10 倍。 |
| 用已保存的逐细胞换标结果重新排序 | `python scripts/review_cells.py --mode frozen --out output/frozen_review.csv` | 8 个配置的前 100 个细胞中，导数法与分差法平均分别找到 70.63 和 58.25 个敏感细胞。 |
| 从四份预处理输入重新求解 | `python scripts/review_cells.py --mode solve --compare-archive --out output/fresh_review.csv` | 8 个配置，每配置 1,000 个细胞、96 种随机参考组成；逐项核对 768,000 个场景标签。 |
| 从预处理输入重新计时 | `python scripts/timing_trial.py --out output/timing_trial` | 两个配置、每个三个方法、每方法三次独立进程，共 18 条记录。历史试跑中 Kang 的 AA2/提出方法耗时比为 0.86，Hao 为 4.78。新机器的时间会变化。 |

完整计时命令为 `python scripts/timing_trial.py --grid --out output/full_timing`。历史调度共 180 条，152 条完成、10 条超时、18 条因前次失败而跳过。2.55 倍对应论文的单线程研究基准。公开 R/Python 接口没有在这次比较中计时。对应的旧编号为图 7B／表 4（换标排序）和图 4（计时），原工作簿保留了这些旧名称。本次复核起点是准备好的嵌入，从 GEO 原始矩阵到 PCA 的步骤尚未验证。
