# scRefAudit 论文复现材料

- 三份结果工作簿见 [workbooks](workbooks/)。工作表 Table_4 对应新稿表 1，其余编号对照见 [英文说明](README.md)。
- 图 7 的肝髓系、肺上皮及肺间充质数据共 78 份，见 [预处理矩阵](matched_tv_20261118/prepared/)。绘图脚本、实验协议和存档结果在 [records](records/)。
- 图 2B／表 1 的复核排序和图 5 的平衡 OT 计时有[独立固定下载](https://github.com/xiaoweichu873-cmyk/scRefAudit/archive/d728fff33f0ac58840f62e2d571764fef42477c3.zip)，含精确输入、双语说明、运行命令和预期结果。
- 探索性非平衡 OT 的[推导和基准 PDF](unbalanced/Supplementary_Addendum_Unbalanced_OT.pdf)以及[完整代码与输入](https://github.com/xiaoweichu873-cmyk/scRefAudit/tree/7969588a5194bf6b8e02256992139ae92e4c22bc/manuscript_data/unbalanced_ot_20261123)可单独查看。

下载整个快照并解压，在根目录执行 `python3 expand_records.py`，可还原研究档案目录并自动核验文件。随后运行 `python3 figure_sources/regenerate_figures.py` 可从存档值重画图表。需要从预处理矩阵重算组织实验时，按 `matched_tv_20261118/README.md` 的命令运行。上游 GEO 表达矩阵生成 PCA 的步骤仍待验证。
