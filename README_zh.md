# scRefAudit：找出该优先复核的单细胞标签

同一个细胞，表达没变，只调整参考库里各类细胞的权重，注释标签也可能变。
scRefAudit 是一道平衡最优传输（OT）注释质控：在用标签做下游比较之前，
先检查哪些标签依赖参考比例。

- **排复核顺序。** `review_cells()` 按估计的翻标签所需比例变化，从小到大给细胞排序。
- **检查指定细胞。** `stress_test()` 针对你关心的细胞提出比例变化，再求解确认标签会不会变。
- **比较指定参考比例。** `audit_compositions()` 输出哪些细胞换了标签、各类得分变了多少。

在两套免疫细胞数据的 8 个配置中，用同一组独立抽样的参考比例来检验，
每复核 100 个细胞，scRefAudit 平均找到 **70.63 个对参考比例敏感的细胞**，
只按前两名分差排序则找到 **58.25 个**。标签会变，说明值得复核，不等于原标签错了。

输入是固定的「查询细胞 × 参考细胞」距离矩阵与参考细胞标签。R 和 Python 共用同一个数值核心。
[在 RStudio 试跑](#在-rstudio-安装) · [English README](README.md)

## 这一版怎么用

这一版有两个入口：

- **Python 包**：`screfaudit`。
- **R 包接口**：`scRefAudit`，通过 reticulate 调用同一份 Python 核心，可以在 RStudio 中使用。

论文原有求解代码逐字节保留。新加的是安装结构、输入检查、用户接口、说明文档、示例和测试。
软件命名没有改动论文的标题、预印本 DOI 或数值结果。

## 在 RStudio 安装

先[下载 GitHub 源码 ZIP](https://github.com/xiaoweichu873-cmyk/scRefAudit/archive/refs/heads/main.zip)并解压。
把下面的路径换成实际解压出来的 `scRefAudit-main` 文件夹：

```r
install.packages("reticulate")
install.packages("/你的路径/scRefAudit-main", repos = NULL, type = "source")
library(scRefAudit)

x <- demo_audit_data()
fit <- audit_fit(x$cost, x$labels, epsilon = x$epsilon)

result <- stress_test(fit, budget = 16L, radius = 0.01,
                      metric = "tv", seed = 37L)
result$comparison
review_cells(fit, top_n = 10L)
```

路径应指向解压后的文件夹，里面应有 `DESCRIPTION` 文件。
首次真正计算时，reticulate 可能需要联网准备 Python、NumPy 和 SciPy。
R 要求至少 4.1，reticulate 至少 1.41。Apple Silicon 应使用同架构的 R 和 Python。
这一版尚未上架 CRAN、Bioconductor 或 PyPI；本机测试情况见 `validation/VALIDATION.md`。

## 接实际数据时需要什么

输入是已经准备好的“查询细胞 × 参考细胞”距离矩阵，以及参考细胞标签。
R 行名可以用真实细胞条形码；返回结果的 `cell_id` 用于对应回原细胞。
它不是从原始 count 矩阵开始的一键单细胞流程，也不会自动加速 Seurat 的其他步骤。

一次模型拟合后，保持细胞、距离和正则化参数不变，可以反复调整参考类别比例。
如果更换参考细胞或重做 PCA，需要重新拟合。

## 几个参数的人话解释

- `radius = 0.01, metric = "tv"`：实际参考组成移动 1% 的总质量。
- `budget = 16L`：定向与随机各检查 16 个场景。
- `order = 0L`：共享曲率，沿用基准解作起点。
- `order = 1L` / `2L`：加入一阶 / 二阶预测，再校正到设定精度。
- `review_cells()`：按局部线性化半径从小到大，列出优先复核的细胞。

标签会变表示依赖参考比例，并不直接说明原来的标签错了。
所有 R 细胞序号从 1 开始；Python 序号从 0 开始。

## 加速与适用范围

反复检查不同比例时，求解器复用基准矩阵分解，也可用泰勒预测作起点。
论文的单线程 CPU 基准实验中，研究代码在与 Anderson 加速 Sinkhorn（AA2）
都完成的 11 个配对配置中有 10 个更快，中位加速 **2.55 倍**。
这是指定数值工作流的结果，不是 R/Python 接口在任意数据上的速度承诺。

## 非平衡 OT 探索性扩展（2026-11-23）

新增了[非平衡 OT 补充说明、结果和完整复现材料](manuscript_data/unbalanced_ot_20261123/README_zh.md)。
这次在双边 KL 边际惩罚下测试共享曲率与一阶预测：96 种配置中 57 种更快，
审计耗时中位加速 **1.34 倍**，计入共同初始化为 **1.26 倍**。
局部扰动中位 **2.26 倍**，较大扰动中位 **0.98 倍**；变慢的配置也全部保留。

首轮五对照比较与回退修订后的复测使用同一组开发数据，不能当作两次独立验证。
这部分作为独立计算原型提供，当前 R/Python 包接口仍针对平衡 OT。

## 放到 GitHub

网站：[https://github.com/](https://github.com/)

1. 登录后新建仓库，建议仓库名 `scRefAudit`。
2. 解压 GitHub 源码 ZIP，进入里面的 `scRefAudit` 文件夹。
3. 上传这个文件夹里的文件和子文件夹，让 `README.md`、`DESCRIPTION` 和 `pyproject.toml` 位于仓库根目录。
4. 不需要把论文正文、投稿信、作者内部清单一起上传。
5. 发布固定版本时可用标签 `v0.1.0`，再通过 GigaDB 或 Zenodo 保存论文对应版本。

网页拖拽通常不会包含隐藏的 `.github` 目录。可用 GitHub Desktop 上传整个目录以保留自动测试配置，
或在网页中另行添加 `.github/workflows/check.yml`。上传后检查 Actions 的 R 与 Python 测试结果。

代码公开前确认两位作者的软件署名、MIT 许可和引用信息。
本仓库提供已发布的 scRefAudit 源码及论文数据下载入口。

## 这份代码包与论文附件的关系

[论文复现材料固定版本](https://github.com/xiaoweichu873-cmyk/scRefAudit/tree/859f41bc4c9752bb02580d557bf2afe2d79b8e58)包含三份工作簿、图 7 的 78 份组织分析输入、绘图与核验代码，以及探索性非平衡 OT 补充说明。[点击下载整个档案](https://github.com/xiaoweichu873-cmyk/scRefAudit/archive/859f41bc4c9752bb02580d557bf2afe2d79b8e58.zip)。

本包有可运行示例、数值核心和三个功能的接口。它不包含论文全部大型数据、计时对照和第三方求解器构建材料。
图 2B、表 1 和图 5 的复核材料见[固定平衡 OT 下载](https://github.com/xiaoweichu873-cmyk/scRefAudit/archive/d728fff33f0ac58840f62e2d571764fef42477c3.zip)。
