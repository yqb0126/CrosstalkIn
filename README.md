# CrosstalkIn

流程：从 STRING PPI 构建基因互作网络、用 PageRank 计算基因权重、构建细胞-细胞网络并按权重融合，最终为每个样本输出融合得分矩阵。

## 功能概览

- STRING PPI 处理：蛋白 ID 映射到基因符号并生成基因-基因加权网络
- PPI 权重计算：在 PPI 图上用 PageRank 计算基因权重
- 网络构建与融合：基于 marker 网络与/或表达衍生网络构建细胞-细胞网络并融合
- 样本得分：在融合网络上用 PageRank 得到每个样本的细胞得分并归一化

## 安装

本项目既可以作为 R 包使用，也可以直接运行工作流脚本。

### 作为 R 包安装

从源码目录安装：

```r
install.packages("data.table")
install.packages("igraph")
install.packages("Matrix")
install.packages("foreach")
install.packages("doParallel")

install.packages("progress")

install.packages("devtools")
devtools::install_local("d:/postgra_code/CrosstalkIn")
```

## 目录结构

- `R/`：R 包的核心实现（可 `library(CrosstalkIn)` 调用）
- `man/`：导出的函数文档
- `model/`：工作流脚本与示例数据（例如 `process-string-ppi.r`、`network_fusion.R`）
- `run_model.Rmd`：示例工作流（推荐从这里快速跑通）

## 快速开始

### 方式 A：直接跑 R Markdown 工作流

打开并 Knit [run_model.Rmd](file:///d:/postgra_code/CrosstalkIn/run_model.Rmd)，修改 YAML 中的 `params` 指向你的输入文件即可。

该工作流默认包含 5 步：

1. 处理 STRING 数据库生成 PPI 网络
2. 加载基因表达矩阵（GEP）并统计与 PPI 的交集
3. 基于 PPI 计算基因权重（PageRank）
4. 读取细胞 marker 并进行网络融合得到每样本得分
5. 保存结果与 `sessionInfo()`

### 方式 B：用 R 包一键跑完整流程

核心入口函数是 [our_model_run](file:///d:/postgra_code/CrosstalkIn/R/our_model_run.R)：

```r
library(CrosstalkIn)

scores <- our_model_run(
  links_file = "data/9606.protein.links.v12.0.txt",
  protein_info_file = "data/9606.protein.info.v12.0.txt",
  expr_rds = "data/sdy311.rds",
  cellmarkers_csv = "data/cellmarkers86.csv",
  shortest_paths_csv = "data/shortest_paths.csv",
  inter_path = "data/surv_inter.csv",
  jacaard_path = "data/surv_jaccard.csv",
  marker_weight = 0.5,
  score_threshold = 700,
  damping_factor = 0.85,
  damping = 0.90,
  networks_dir = "networks",
  cl.cores = 8,
  output_dir = "result",
  save_scores_to = "result/scores.rds"
)
```

## 输入文件说明

### 1) STRING 相关

- `links_file`：STRING `protein.links` 文件，至少包含列 `protein1`、`protein2`、`combined_score`
- `protein_info_file`：STRING `protein.info` 文件，至少包含 `preferred_name`，以及 ID 列之一：
  - `protein_external_id` / `X.string_protein_id` / `#string_protein_id` / `string_protein_id`

处理逻辑见 [process_string_ppi_dt](file:///d:/postgra_code/CrosstalkIn/R/process_string_ppi_dt.R)。

### 2) 表达矩阵（GEP）

- `expr_rds`：RDS 文件，内容可以是
  - `matrix/data.frame`（基因×样本），或
  - `list(expr = <matrix>)`
- 必须满足：
  - `rownames(GEP)` 是基因名（建议与 STRING 映射到的基因符号一致）
  - `colnames(GEP)` 是样本名

### 3) 细胞 marker

- `cellmarkers_csv`：CSV 文件，行名为细胞类型；第 1 列为 marker 基因，逗号分隔
  - 仓库示例：[model/cellmarkers86.csv](file:///d:/postgra_code/CrosstalkIn/model/cellmarkers86.csv)

### 4) 网络融合相关

`marker_weight` 控制 marker 网络与表达网络的融合比例：

- `marker_weight = 1`：只用 marker 网络
- `marker_weight = 0`：只用表达衍生网络
- 介于 0 与 1：两者加权融合

对应地，以下输入在不同设置下会被用到：

- `shortest_paths_csv`：当 `marker_weight > 0` 时需要，用于 marker 网络构建（基因-基因最短路距离矩阵）
- `inter_path` 与 `jacaard_path`：当 `marker_weight < 1` 时需要，用于表达衍生网络构建
- `ppi_weights_file`：PPI 权重文件（通常由流程自动生成）

实现见 [fuse_networks](file:///d:/postgra_code/CrosstalkIn/R/fuse_networks.R)。

## 输出

- `our_model_run()` 返回：`matrix`（行=细胞类型，列=样本），数值为归一化后的融合得分
- `output_dir` 下的中间产物（示例）：
  - `ppi_network.txt`：基因-基因加权互作表
  - `ppi_weight.txt`：基因权重（PageRank）
- `networks_dir` 下会缓存每个样本的中间网络，便于断点续跑：
  - `marker_<sample>.rds`
  - `expression_<sample>.rds`

## 进度条与性能提示

如果你看到类似：

```
format = Progress [--------------------]   1% (elapsed:  4m)
```

这来自 [fuse_networks](file:///d:/postgra_code/CrosstalkIn/R/fuse_networks.R) 内部对 `sample_names` 的循环进度条（依赖 `progress` 包）。`1%` 表示“按样本处理”的循环只完成了约 1%。

加速建议：

- 合理设置 `cl.cores`（仅影响表达衍生网络的并行分发部分）
- 使用固定的 `networks_dir`，确保缓存命中，避免重复计算

## License

GPL-3（见 [DESCRIPTION](file:///d:/postgra_code/CrosstalkIn/DESCRIPTION)）。
