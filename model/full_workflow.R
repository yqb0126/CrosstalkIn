#----------------------------------------------------
# 第1步: 加载必要的库
#----------------------------------------------------
library(igraph)
library(parallel)

#----------------------------------------------------
# 第2步: 设置文件路径
#----------------------------------------------------
# STRING数据库文件路径
links_file <- "./model/9606.protein.links.v12.0.txt"  # 蛋白质互作关系
protein_info_file <- "./model/9606.protein.info.v12.0.txt"  # 蛋白质信息
# 输出文件
ppi_network_file <- "./model/ppi_network.txt"  # 处理后的PPI网络
ppi_weights_file <- "./model/ppi_weight.txt"  # 计算后的基因权重

#----------------------------------------------------
# 第3步: 处理STRING数据库文件
#----------------------------------------------------
# 首先将蛋白质ID转换为基因符号并提取PPI网络
source("./model/process-string-ppi.r")
ppi_network <- process_string_ppi_dt(
  links_file = links_file,
  protein_info_file = protein_info_file,
  output_file = ppi_network_file,
  score_threshold = 700  # 使用700分作为阈值，可以根据需要调整
)

# 显示网络基本信息
message("PPI网络基本信息:")
message("- 节点数量: ", length(unique(c(ppi_network$gene1, ppi_network$gene2))))
message("- 边数量: ", nrow(ppi_network))
message("- 平均权重: ", mean(ppi_network$weight))

#----------------------------------------------------
# 第4步: 获取基因表达数据
#----------------------------------------------------
GEP <- readRDS("./data/sdy311.rds")
GEP <- GEP$expr

# 检查基因交集
genes_in_expr <- rownames(GEP)
genes_in_ppi <- unique(c(ppi_network$gene1, ppi_network$gene2))
common_genes <- intersect(genes_in_expr, genes_in_ppi)

message("基因交集信息:")
message("- 表达数据中的基因数: ", length(genes_in_expr))
message("- PPI网络中的基因数: ", length(genes_in_ppi))
message("- 共有基因数: ", length(common_genes))
message("- 交集比例: ", round(length(common_genes) / length(genes_in_expr) * 100, 2), "%")

#----------------------------------------------------
# 第5步: 计算PPI网络权重
#----------------------------------------------------
# 计算基因在PPI网络中的权重
source("./model/calculate-ppi-weights.r")
gene_weights <- calculate_ppi_weights(
  PPI_network = ppi_network,
  gene_list = genes_in_expr,
  output_file = ppi_weights_file,
  damping_factor = 0.85  # PageRank默认重启概率
)

#----------------------------------------------------
# 第6步: 运行CrosstalkIn_PPI
#----------------------------------------------------
# 确定要使用的CPU核心数
cores_to_use <- min(8, parallel::detectCores() - 1)
message("使用 ", cores_to_use, " 个CPU核心进行分析")
message("开始分析...")

# 读取 cellmarkers.csv 文件
cellmarkers_df <- read.csv("./model/cellmarkers86.csv", row.names = 1, header = FALSE)
# 初始化一个空列表用于存储细胞标记基因
cell_markers <- list()
# 遍历数据框的每一行
for (cell_type in rownames(cellmarkers_df)) {
  # 获取当前细胞类型对应的标记基因字符串
  markers_str <- cellmarkers_df[cell_type, 1]
  # 使用逗号分隔标记基因字符串，转换为向量
  markers <- strsplit(markers_str, ",")[[1]]
  # 将标记基因向量添加到 cell_markers 列表中
  cell_markers[[cell_type]] <- markers
}
# 计算网络融合得分
message("开始计算网络融合得分...")

source("./model/network_fusion.R")
fused_scores <- fuse_networks(GEP = GEP,
                              cell_markers = cell_markers,
                              sample_names = colnames(GEP),
                              marker_weight = 0.5,
                              networks_dir = file.path(getwd(), "networks"),
                              # networks_dir = file.path("D:/postgra_code/CITMIC_3/sdy311_newnew/networks78/"),
                              ppi_weights_file = file.path("./model/ppi_weight.txt"),
                              inter_path = "./model/surv_inter.csv",
                              jacaard_path = "./model/surv_jaccard.csv",
                              cl.cores = cores_to_use)

#保存融合得分结果
# message("保存融合得分结果...")
saveRDS(fused_scores, "score_sdy311_0.5_923.rds")

