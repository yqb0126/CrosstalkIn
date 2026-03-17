#' @title process_string_ppi_dt
#' @description 使用data.table优化的STRING数据库PPI网络处理函数，显著提高大数据处理性能
#' @param links_file STRING数据库的蛋白质互作关系文件路径(例如"9606.protein.links.v12.0.txt.gz")
#' @param protein_info_file STRING数据库的蛋白质信息文件路径(例如"9606.protein.info.v12.0.txt")
#' @param output_file 输出的PPI网络文件路径，包含基因符号(默认: "ppi_network.txt")
#' @param score_threshold 蛋白质互作关系的分数阈值(默认: 700，分数范围0-1000)
#' @return 包含gene1、gene2和weight列的数据框，同时将结果保存到指定的输出文件
#' @importFrom data.table fread
#' @importFrom data.table setDT
#' @importFrom data.table setnames
#' @importFrom data.table setkey
#' @importFrom utils read.delim
#' @importFrom utils write.table
#' @export
#' @examples
#' # 处理STRING数据库文件
#' ppi_network <- process_string_ppi_dt(
#'   links_file = "9606.protein.links.v12.0.txt.gz",
#'   protein_info_file = "9606.protein.info.v12.0.txt",
#'   output_file = "ppi_network.txt",
#'   score_threshold = 700
#' )

process_string_ppi_dt <- function(links_file, protein_info_file, output_file = "ppi_network.txt", score_threshold = 700) {
  # 确保加载data.table包
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("请安装data.table包: install.packages('data.table')")
  }
  # 导入命名空间，不需要显式加载包
  library(data.table)
  
  # 输入验证
  if (!file.exists(links_file)) {
    stop("蛋白质互作关系文件不存在: ", links_file)
  }
  
  if (!file.exists(protein_info_file)) {
    stop("蛋白质信息文件不存在: ", protein_info_file)
  }
  #links_file='9606.protein.links.v12.0.txt.gz'
  #protein_info_file='9606.protein.info.v12.0.txt'
  #score_threshold = 900
  output_file="ppi_network.txt"
  # 显示进度信息
  message("正在读取蛋白质互作关系文件...")
  
  # 使用fread高效读取大型文件
  links_data <- fread(links_file, showProgress = TRUE)
  
  # 过滤低分互作关系
  message("过滤分数低于 ", score_threshold, " 的互作关系...")
  links_data <- links_data[combined_score >= score_threshold]
  message("保留了 ", nrow(links_data), " 条蛋白质互作关系")
  
  # 读取蛋白质信息文件，获取蛋白质ID到基因符号的映射
  message("正在读取蛋白质信息文件...")
  protein_info <- fread(protein_info_file)
  
  # 检查必需列是否存在
  if (!all(c("protein_external_id", "preferred_name") %in% colnames(protein_info))) {
    # 尝试其他可能的列名
    if ("X.string_protein_id" %in% colnames(protein_info) && "preferred_name" %in% colnames(protein_info)) {
      # 使用替代列名
      protein_mapping <- protein_info[, .(protein_id = `X.string_protein_id`, preferred_name)]
    } else if ("#string_protein_id" %in% colnames(protein_info) && "preferred_name" %in% colnames(protein_info)) {
      # 另一种可能的列名
      protein_mapping <- protein_info[, .(protein_id = `#string_protein_id`, preferred_name)]
    } else {
      stop("蛋白质信息文件缺少必需的列: protein_external_id 和 preferred_name")
    }
  } else {
    # 创建蛋白质ID到基因符号的映射
    protein_mapping <- protein_info[, .(protein_id = protein_external_id, preferred_name)]
  }
  
  message("读取到 ", nrow(protein_mapping), " 个蛋白质的基因符号信息")
  
  # 创建映射表，提高连接性能
  setkey(protein_mapping, protein_id)
  
  # 提取物种ID前缀 (如果需要)
  # 通常不需要，因为我们直接在连接时处理这个问题
  
  # 使用data.table的连接操作进行快速映射转换
  message("将蛋白质ID转换为基因符号...")
  
  # 首先处理protein1映射
  links_data[, gene1 := protein_mapping[protein1, preferred_name]]
  
  # 然后处理protein2映射
  links_data[, gene2 := protein_mapping[protein2, preferred_name]]
  
  # 检查并移除未能映射的蛋白质对
  na_count <- sum(is.na(links_data$gene1) | is.na(links_data$gene2))
  if (na_count > 0) {
    message("警告: ", na_count, " 条互作关系中的蛋白质ID无法映射到基因符号")
    links_data <- links_data[!is.na(gene1) & !is.na(gene2)]
    message("保留了 ", nrow(links_data), " 条有效的互作关系")
  }
  
  # 计算权重 (将combined_score转换为0-1范围)
  links_data[, weight := combined_score / 1000]
  
  # 处理重复的基因对，确保基因对的统一顺序
  message("处理基因对顺序和重复关系...")
  
  # 创建临时列来标记需要交换的行
  links_data[, need_swap := gene1 > gene2]
  
  # 对需要交换的行进行处理
  swap_rows <- links_data[need_swap == TRUE]
  if (nrow(swap_rows) > 0) {
    # 交换gene1和gene2
    links_data[need_swap == TRUE, `:=`(
      gene1 = gene2,
      gene2 = gene1
    )]
  }
  
  # 删除临时标记列
  links_data[, need_swap := NULL]
  
  # 保留必要的列
  links_data <- links_data[, .(gene1, gene2, weight)]
  
  # 使用data.table的优化聚合函数，这比base R的aggregate快很多倍
  message("对重复的基因对取最高权重值...")
  ppi_network <- links_data[, .(weight = max(weight)), by = .(gene1, gene2)]
  
  message("最终网络包含 ", nrow(ppi_network), " 条基因互作关系")
  
  # 保存结果
  message("将结果保存到文件: ", output_file)
  fwrite(ppi_network, file = output_file, sep = "\t", 
         row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # 返回结果
  return(as.data.frame(ppi_network))
}
