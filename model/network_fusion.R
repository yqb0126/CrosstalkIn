# 必要的库
library(data.table)
library(Matrix)
library(igraph)
library(parallel)
library(doParallel)
library(progress)


fuse_networks <- function(expr_networks = NULL, marker_networks = NULL,
                          GEP = NULL, cell_markers = NULL,
                          sample_names = NULL, marker_weight = 1,
                          weighted = TRUE, damping = 0.90,
                          networks_dir = "./sdy420/sdy311_networks",
                          ppi_weights_file = "ppi_weight.txt",
                          inter_path = "surv_inter.csv",
                          jacaard_path = "surv_jaccard.csv",
                          cl.cores = parallel::detectCores() - 1) {
  message("开始计算！")

  # 参数检查
  if (is.null(sample_names)) {
    if (!is.null(GEP)) {
      sample_names <- colnames(GEP)
      message(sprintf("使用GEP矩阵的列名作为sample_names: %s",
                      paste(sample_names[1:min(3, length(sample_names))], collapse=", "),
                      if(length(sample_names) > 3) "..." else ""))
    } else {
      stop("必须提供sample_names参数或包含样本名的GEP矩阵")
    }
  }

  if (length(sample_names) == 0) {
    stop("sample_names不能为空")
  }

  # 创建网络存储目录
  if (!dir.exists(networks_dir)) {
    dir.create(networks_dir, recursive = TRUE)
  }
  # 如果没有提供网络，则尝试从文件加载
  if (is.null(sample_names)) {
    stop("如果不提供网络，必须提供sample_names参数")
  }

  # 设置并行计算
  cl <- makeCluster(cl.cores)
  on.exit(stopCluster(cl))  # 确保函数退出时关闭集群
  registerDoParallel(cl)

  # 创建进度条
  pb <- progress_bar$new(
    total = length(sample_names),
    format = "计算网络进度 [:bar] :percent 已完成 (用时: :elapsed)",
    clear = FALSE,
    width = 60
  )

  # 并行计算marker网络
  if (marker_weight > 0) {
    marker_networks <- list()
    library(readr)
    # 读取并处理最短路径矩阵
    shortest_path_matrix <- read_csv("./model/shortest_paths_optimized.csv")
    names_vec <- shortest_path_matrix[[1]]
    shortest_path_matrix <- as.matrix(shortest_path_matrix[, -1])  # 直接转换为矩阵
    rownames(shortest_path_matrix) <- colnames(shortest_path_matrix) <- names_vec

    # 确保最短路径矩阵是数值型
    shortest_path_matrix <- matrix(as.numeric(shortest_path_matrix),
                                   nrow = nrow(shortest_path_matrix),
                                   dimnames = list(rownames(shortest_path_matrix),
                                                   colnames(shortest_path_matrix)))

    # 读取PPI权重
    ppi_data <- fread(ppi_weights_file, sep = "\t")
    ppi_weights <- setNames(ppi_data$weight, ppi_data$gene)

    # 预处理权重矩阵
    all_genes <- unique(c(names(ppi_weights), rownames(shortest_path_matrix)))
    weight_matrix <- matrix(1, nrow = length(all_genes), ncol = length(all_genes))
    rownames(weight_matrix) <- colnames(weight_matrix) <- all_genes

    # 使用向量化操作设置权重
    weight_indices <- match(names(ppi_weights), all_genes)
    weight_matrix[weight_indices, ] <- ppi_weights
    weight_matrix[, weight_indices] <- ppi_weights

    epsilon <- 1e-6
    cell_types <- names(cell_markers)

    for (sample in sample_names) {
      marker_file <- file.path(networks_dir, paste0("marker_", sample, ".rds"))
      message(sprintf("处理PPI网络: %s", sample))

      if (file.exists(marker_file)) {
        marker_networks[[sample]] <- readRDS(marker_file)
        next
      }
      sample_expr <- GEP[, sample, drop = FALSE]
      if (ncol(sample_expr) == 0 || nrow(sample_expr) == 0) {
        stop(sprintf("样本 %s 的表达数据为空", sample))
      }

      # 初始化结果矩阵
      path_score_mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types))
      rownames(path_score_mat) <- colnames(path_score_mat) <- cell_types
      # 预计算所有基因的表达值
      expr_values <- as.numeric(sample_expr[rownames(sample_expr), ])
      names(expr_values) <- rownames(sample_expr)
      # 使用矩阵运算计算得分
      for (i in seq_along(cell_types)) {
        cell_i <- cell_types[i]
        genes_ct1 <- intersect(cell_markers[[cell_i]], rownames(sample_expr))
        genes_ct1 <- intersect(genes_ct1, rownames(shortest_path_matrix))

        if (length(genes_ct1) == 0) next

        for (j in seq_along(cell_types)) {
          cell_j <- cell_types[j]
          genes_ct2 <- intersect(cell_markers[[cell_j]], rownames(sample_expr))
          genes_ct2 <- intersect(genes_ct2, rownames(shortest_path_matrix))

          if (length(genes_ct2) == 0) next

          # 提取子矩阵
          dist_submat <- shortest_path_matrix[genes_ct1, genes_ct2, drop = FALSE]
          weight_submat <- weight_matrix[genes_ct1, genes_ct2, drop = FALSE]
          weight_submat_t <- weight_matrix[genes_ct2, genes_ct1, drop = FALSE]

          # 计算表达值外积
          expr_prod <- outer(expr_values[genes_ct1], expr_values[genes_ct2])

          # 计算权重和
          weight_sum <- outer(rep(1, length(genes_ct1)), rep(1, length(genes_ct2))) *
            (weight_submat + t(weight_submat_t))

          # 合并计算
          valid_idx <- !is.infinite(dist_submat) & !is.na(dist_submat)
          if (any(valid_idx)) {
            score_mat <- (expr_prod * weight_sum) / (dist_submat + epsilon)
            path_score_mat[i, j] <- sum(score_mat[valid_idx])
          }
        }
      }

      # 归一化处理
      if (sum(path_score_mat) > 0) {
        path_score_mat <- log1p(path_score_mat)  # 对数转换
        path_score_mat <- (path_score_mat - min(path_score_mat)) /
          (max(path_score_mat) - min(path_score_mat))  # 归一化
      }

      saveRDS(path_score_mat, marker_file)
      marker_networks[[sample]] <- path_score_mat
    }
  }
  # 并行计算表达网络
  if (marker_weight < 1) {
    # 读取PPI权重
    ppi_data <- fread(ppi_weights_file, sep = "\t")
    ppi_weights <- setNames(ppi_data$weight, ppi_data$gene)

    cell_types <- names(cell_markers)

    required_files <- c(ppi_weights_file, inter_path, jacaard_path)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0) {
      stop("以下文件未找到：", paste(missing_files, collapse = ", "))
    }

    # 读取并处理交互矩阵
    matrix_cell_go_inter <- read.csv(inter_path)
    matrix_cell_go_inter <- as.data.frame(matrix_cell_go_inter)
    rownames(matrix_cell_go_inter) <- matrix_cell_go_inter[[1]]
    matrix_cell_go_inter <- matrix_cell_go_inter[,-1]
    matrix_cell_go_inter <- t(matrix_cell_go_inter)

    # 读取并处理Jaccard系数矩阵
    matrix_cell_go_jaccard <- read.csv(jacaard_path)
    matrix_cell_go_jaccard <- as.data.frame(matrix_cell_go_jaccard)
    rownames(matrix_cell_go_jaccard) <- matrix_cell_go_jaccard[[1]]
    matrix_cell_go_jaccard <- matrix_cell_go_jaccard[,-1]
    matrix_cell_go_jaccard <- t(matrix_cell_go_jaccard)

    # 数据预处理
    matrix_cell_go_inter <- apply(matrix_cell_go_inter, c(1, 2), as.character)
    matrix_cell_go_jaccard <- apply(matrix_cell_go_jaccard, c(1, 2), as.numeric)

    # 处理缺失值
    matrix_cell_go_jaccard[is.na(matrix_cell_go_jaccard)] <- 0
    matrix_cell_go_inter[is.na(matrix_cell_go_inter)] <- ""

    # 定义加权基因评分计算函数
    weighted_gene_score_row <- function(a, table, expr, weights, c) {
      score_row <- numeric(c)
      weights_lookup <- weights[names(weights) %in% table]

      for (j in seq_along(a)) {
        if (a[j] != "") {
          genes <- unlist(strsplit(gsub("^\"|\"$", "", a[j]), split = ","))
          location <- match(genes, table)

          valid_indices <- !is.na(location)
          if (any(valid_indices)) {
            valid_genes <- genes[valid_indices]
            valid_expr <- as.numeric(expr[location[valid_indices]])

            valid_weights <- ifelse(valid_genes %in% names(weights_lookup),
                                  weights_lookup[valid_genes],
                                  1)

            if (length(valid_expr) > 0 && sum(valid_weights) > 0) {
              score_row[j] <- sum(valid_expr * (valid_weights / sum(valid_weights)))
            }
          }
        }
      }
      return(score_row)
    }

    # 并行处理每个样本
    expr_networks <- foreach::foreach(sample = sample_names, .packages = c("Matrix", "data.table")) %dopar% {
      # 检查缓存
      expr_file <- file.path(networks_dir, paste0("expression_", sample, ".rds"))
      if (file.exists(expr_file)) {
        cached <- readRDS(expr_file)
        if (!is.matrix(cached)) {
          cached <- as.matrix(cached)
        }
        return(cached)
      }

      # 获取样本数据
      sample_expr <- GEP[, sample, drop = FALSE]
      
      # 预处理GEP数据
      GEPscore <- data.table(gene = rownames(sample_expr),
                            value = as.numeric(sample_expr[, 1]))
      setkey(GEPscore, gene)

      # 计算加权得分
      weighted_score <- Matrix(0, nrow = length(rownames(matrix_cell_go_inter)),
                              ncol = length(colnames(matrix_cell_go_inter)),
                              sparse = TRUE)

      # 批量处理
      batch_size <- 200
      total_rows <- length(rownames(matrix_cell_go_inter))
      num_batches <- ceiling(total_rows / batch_size)

      for (b in 1:num_batches) {
        start_idx <- (b-1) * batch_size + 1
        end_idx <- min(b * batch_size, total_rows)

        for (k in start_idx:end_idx) {
          genes_vector <- matrix_cell_go_inter[k, ]
          row <- weighted_gene_score_row(genes_vector, GEPscore$gene,
                                       GEPscore$value, ppi_weights,
                                       length(colnames(matrix_cell_go_inter)))
          weighted_score[k, ] <- row
        }
      }

      # 计算最终矩阵
      matrix_weighted_genes <- weighted_score * matrix_cell_go_jaccard
      if (!(is.matrix(matrix_weighted_genes) || inherits(matrix_weighted_genes, "Matrix"))) {
        stop(sprintf(
          "expression network build failed for sample %s: matrix_weighted_genes is not a matrix (class=%s)",
          sample, paste(class(matrix_weighted_genes), collapse = ",")
        ))
      }
      expr_cell_matrix <- Matrix::crossprod(matrix_weighted_genes)
      expr_cell_matrix <- as.matrix(expr_cell_matrix)
      expr_cell_matrix[is.na(expr_cell_matrix)] <- 0
      diag(expr_cell_matrix) <- 0
      # 归一化处理
      if (sum(expr_cell_matrix) > 0) {
        expr_cell_matrix <- log1p(expr_cell_matrix)  # 对数转换
        expr_cell_matrix <-  as.matrix((expr_cell_matrix - min(expr_cell_matrix)) /
          (max(expr_cell_matrix) - min(expr_cell_matrix)))  # 归一化
      }
      message(sprintf("处理GO网络: %s", sample))
      # 保存结果
      saveRDS(expr_cell_matrix, expr_file)
      expr_cell_matrix
    }

    # 将结果转换为命名列表
    names(expr_networks) <- sample_names
  }

  message("开始读取并融合网络...")
  # 读取第一个样本的网络以获取维度信息
  first_sample <- sample_names[1]
  first_marker_net <- readRDS(file.path(networks_dir, paste0("marker_", first_sample, ".rds")))
  message(sprintf("第一个样本的标记网络维度: %s", paste(dim(first_marker_net), collapse = "x")))
  message(sprintf("第一个样本的标记网络和: %s", sum(first_marker_net)))

  # 初始化得分矩阵
  Score_rankwalk <- matrix(0, nrow = nrow(first_marker_net), ncol = length(sample_names))
  rownames(Score_rankwalk) <- rownames(first_marker_net)
  colnames(Score_rankwalk) <- sample_names

  # 对每个样本进行处理
  for (i in seq_along(sample_names)) {
    sample <- sample_names[i]
    message(sprintf("\n处理样本: %s", sample))

    # 读取网络文件
    marker_net <- readRDS(file.path(networks_dir, paste0("marker_", sample, ".rds")))
    expr_net <- readRDS(file.path(networks_dir, paste0("expression_", sample, ".rds")))

    message(sprintf("标记网络维度: %s, 和: %s",
                    paste(dim(marker_net), collapse = "x"), sum(marker_net)))
    message(sprintf("表达网络维度: %s, 和: %s",
                    paste(dim(expr_net), collapse = "x"), sum(expr_net)))

    # 融合网络
    fused_network <- if (marker_weight == 1) {
      marker_net
    } else if (marker_weight == 0) {
      expr_net
    } else {
      # 确保两个网络具有相同的维度
      common_cells <- intersect(rownames(marker_net), rownames(expr_net))
      message(sprintf("共同细胞类型数量: %d", length(common_cells)))

      marker_submat <- marker_net[common_cells, common_cells]
      expr_submat <- expr_net[common_cells, common_cells]

      message(sprintf("子矩阵维度: %s", paste(dim(marker_submat), collapse = "x")))
      message(sprintf("标记子网络和: %s", sum(marker_submat)))
      message(sprintf("表达子网络和: %s", sum(expr_submat)))

      # 融合两个网络
      combined_net <- (1 - marker_weight) * expr_submat + marker_weight * marker_submat
      message(sprintf("融合网络和: %s", sum(combined_net)))

      # 确保网络值在合理范围内
      if (sum(combined_net) > 0) {
        combined_net <- combined_net / max(combined_net)
        message(sprintf("归一化后融合网络和: %s", sum(combined_net)))
      }
      combined_net
    }

    # 随机游走计算得分
    tryCatch({
      # 创建图对象前检查网络
      message(sprintf("融合网络维度: %s", paste(dim(fused_network), collapse = "x")))
      message(sprintf("融合网络非零元素数量: %d", sum(fused_network != 0)))

      graph <- graph_from_adjacency_matrix(
        fused_network,
        mode = "undirected",
        weighted = TRUE
      )

      message(sprintf("图的边数: %d", ecount(graph)))
      message(sprintf("图的节点数: %d", vcount(graph)))

      # 计算PageRank
      pr <- page_rank(graph, directed = FALSE, damping = damping)
      message(sprintf("PageRank得分和: %s", sum(pr$vector)))
      message(sprintf("PageRank得分范围: [%s, %s]",
                      min(pr$vector), max(pr$vector)))

      Score_rankwalk[rownames(fused_network), i] <- pr$vector

    }, error = function(e) {
      warning(sprintf("处理样本 %s 时发生错误: %s", sample, e$message))
    })

    # 清理内存
    rm(marker_net, expr_net, fused_network, graph)
    gc()
  }

  # 创建进度条
  pb <- progress_bar$new(
    total = length(sample_names),
    format = "计算融合得分进度 [:bar] :percent 已完成 (用时: :elapsed)",
    clear = FALSE,
    width = 60
  )

  # 检查最终得分
  if (all(Score_rankwalk == 0)) {
    stop("所有样本的得分都为0，请检查网络计算过程")
  }

  # 对数转换和归一化
  Score_rankwalk <- apply(Score_rankwalk, 2, function(x) {
    if (all(x <= 0)) {
      warning("检测到非正值，将使用小的正数替代")
      x <- rep(.Machine$double.eps, length(x))
    } else {
      x[x <= 0] <- min(x[x > 0]) / 2
    }
    return(log(x + 1))  # 使用自然对数，也可以指定其他底数
  })

  # 最小-最大归一化
  score_range <- range(Score_rankwalk, na.rm = TRUE, finite = TRUE)
  if (diff(score_range) > 0) {
    Score_rankwalk <- (Score_rankwalk - score_range[1]) / diff(score_range)
  }

  return(Score_rankwalk)

}
