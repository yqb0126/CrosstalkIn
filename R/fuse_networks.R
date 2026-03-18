#' Fuse marker-derived and expression-derived cell-cell networks
#'
#' Construct marker-based and/or expression-based cell-cell interaction
#' networks for each sample, fuse them with a user-defined weight, and
#' compute sample-specific infiltration scores using PageRank.
#'
#' @param expr_networks Optional precomputed expression-derived networks.
#' @param marker_networks Optional precomputed marker-derived networks.
#' @param GEP Gene expression matrix with genes in rows and samples in columns.
#' @param cell_markers A named list of marker genes for each cell type.
#' @param sample_names Optional sample names. Defaults to \code{colnames(GEP)}.
#' @param marker_weight Numeric fusion weight in \code{[0, 1]}.
#' @param weighted Logical; currently retained for interface compatibility.
#' @param damping Damping factor used in PageRank on fused networks.
#' @param networks_dir Directory for caching intermediate sample-specific networks.
#' @param ppi_weights_file Path to the PPI weight file.
#' @param inter_path Path to the GO interaction file.
#' @param jacaard_path Path to the GO Jaccard file.
#' @param shortest_paths_csv Path to the shortest path matrix file.
#' @param cl.cores Number of CPU cores for parallel computation.
#'
#' @return A numeric matrix of normalized sample-specific scores, with
#'   cell types in rows and samples in columns.
#' @importFrom foreach foreach %dopar%
#' @export
fuse_networks <- function(expr_networks = NULL, marker_networks = NULL,
                          GEP = NULL, cell_markers = NULL,
                          sample_names = NULL, marker_weight = 1,
                          weighted = TRUE, damping = 0.90,
                          networks_dir = file.path(tempdir(), "CrosstalkIn_networks"),
                          ppi_weights_file = NULL,
                          inter_path = NULL,
                          jacaard_path = NULL,
                          shortest_paths_csv = NULL,
                          cl.cores = max(1, parallel::detectCores() - 1)) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it with install.packages('data.table').")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required. Please install it with install.packages('Matrix').")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Please install it with install.packages('igraph').")
  }
  if (is.null(GEP) || is.null(cell_markers)) {
    stop("Both 'GEP' and 'cell_markers' must be provided.")
  }
  if (is.null(sample_names)) {
    sample_names <- colnames(GEP)
  }
  if (length(sample_names) == 0) {
    stop("sample_names is empty")
  }
  if (!is.numeric(marker_weight) || length(marker_weight) != 1 || is.na(marker_weight) ||
      marker_weight < 0 || marker_weight > 1) {
    stop("marker_weight in [0, 1]")
  }
  if (!dir.exists(networks_dir)) {
    dir.create(networks_dir, recursive = TRUE, showWarnings = FALSE)
  }

  make_progress <- function(total) {
    if (requireNamespace("progress", quietly = TRUE)) {
      return(progress::progress_bar$new(
        total = total,
        format = "format = Progress [:bar] :percent (elapsed: :elapsed)",
        clear = FALSE,
        width = 60
      ))
    }
    NULL
  }

  if (marker_weight > 0) {
    if (is.null(shortest_paths_csv) || !file.exists(shortest_paths_csv)) {
      stop("marker_weight > 0 need shortest_paths_csv pathway")
    }
    if (is.null(ppi_weights_file) || !file.exists(ppi_weights_file)) {
      stop("marker_weight > 0 need ppi_weights_file pathway")
    }

    shortest_path_df <- utils::read.csv(shortest_paths_csv, check.names = FALSE, stringsAsFactors = FALSE)
    names_vec <- shortest_path_df[[1]]
    shortest_path_matrix <- as.matrix(shortest_path_df[, -1, drop = FALSE])
    rownames(shortest_path_matrix) <- colnames(shortest_path_matrix) <- names_vec
    shortest_path_matrix <- matrix(
      as.numeric(shortest_path_matrix),
      nrow = nrow(shortest_path_matrix),
      dimnames = list(rownames(shortest_path_matrix), colnames(shortest_path_matrix))
    )

    ppi_data <- data.table::fread(ppi_weights_file, sep = "\t", showProgress = FALSE)
    if (!all(c("gene", "weight") %in% names(ppi_data))) {
      stop("ppi_weights_file need include gene and weight ")
    }
    ppi_weights <- stats::setNames(ppi_data$weight, ppi_data$gene)

    all_genes <- unique(c(names(ppi_weights), rownames(shortest_path_matrix)))
    weight_matrix <- matrix(1, nrow = length(all_genes), ncol = length(all_genes))
    rownames(weight_matrix) <- colnames(weight_matrix) <- all_genes
    weight_indices <- match(names(ppi_weights), all_genes)
    weight_matrix[weight_indices, ] <- ppi_weights
    weight_matrix[, weight_indices] <- ppi_weights

    epsilon <- 1e-6
    cell_types <- names(cell_markers)
    marker_networks <- list()
    pb <- make_progress(length(sample_names))

    for (sample in sample_names) {
      marker_file <- file.path(networks_dir, paste0("marker_", sample, ".rds"))
      if (file.exists(marker_file)) {
        marker_networks[[sample]] <- readRDS(marker_file)
        if (!is.null(pb)) pb$tick()
        next
      }

      sample_expr <- GEP[, sample, drop = FALSE]
      if (ncol(sample_expr) == 0 || nrow(sample_expr) == 0) {
        stop(sprintf("expression is empty", sample))
      }

      path_score_mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types))
      rownames(path_score_mat) <- colnames(path_score_mat) <- cell_types

      expr_values <- as.numeric(sample_expr[rownames(sample_expr), 1])
      names(expr_values) <- rownames(sample_expr)

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

          dist_submat <- shortest_path_matrix[genes_ct1, genes_ct2, drop = FALSE]
          weight_submat <- weight_matrix[genes_ct1, genes_ct2, drop = FALSE]
          weight_submat_t <- weight_matrix[genes_ct2, genes_ct1, drop = FALSE]

          expr_prod <- outer(expr_values[genes_ct1], expr_values[genes_ct2])
          weight_sum <- outer(rep(1, length(genes_ct1)), rep(1, length(genes_ct2))) * (weight_submat + t(weight_submat_t))

          valid_idx <- !is.infinite(dist_submat) & !is.na(dist_submat)
          if (any(valid_idx)) {
            score_mat <- (expr_prod * weight_sum) / (dist_submat + epsilon)
            path_score_mat[i, j] <- sum(score_mat[valid_idx])
          }
        }
      }

      if (sum(path_score_mat) > 0) {
        path_score_mat <- log1p(path_score_mat)
        rng <- range(path_score_mat)
        if (diff(rng) > 0) {
          path_score_mat <- (path_score_mat - rng[1]) / diff(rng)
        }
      }

      saveRDS(path_score_mat, marker_file)
      marker_networks[[sample]] <- path_score_mat
      if (!is.null(pb)) pb$tick()
    }
  }

  if (marker_weight < 1) {
    if (is.null(ppi_weights_file) || !file.exists(ppi_weights_file)) {
      stop("marker_weight < 1 need ppi_weights_file")
    }
    if (is.null(inter_path) || !file.exists(inter_path)) {
      stop("marker_weight < 1 need inter_path")
    }
    if (is.null(jacaard_path) || !file.exists(jacaard_path)) {
      stop("marker_weight < 1 need jacaard_path")
    }

    ppi_data <- data.table::fread(ppi_weights_file, sep = "\t", showProgress = FALSE)
    if (!all(c("gene", "weight") %in% names(ppi_data))) {
      stop("ppi_weights_file need gene and weight")
    }
    ppi_weights <- stats::setNames(ppi_data$weight, ppi_data$gene)

    matrix_cell_go_inter <- utils::read.csv(inter_path, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(matrix_cell_go_inter) <- matrix_cell_go_inter[[1]]
    matrix_cell_go_inter <- matrix_cell_go_inter[, -1, drop = FALSE]
    matrix_cell_go_inter <- t(as.matrix(matrix_cell_go_inter))

    matrix_cell_go_jaccard <- utils::read.csv(jacaard_path, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(matrix_cell_go_jaccard) <- matrix_cell_go_jaccard[[1]]
    matrix_cell_go_jaccard <- matrix_cell_go_jaccard[, -1, drop = FALSE]
    matrix_cell_go_jaccard <- t(as.matrix(matrix_cell_go_jaccard))

    matrix_cell_go_inter <- apply(matrix_cell_go_inter, c(1, 2), as.character)
    matrix_cell_go_jaccard <- apply(matrix_cell_go_jaccard, c(1, 2), as.numeric)
    matrix_cell_go_jaccard[is.na(matrix_cell_go_jaccard)] <- 0
    matrix_cell_go_inter[is.na(matrix_cell_go_inter)] <- ""

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
            valid_weights <- ifelse(valid_genes %in% names(weights_lookup), weights_lookup[valid_genes], 1)
            if (length(valid_expr) > 0 && sum(valid_weights) > 0) {
              score_row[j] <- sum(valid_expr * (valid_weights / sum(valid_weights)))
            }
          }
        }
      }
      score_row
    }

    expr_file_for <- function(sample) file.path(networks_dir, paste0("expression_", sample, ".rds"))

    run_one_sample_expr <- function(sample) {
      expr_file <- expr_file_for(sample)
      if (file.exists(expr_file)) {
        cached <- readRDS(expr_file)
        if (!is.matrix(cached)) {
          cached <- as.matrix(cached)
        }
        return(cached)
      }

      sample_expr <- GEP[, sample, drop = FALSE]
      GEPscore <- data.table::data.table(gene = rownames(sample_expr), value = as.numeric(sample_expr[, 1]))
      data.table::setkey(GEPscore, gene)

      weighted_score <- Matrix::Matrix(0, nrow = nrow(matrix_cell_go_inter), ncol = ncol(matrix_cell_go_inter), sparse = TRUE)

      batch_size <- 200
      total_rows <- nrow(matrix_cell_go_inter)
      num_batches <- ceiling(total_rows / batch_size)

      for (b in 1:num_batches) {
        start_idx <- (b - 1) * batch_size + 1
        end_idx <- min(b * batch_size, total_rows)
        for (k in start_idx:end_idx) {
          genes_vector <- matrix_cell_go_inter[k, ]
          row <- weighted_gene_score_row(genes_vector, GEPscore$gene, GEPscore$value, ppi_weights, ncol(matrix_cell_go_inter))
          weighted_score[k, ] <- row
        }
      }

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

      if (sum(expr_cell_matrix) > 0) {
        expr_cell_matrix <- log1p(expr_cell_matrix)
        rng <- range(expr_cell_matrix)
        if (diff(rng) > 0) {
          expr_cell_matrix <- as.matrix((expr_cell_matrix - rng[1]) / diff(rng))
        }
      }

      saveRDS(expr_cell_matrix, expr_file)
      expr_cell_matrix
    }

    can_parallel <- cl.cores > 1 &&
      requireNamespace("foreach", quietly = TRUE) &&
      requireNamespace("doParallel", quietly = TRUE)

    if (can_parallel) {
      cl <- parallel::makeCluster(cl.cores)
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      doParallel::registerDoParallel(cl)

      expr_networks <- foreach::foreach(sample = sample_names, .packages = c("Matrix", "data.table")) %dopar% {
        run_one_sample_expr(sample)
      }
    } else {
      expr_networks <- lapply(sample_names, run_one_sample_expr)
    }
    names(expr_networks) <- sample_names
  }

  first_sample <- sample_names[1]
  first_net <- NULL
  if (marker_weight > 0) {
    first_path <- file.path(networks_dir, paste0("marker_", first_sample, ".rds"))
    if (!file.exists(first_path)) stop("need marker network file: ", first_path)
    first_net <- readRDS(first_path)
  } else {
    first_path <- file.path(networks_dir, paste0("expression_", first_sample, ".rds"))
    if (!file.exists(first_path)) stop("need expression network file: ", first_path)
    first_net <- readRDS(first_path)
  }

  Score_rankwalk <- matrix(0, nrow = nrow(first_net), ncol = length(sample_names))
  rownames(Score_rankwalk) <- rownames(first_net)
  colnames(Score_rankwalk) <- sample_names

  pb <- make_progress(length(sample_names))
  for (i in seq_along(sample_names)) {
    sample <- sample_names[i]

    marker_net <- NULL
    expr_net <- NULL
    if (marker_weight > 0) {
      marker_net <- readRDS(file.path(networks_dir, paste0("marker_", sample, ".rds")))
    }
    if (marker_weight < 1) {
      expr_net <- readRDS(file.path(networks_dir, paste0("expression_", sample, ".rds")))
    }

    fused_network <- if (marker_weight == 1) {
      marker_net
    } else if (marker_weight == 0) {
      expr_net
    } else {
      common_cells <- intersect(rownames(marker_net), rownames(expr_net))
      marker_submat <- marker_net[common_cells, common_cells, drop = FALSE]
      expr_submat <- expr_net[common_cells, common_cells, drop = FALSE]

      combined_net <- (1 - marker_weight) * expr_submat + marker_weight * marker_submat
      if (sum(combined_net) > 0) {
        combined_net <- combined_net / max(combined_net)
      }
      combined_net
    }

    tryCatch({
      graph <- igraph::graph_from_adjacency_matrix(fused_network, mode = "undirected", weighted = TRUE)
      pr <- igraph::page_rank(graph, directed = FALSE, damping = damping)
      Score_rankwalk[rownames(fused_network), i] <- pr$vector
    }, error = function(e) {
      warning(sprintf("sample %s is error: %s", sample, e$message))
    })

    if (!is.null(pb)) pb$tick()
  }

  if (all(Score_rankwalk == 0)) {
    stop("all score is 0, please cheak")
  }

  Score_rankwalk <- apply(Score_rankwalk, 2, function(x) {
    if (all(x <= 0)) {
      x <- rep(.Machine$double.eps, length(x))
    } else {
      x[x <= 0] <- min(x[x > 0]) / 2
    }
    log(x + 1)
  })

  score_range <- range(Score_rankwalk, na.rm = TRUE, finite = TRUE)
  if (diff(score_range) > 0) {
    Score_rankwalk <- (Score_rankwalk - score_range[1]) / diff(score_range)
  }

  Score_rankwalk
}
