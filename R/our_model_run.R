#' Run the full CrosstalkIn workflow
#'
#' Execute the complete CrosstalkIn pipeline, including PPI processing,
#' PageRank-based PPI weight calculation, network fusion, and final score
#' generation.
#'
#' @param links_file Path to the STRING links file.
#' @param protein_info_file Path to the STRING protein info file.
#' @param expr_rds Path to an RDS file containing the expression matrix, or
#'   an object with an \code{expr} element.
#' @param cellmarkers_csv Path to the cell marker CSV file.
#' @param shortest_paths_csv Path to the shortest path matrix CSV file.
#' @param inter_path Path to the GO interaction file.
#' @param jacaard_path Path to the GO Jaccard file.
#' @param marker_weight Numeric fusion weight in \code{[0, 1]}.
#' @param score_threshold Minimum STRING interaction score to retain.
#' @param damping_factor Damping parameter for PPI PageRank weighting.
#' @param damping Damping parameter for PageRank on fused cell-cell networks.
#' @param networks_dir Directory for intermediate network files.
#' @param cl.cores Number of CPU cores for parallel computation.
#' @param output_dir Directory for output files.
#' @param save_scores_to Optional path to save final scores as an RDS file.
#'
#' @return A numeric matrix of fused CrosstalkIn scores.
#' @export
our_model_run <- function(links_file,
                          protein_info_file,
                          expr_rds,
                          cellmarkers_csv,
                          shortest_paths_csv,
                          inter_path,
                          jacaard_path,
                          marker_weight = 0.5,
                          score_threshold = 700,
                          damping_factor = 0.85,
                          damping = 0.90,
                          networks_dir = file.path(tempdir(), "CrosstalkIn_networks"),
                          cl.cores = max(1, parallel::detectCores() - 1),
                          output_dir = file.path(tempdir(), "ourmodel_outputs"),
                          save_scores_to = NULL) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!file.exists(expr_rds)) {
    stop("expr_rds not exist: ", expr_rds)
  }
  expr_obj <- readRDS(expr_rds)
  if (is.list(expr_obj) && "expr" %in% names(expr_obj)) {
    GEP <- expr_obj$expr
  } else {
    GEP <- expr_obj
  }
  if (!is.matrix(GEP) && !is.data.frame(GEP)) {
    stop("expr_rds need dataframe or have $expr")
  }
  GEP <- as.matrix(GEP)

  if (is.null(rownames(GEP)) || any(rownames(GEP) == "")) {
    stop("expression need gene names(rownames)")
  }
  if (is.null(colnames(GEP)) || any(colnames(GEP) == "")) {
    stop("expression need sample names(colnames)")
  }

  if (!file.exists(cellmarkers_csv)) {
    stop("cellmarkers_csv not exist: ", cellmarkers_csv)
  }
  cellmarkers_df <- utils::read.csv(cellmarkers_csv, row.names = 1, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(cellmarkers_df) < 1) {
    stop("cellmarkers_csv need one col")
  }
  cell_markers <- lapply(rownames(cellmarkers_df), function(cell_type) {
    markers_str <- as.character(cellmarkers_df[cell_type, 1])
    markers <- unlist(strsplit(markers_str, ",", fixed = TRUE))
    markers <- trimws(markers)
    markers[nzchar(markers)]
  })
  names(cell_markers) <- rownames(cellmarkers_df)

  ppi_network_file <- file.path(output_dir, "ppi_network.txt")
  ppi_weights_file <- file.path(output_dir, "ppi_weight.txt")

  ppi_network <- process_string_ppi_dt(
    links_file = links_file,
    protein_info_file = protein_info_file,
    output_file = ppi_network_file,
    score_threshold = score_threshold
  )

  calculate_ppi_weights(
    PPI_network = ppi_network,
    gene_list = rownames(GEP),
    output_file = ppi_weights_file,
    damping_factor = damping_factor
  )

  fused_scores <- fuse_networks(
    GEP = GEP,
    cell_markers = cell_markers,
    sample_names = colnames(GEP),
    marker_weight = marker_weight,
    damping = damping,
    networks_dir = networks_dir,
    ppi_weights_file = ppi_weights_file,
    inter_path = inter_path,
    jacaard_path = jacaard_path,
    shortest_paths_csv = shortest_paths_csv,
    cl.cores = cl.cores
  )

  if (!is.null(save_scores_to)) {
    out_dir <- dirname(save_scores_to)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    saveRDS(fused_scores, save_scores_to)
  }

  fused_scores
}
