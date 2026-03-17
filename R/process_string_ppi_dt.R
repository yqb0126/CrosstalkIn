#' Process STRING PPI data and map proteins to genes
#'
#' Read STRING interaction and protein annotation files, filter interactions
#' by score, map protein identifiers to gene symbols, and generate a
#' gene-gene weighted interaction table.
#'
#' @param links_file Path to the STRING links file.
#' @param protein_info_file Path to the STRING protein info file.
#' @param output_file Optional output path for saving the processed network.
#' @param score_threshold Numeric threshold for filtering STRING combined scores.
#'
#' @return A data.frame with columns \code{gene1}, \code{gene2}, and \code{weight}.
#' @export
process_string_ppi_dt <- function(links_file, protein_info_file, output_file = NULL, score_threshold = 700) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("need data.table：install.packages('data.table')")
  }
  if (!file.exists(links_file)) {
    stop("PPI network not exist: ", links_file)
  }
  if (!file.exists(protein_info_file)) {
    stop("PPI imformation not exist: ", protein_info_file)
  }
  if (!is.numeric(score_threshold) || length(score_threshold) != 1 || is.na(score_threshold)) {
    stop("score_threshold must number")
  }

  message("Reading STRING interaction file...")
  links_data <- data.table::fread(links_file, showProgress = TRUE)
  if (!"combined_score" %in% names(links_data)) {
    stop("links_file need combined_score colmumn")
  }
  if (!all(c("protein1", "protein2") %in% names(links_data))) {
    stop("links_file need protein1/protein2 colmumn")
  }

  message("Filtering interactions with score >= ", score_threshold, " ...")
  links_data <- links_data[links_data[["combined_score"]] >= score_threshold, ]
  message("Retained ", nrow(links_data), " PPI interactions")

  message("read PPI imformation...")
  protein_info <- data.table::fread(protein_info_file, showProgress = TRUE)

  id_col <- NULL
  if ("protein_external_id" %in% names(protein_info)) {
    id_col <- "protein_external_id"
  } else if ("X.string_protein_id" %in% names(protein_info)) {
    id_col <- "X.string_protein_id"
  } else if ("#string_protein_id" %in% names(protein_info)) {
    id_col <- "#string_protein_id"
  } else if ("string_protein_id" %in% names(protein_info)) {
    id_col <- "string_protein_id"
  }

  if (is.null(id_col) || !"preferred_name" %in% names(protein_info)) {
    stop("protein_info_file must contain one of protein_external_id, X.string_protein_id, #string_protein_id, or string_protein_id, together with preferred_name")
  }

  protein_mapping <- data.frame(
    protein_id = protein_info[[id_col]],
    preferred_name = protein_info[["preferred_name"]],
    stringsAsFactors = FALSE
  )

  message("we read: ", nrow(protein_mapping), " gene information")

  message("change PPI ID to gene ID...")

  id_to_gene <- protein_mapping[["preferred_name"]]
  names(id_to_gene) <- protein_mapping[["protein_id"]]

  links_data[["gene1"]] <- id_to_gene[links_data[["protein1"]]]
  links_data[["gene2"]] <- id_to_gene[links_data[["protein2"]]]

  na_count <- sum(is.na(links_data[["gene1"]]) | is.na(links_data[["gene2"]]))
  if (na_count > 0) {
    message("Warning: ", na_count, " interactions could not be mapped to gene symbols")
    keep_idx <- !is.na(links_data[["gene1"]]) & !is.na(links_data[["gene2"]])
    links_data <- links_data[keep_idx, , drop = FALSE]
    message("Retained: ", nrow(links_data), " mapped interactions")
  }

  links_data[["weight"]] <- links_data[["combined_score"]] / 1000
  links_data[["need_swap"]] <- links_data[["gene1"]] > links_data[["gene2"]]

  if (any(links_data[["need_swap"]], na.rm = TRUE)) {
    idx <- which(links_data[["need_swap"]])
    tmp1 <- links_data[["gene1"]][idx]
    tmp2 <- links_data[["gene2"]][idx]
    links_data[["gene1"]][idx] <- tmp2
    links_data[["gene2"]][idx] <- tmp1
  }

  links_data[["need_swap"]] <- NULL
  links_data <- links_data[, c("gene1", "gene2", "weight"), drop = FALSE]

  ppi_network <- aggregate(weight ~ gene1 + gene2, data = links_data, FUN = max)
  message("Finally network save: ", nrow(ppi_network), " gene relationship")

  if (!is.null(output_file)) {
    out_dir <- dirname(output_file)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    message("save network to file: ", output_file)
    data.table::fwrite(ppi_network, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  as.data.frame(ppi_network)
}
