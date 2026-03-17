#' Calculate PPI-based gene weights using PageRank
#'
#' Filter a protein-protein interaction network to a target gene set and
#' compute gene weights using the PageRank algorithm.
#'
#' @param PPI_network A data.frame containing columns \code{gene1},
#'   \code{gene2}, and \code{weight}.
#' @param gene_list A character vector of genes of interest.
#' @param output_file Optional output path for saving the weights.
#' @param damping_factor Damping parameter used in PageRank.
#'
#' @return A named numeric vector of normalized gene weights.
#' @export
calculate_ppi_weights <- function(PPI_network, gene_list, output_file = NULL, damping_factor = 0.85) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("need igraph package：install.packages('igraph')")
  }
  if (!is.data.frame(PPI_network)) {
    stop("PPI_network must data.frame")
  }
  if (!all(c("gene1", "gene2", "weight") %in% colnames(PPI_network))) {
    stop("PPI_network must need col: gene1, gene2, weight")
  }
  if (!is.character(gene_list)) {
    stop("gene_list must number must be a character vector")
  }
  if (!is.numeric(damping_factor) || length(damping_factor) != 1 || is.na(damping_factor)) {
    stop("damping_factor must be a single value.")
  }

  message("Filtering PPI network to genes in gene list...")
  ppi_filtered <- PPI_network[PPI_network$gene1 %in% gene_list & PPI_network$gene2 %in% gene_list, , drop = FALSE]

  if (nrow(ppi_filtered) == 0) {
    message("No interactions found between genes in the gene list. Using uniform weights.")
    weights <- rep(1, length(gene_list))
    names(weights) <- gene_list
    weights <- weights / sum(weights)

    if (!is.null(output_file)) {
      utils::write.table(
        data.frame(gene = names(weights), weight = weights),
        file = output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
    return(weights)
  }

  message("Creating graph from PPI network...")
  ppi_graph <- igraph::graph_from_data_frame(ppi_filtered, directed = FALSE)

  components <- igraph::components(ppi_graph)
  message("Graph contains ", components$no, " connected components")
  message("Largest component has ", max(components$csize), " genes")

  message("Running PageRank algorithm to calculate gene weights...")
  pr_result <- igraph::page_rank(
    ppi_graph,
    directed = FALSE,
    damping = damping_factor,
    weights = igraph::E(ppi_graph)$weight
  )

  all_weights <- rep(0, length(gene_list))
  names(all_weights) <- gene_list
  genes_in_network <- names(pr_result$vector)
  all_weights[genes_in_network] <- pr_result$vector

  positive_pr <- pr_result$vector[pr_result$vector > 0]
  baseline <- if (length(positive_pr) > 0) min(positive_pr) / 10 else .Machine$double.eps
  all_weights[all_weights == 0] <- baseline
  all_weights <- all_weights / sum(all_weights)

  message("Weight calculation complete.")
  message("Min weight: ", min(all_weights))
  message("Max weight: ", max(all_weights))
  message("Mean weight: ", mean(all_weights))

  if (!is.null(output_file)) {
    utils::write.table(
      data.frame(gene = names(all_weights), weight = all_weights),
      file = output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message("Weights saved to ", output_file)
  }

  all_weights
}
