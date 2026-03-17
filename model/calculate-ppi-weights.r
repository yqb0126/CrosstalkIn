#' @title calculate_ppi_weights
#' @description Calculate gene weights based on PPI network using random walk algorithm and save them to file
#' @param PPI_network A data frame containing the PPI network information with columns: gene1, gene2, and weight
#' @param gene_list Character vector containing gene names for which weights should be calculated
#' @param output_file Path to save the resulting gene weights (default: "ppi_weight.txt")
#' @param damping_factor Restart probability for the random walk algorithm (default: 0.85)
#' @return A named vector of gene weights (also saved to the specified file)
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph page_rank
#' @export
#' @examples
#' # Load PPI network from STRING database
#' PPI_network <- read.csv("string_human_ppi.csv")
#' # Get gene list from expression data
#' gene_list <- rownames(GEP)
#' # Calculate and save weights
#' gene_weights <- calculate_ppi_weights(PPI_network, gene_list, "ppi_weight.txt", 0.85)

calculate_ppi_weights <- function(PPI_network, gene_list, output_file = "ppi_weight.txt", damping_factor = 0.85) {
  # Input validation
  #检查 PPI_network 是否为数据框，并且包含 gene1、gene2 和 weight 列。
  if (!is.data.frame(PPI_network)) {
    stop("PPI_network must be a data frame")
  }
  if (!all(c("gene1", "gene2", "weight") %in% colnames(PPI_network))) {
    stop("PPI_network must contain columns: gene1, gene2, and weight")
  }
  #检查 gene_list 是否为字符向量
  if (!is.character(gene_list)) {
    stop("gene_list must be a character vector")
  }
  # 仅保留基因列表中的基因对。
  # Filter PPI network to include only genes in our gene list
  message("Filtering PPI network to genes in gene list...")
  ppi_filtered <- PPI_network[PPI_network$gene1 %in% gene_list & 
                             PPI_network$gene2 %in% gene_list, ]
  
  # If no edges remain, return uniform weights
  #如果过滤后没有剩余的互作关系，则返回均匀权重。
  if (nrow(ppi_filtered) == 0) {
    message("No interactions found between genes in the gene list. Using uniform weights.")
    weights <- rep(1, length(gene_list))
    names(weights) <- gene_list
    
    # Save to file
    write.table(data.frame(gene = names(weights), weight = weights), 
                file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    return(weights)
  }
  
  # Create graph from PPI network
  message("Creating graph from PPI network...")
  # 使用 igraph 
  ppi_graph <- igraph::graph_from_data_frame(ppi_filtered, directed = FALSE)
  
  # Check for disconnected components and print information
  #检查图中的连通分量数量和最大连通分量的大小。
  components <- igraph::components(ppi_graph)
  message("Graph contains ", components$no, " connected components")
  message("Largest component has ", max(components$csize), " genes")
  
  # Run PageRank algorithm to get importance weights
  message("Running PageRank algorithm to calculate gene weights...")
  # 对PPI进行随机游走，得到每个节点基因的权重
  pr_result <- igraph::page_rank(ppi_graph, directed = FALSE, 
                                damping = damping_factor, 
                                weights = E(ppi_graph)$weight)
  
  # Get weights for all genes in the gene list
  #为所有基因初始化权重为 0。
  all_weights <- rep(0, length(gene_list))
  names(all_weights) <- gene_list
  
  # Fill in weights for genes in the PPI network
  #为 PPI 网络中的基因填充 PageRank 权重
  genes_in_network <- names(pr_result$vector)
  all_weights[genes_in_network] <- pr_result$vector
  
  # For genes not in network, use a small baseline weight
  #对于未出现在网络中的基因，分配一个较小的基线权重。
  baseline <- min(pr_result$vector[pr_result$vector > 0]) / 10
  all_weights[all_weights == 0] <- baseline
  
  # Normalize weights to sum to 1
  #将权重归一化，使其总和为 1。
  all_weights <- all_weights / sum(all_weights)
  
  # Print summary statistics
  message("Weight calculation complete.")
  message("Min weight: ", min(all_weights))
  message("Max weight: ", max(all_weights))
  message("Mean weight: ", mean(all_weights))
  
  # Save to file
  #保存结果
  write.table(data.frame(gene = names(all_weights), weight = all_weights), 
              file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Weights saved to ", output_file)
  
  return(all_weights)
}
