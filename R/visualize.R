toBin <- function(n, nG) {
  x = paste(as.integer(rev(intToBits(n))), collapse = "")
  substr(x, nchar(x) - nG + 1, nchar(x))
}

#' Visualize CBN Model
#'
#' @param fitness Fitness vectors for each genotype provided in selectNodes or for all genotypes if none selected
#' @param selectNodes Select genotypes to visualize
#' @param nGenes Length of each genotype
#' @param lowColor Color for wild type genotype
#' @param highColor Color for fully mutated genotype
#'
#' @return Plot visualization of CBN model
#' @export
#'
#' @examples
#' Genotypes<-c(
#' "0000",
#' "1000",
#' "0100",
#' "0010",
#' "0001",
#' "1100",
#' "1010",
#' "1001",
#' "0110",
#' "0101",
#' "0011",
#' "1110",
#' "1101",
#' "1011",
#' "0111",
#' "1111"
#' )
# 
#' COLintensity<-c(0,rep(0.25,4),rep(0.5,6),rep(0.75,4),1)
#' visualize(COLintensity)
visualize_model <- function(fitness, selectNodes = NULL, nGenes = 4, lowColor = "lightblue", highColor = "pink") {
  allStrings = lapply(0:(2^nGenes - 1), function(x) toBin(x, nGenes))
  
  count_ones <- function(s) sum(unlist(gregexpr("1", s)) > 0)
  getColumn <- function(num) rev(allStrings[sapply(allStrings, count_ones) == num])
  columns = lapply(0:nGenes, getColumn)
  
  # Construct edges
  edges = list()
  for (i in 1:(length(allStrings) - 1)) {
    for (j in (i + 1):length(allStrings)) {
      node1 = allStrings[[i]]
      node2 = allStrings[[j]]
      if (sum(strsplit(node1, NULL)[[1]] != strsplit(node2, NULL)[[1]]) == 1 & ((node1 %in% selectNodes & node2 %in% selectNodes) || is.null(selectNodes))) {
        edges = append(edges, list(c(node1, node2)))
      }
    }
  }
  
  nodes <- unlist(columns)
  
  # Create data frame for nodes and positions
  layout_df <- data.frame(name = nodes, stringsAsFactors = FALSE)
  
  # Assign x/y layout
  col_start <- 1
  x_vals <- numeric(length(nodes))
  y_vals <- numeric(length(nodes))
  
  for (column_nodes in columns) {
    column_height <- length(column_nodes)
    
    ys <- seq((column_height - 1) / 2, -(column_height - 1) / 2)
    
    for (j in seq_along(column_nodes)) {
      idx <- which(layout_df$name == column_nodes[j])
      x_vals[idx] <- col_start
      y_vals[idx] <- ys[j]
    }
    
    col_start <- col_start + 1
  }
  
  getFitness <- function(name) {
    if (is.null(selectNodes)) {
      fitness[match(name, nodes)]
    } else if (name %in% selectNodes) {
      fitness[match(name, selectNodes)]
    } else {
      NA
    }
  }
  
  layout_df$x <- x_vals
  layout_df$y <- y_vals
  layout_df$fitness <- unlist(lapply(layout_df$name, getFitness))
  
  # Convert edge list to data frame
  edge_df <- do.call(rbind, edges)
  colnames(edge_df) <- c("from", "to")
  edge_df <- as.data.frame(edge_df, stringsAsFactors = FALSE)
  
  # Create tidygraph object
  g_tbl <- tbl_graph(nodes = layout_df, edges = edge_df, directed = FALSE)
  
  # Plot with ggraph
  ggraph(g_tbl, layout = "manual", x = x, y = y) +
    geom_edge_link(color = "gray80") +
    geom_node_point(aes(fill = fitness), shape = 21, size = 12, stroke = 0.3, color="black") +
    geom_node_text(aes(label = name), size = 3) +
    scale_fill_gradient(low = lowColor, high = highColor, na.value = "white", name = "Fitness") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title.align = 0.5,
      legend.title.position = "top"
    ) +
    ggtitle("Fitness Landscape")
}