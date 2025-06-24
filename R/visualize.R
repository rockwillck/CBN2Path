toBin <- function(n, nG) {
  x = paste(as.integer(rev(intToBits(n))), collapse = "")
  substr(x, nchar(x) - nG + 1, nchar(x))
}

#' Visualize Fitness Landscape
#'
#' @param fitness Fitness vectors for each genotype provided in selectNodes or for all genotypes if none selected
#' @param selectNodes Select genotypes to visualize
#' @param nGenes Length of each genotype
#' @param lowColor Color for wild type genotype
#' @param highColor Color for fully mutated genotype
#'
#' @return Plot (gg object) visualization of fitness landscape
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
#' visualize_fitness_landscape(COLintensity)
visualize_fitness_landscape <- function(fitness,
                                        selectNodes = NULL,
                                        nGenes = 4,
                                        lowColor = "lightblue",
                                        highColor = "blue") {
  allStrings = lapply(0:(2 ^ nGenes - 1), function(x)
    toBin(x, nGenes))
  
  count_ones = function(s)
    sum(unlist(gregexpr("1", s)) > 0)
  getColumn = function(num)
    rev(allStrings[sapply(allStrings, count_ones) == num])
  columns = lapply(0:nGenes, getColumn)
  
  # Construct edges
  edges = list()
  for (i in 1:(length(allStrings) - 1)) {
    for (j in (i + 1):length(allStrings)) {
      node1 = allStrings[[i]]
      node2 = allStrings[[j]]
      if (sum(strsplit(node1, NULL)[[1]] != strsplit(node2, NULL)[[1]]) == 1 &
          ((node1 %in% selectNodes &
            node2 %in% selectNodes) || is.null(selectNodes))) {
        edges = append(edges, list(c(node1, node2)))
      }
    }
  }
  
  nodes = unlist(columns)
  
  # Create data frame for nodes and positions
  layout_df = data.frame(name = nodes, stringsAsFactors = FALSE)
  
  # Assign x/y layout
  col_start = 1
  x_vals = numeric(length(nodes))
  y_vals = numeric(length(nodes))
  
  for (column_nodes in columns) {
    column_height = length(column_nodes)
    
    ys = seq((column_height - 1) / 2, -(column_height - 1) / 2)
    
    for (j in seq_along(column_nodes)) {
      idx = which(layout_df$name == column_nodes[j])
      x_vals[idx] = col_start
      y_vals[idx] = ys[j]
    }
    
    col_start = col_start + 1
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
  
  layout_df$x = x_vals
  layout_df$y = y_vals
  layout_df$fitness = unlist(lapply(layout_df$name, getFitness))
  
  # Convert edge list to data frame
  edge_df = do.call(rbind, edges)
  colnames(edge_df) = c("from", "to")
  edge_df = as.data.frame(edge_df, stringsAsFactors = FALSE)
  
  # Create tidygraph object
  g_tbl = tbl_graph(nodes = layout_df,
                    edges = edge_df,
                    directed = FALSE)
  
  # Plot with ggraph
  ggraph(g_tbl,
         layout = "manual",
         x = x,
         y = y) +
    geom_edge_link(color = "black") +
    geom_node_point(
      aes(fill = fitness),
      shape = 21,
      size = 12,
      stroke = 0.3,
      color = "black"
    ) +
    geom_node_text(aes(label = name), size = 3) +
    scale_fill_gradient(
      low = lowColor,
      high = highColor,
      na.value = "white",
      name = "Fitness"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title.position = "top"
    ) +
    ggtitle("Fitness Landscape")
}

#' Visualize CBN Model
#'
#' @param poset Poset object to visualize
#' @param nodeColor Color of nodes in resulting graph
#'
#' @return Plot (gg object) visualization of CBN model
#' @export
#'
#' @examples
#' poset = read_poset(get_examples()[1])
#' visualize_cbn_model(poset)
visualize_cbn_model <- function(poset, nodeColor = "darkgreen") {
  nodes = data.frame(name = sort(unlist(unique(as.list(
    poset[[2]]
  )))))
  edges = as.data.frame(poset[[2]])
  colnames(edges) = c("from", "to")
  
  g_tbl = tbl_graph(nodes = nodes,
                    edges = edges,
                    directed = TRUE)
  ggraph(g_tbl) +
    geom_edge_link(
      colour = "black",
      arrow = arrow(length = unit(2, 'mm')),
      end_cap = circle(3, 'mm')
    ) +
    geom_node_point(
      fill = nodeColor,
      shape = 21,
      size = 6,
      stroke = 0.3,
      color = "black"
    ) +
    geom_node_text(aes(label = name), color = "white", size = 3) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("CBN Model")
}

inverse_factorial <- function(n) {
  if (n < 1)
    return(NA)
  
  log_n <- log(n)
  log_fact <- 0
  k <- 1
  
  while (log_fact <= log_n) {
    k <- k + 1
    log_fact <- log_fact + log(k)
  }
  
  return(k - 1)
}

generate_gg_text <- function(text, bg) {
  data.frame(label = text, x = 0, y = 0) %>%
    ggplot(aes(x, y, label = label)) +
    geom_text(parse = TRUE, family = "serif") +
    theme_void() +
    theme(panel.background = element_rect(fill = bg))
}

generate_geom_node_point <- function(gra, fill, color, name) {
  node_names = gra %>%
    activate(nodes) %>%
    pull(name)
  if (length(fill) > 1) {
    color = c(NA, rep(color, length(fill)), NA)
    strokes = c(NA, rep(0.3, length(fill)), NA)
    fill = c(NA, fill, NA)
  }
  variable_cap_size(gra, rep(2, length(node_names))) +
    geom_node_point(
      fill = fill,
      shape = 21,
      size = 6,
      stroke = strokes
    ) +
    geom_node_text(aes(label = name), color = color, size = 3)
}

generate_geom_node_text <- function(gra, color, name) {
  node_names = gra %>%
    activate(nodes) %>%
    pull(name)
  if (length(color) > 1) {
    color = c(NA, color, NA)
  }
  variable_cap_size(gra, nchar(node_names)) +
    geom_node_text(aes(label = name), color = color, size = 3)
}

variable_cap_size <- function(g_tbl, node_sizes) {
  g_tbl = g_tbl %>%
    activate(edges) %>%
    mutate(cap_size = 1:gsize(g_tbl))
  
  graph = ggraph(g_tbl,
                 layout = "manual",
                 x = x,
                 y = y) + theme_void()
  
  for (i in 2:(gsize(g_tbl) - 1)) {
    filter_expr = call2("==", sym("cap_size"), i)
    
    graph = graph + geom_edge_link(
      aes(filter = !!filter_expr),
      colour = "black",
      arrow = arrow(length = unit(2, 'mm')),
      end_cap = circle(node_sizes[[i + 1]] + 2, 'mm'),
      start_cap = circle(node_sizes[[i]] + 2, 'mm')
    )
  }
  
  graph
}

#' Visualize Pathway Probabilities
#'
#' @param probabilities List or matrix of probabilities for each pathway (matrix if multiple models)
#' @param outputFile File to output to; if none provided, a plot will be returned
#' @param geneNames Gene names; if single character, rendered in circles
#' @param geneColors Gene colors
#' @param columnTitles Include column titles
#'
#' @return Plot or file name
#' @export
#'
#' @examples
#' visualize_probabilities(c(0.05, 0.03, 0.12, 0.04, 0.02, 0, 0.05, 0.04, 0.05, 0.06, 0.04, 0.02, 0.03, 0.02, 0.05, 0.03, 0.01, 0.09, 0.06, 0.04, 0, 0.08, 0.05, 0.02))
#'
#' visualize_probabilities(c(0.05, 0.03, 0.12, 0.04, 0.02, 0, 0.05, 0.04, 0.05, 0.06, 0.04, 0.02, 0.03, 0.02, 0.05, 0.03, 0.01, 0.09, 0.06, 0.04, 0, 0.08, 0.05, 0.02), geneNames = c("AAAA", "BBBB", "CCCC", "DDDD"))
#'
#' mat = matrix(c(0.1, 0.3, 0, 0.2, 0.4, 0, 0.2, 0.2, 0.1, 0, 0.2, 0.3), ncol = 2)
#' visualize_probabilities(mat, columnTitles = TRUE)
visualize_probabilities <- function(probabilities,
                                    outputFile = NULL,
                                    geneNames = as.character(1:inverse_factorial(length(probabilities))),
                                    geneColors = rainbow(length(geneNames), v = 0.5),
                                    columnTitles = TRUE) {
  numCol = 1
  if (!is.matrix(probabilities)) {
    probabilities = matrix(probabilities, ncol = 1)
  } else {
    numCol = ncol(probabilities)
  }
  
  pathway_length = inverse_factorial(nrow(probabilities))
  
  if (factorial(pathway_length) != nrow(probabilities)) {
    stop("Length of probabilities is not a factorial.")
  }
  
  perms = permutations(pathway_length, pathway_length)
  labels = sprintf("Pi[%d]", 1:length(probabilities))
  if (numCol == 1) {
    perms = perms[order(probabilities, decreasing = TRUE), ]
    # labels = labels[order(probabilities, decreasing = TRUE)]
    probabilities = matrix(sort(probabilities, decreasing = TRUE), ncol = 1)
  }
  
  generate_row <- function(row, padding = FALSE) {
    row = unlist(lapply(row, function(x) {
      geneNames[[x]]
    }))
    if (padding) {
      row = c(" " , row, " ")
      nodes = data.frame(
        name = row,
        x = c(0, 0.1:(length(row) - 2), length(row) - 2.8),
        y = rep(0, length(row))
      )
    } else {
      nodes = data.frame(
        name = row,
        x = 1:length(row),
        y = rep(0, length(row))
      )
    }
    edges = data.frame(from = row[-length(row)], to = row[-1])
    
    tbl_graph(nodes = nodes,
              edges = edges,
              directed = TRUE)
  }
  
  elements = vector("list", (factorial(pathway_length) + 1) * (2 + numCol))
  
  if (columnTitles) {
    elements[[1]] = generate_gg_text("pi", "lightgray")
    for (i in 1:numCol) {
      elements[[2+i]] = generate_gg_text("P(pi)", "lightgray")
    }
    elements[[2]] = generate_gg_text("Pathways", "lightgray")
  }
  
  for (i in 1:factorial(pathway_length)) {
    if (i %% 2 == 0) {
      bg_col = "floralwhite"
    } else {
      bg_col = "white"
    }
    
    elements[[(i) * (2 + numCol) + 1]] = generate_gg_text(labels[[i]], bg_col)
    
    gra = generate_row(perms[i, ], TRUE)
    
    if (all(unlist(lapply(as.list(geneNames), function(x) {
      nchar(x) == 1
    })))) {
      if (numCol == 1 & probabilities[[i]] == 0) {
        gra = generate_geom_node_point(gra, "white", "black", name)
      } else {
        gra = generate_geom_node_point(gra, unlist(lapply(perms[i, ], function(x) {
          geneColors[[x]]
        })), "white", name)
      }
    } else {
      if (probabilities[[i]] == 0) {
        gra = generate_geom_node_text(gra, "gray", name)
      } else {
        gra = generate_geom_node_text(gra, unlist(lapply(perms[i, ], function(x) {
          geneColors[[x]]
        })), name)
      }
      
    }
    
    elements[[(i) * (2 + numCol) + 2]] = gra +
      theme(panel.background = element_rect(fill = bg_col))
    
    for (colI in 1:numCol) {
      elements[[(i) * (2 + numCol) + 2 + colI]] = generate_gg_text(sprintf("%.2f", probabilities[i, colI]), bg_col)
    }
  }
  
  if (!columnTitles) {
    elements = elements[-(1:(2 + numCol))]
  }
  
  out = wrap_plots(elements,
                   ncol = 2 + numCol,
                   widths = c(3, (sum(nchar(gra$data$name)) + length(gra$data$name) * 4)/1.5, rep(4, numCol)))
  
  if (is.null(outputFile)) {
    plot(out)
  } else {
    ggsave(
      outputFile,
      out,
      width = 40 * (3 + (sum(nchar(gra$data$name)) + length(gra$data$name) * 4)/1.5 + 4 * numCol),
      height = 33 * length(perms),
      limitsize = FALSE,
      units = "px"
    )
    return(outputFile)
  }
}