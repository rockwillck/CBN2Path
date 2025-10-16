toBin <- function(n, nG) {
  x <- paste(as.integer(rev(intToBits(n))), collapse = "")
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
#' genotypes <- c(
#'     "0000",
#'     "1000",
#'     "0100",
#'     "0010",
#'     "0001",
#'     "1100",
#'     "1010",
#'     "1001",
#'     "0110",
#'     "0101",
#'     "0011",
#'     "1110",
#'     "1101",
#'     "1011",
#'     "0111",
#'     "1111"
#' )
#' #
#' colIntensity <- c(0, rep(0.25, 4), rep(0.5, 6), rep(0.75, 4), 1)
#' visualizeFitnessLandscape(colIntensity)
visualizeFitnessLandscape <- function(fitness,
                                      selectNodes = NULL,
                                      nGenes = 4,
                                      lowColor = "white",
                                      highColor = "blue") {
  allStrings <- lapply(0:(2^nGenes - 1), function(x) {
    toBin(x, nGenes)
  })

  countOnes <- function(s) {
    sum(unlist(gregexpr("1", s)) > 0)
  }
  getColumn <- function(num) {
    rev(allStrings[sapply(allStrings, countOnes) == num])
  }
  columns <- lapply(0:nGenes, getColumn)

  # Construct edges
  edges <- list()
  for (i in 1:(length(allStrings) - 1)) {
    for (j in (i + 1):length(allStrings)) {
      node1 <- allStrings[[i]]
      node2 <- allStrings[[j]]
      if (sum(strsplit(node1, NULL)[[1]] != strsplit(node2, NULL)[[1]]) == 1 &
          ((node1 %in% selectNodes &
            node2 %in% selectNodes) || is.null(selectNodes))) {
        edges <- append(edges, list(c(node1, node2)))
      }
    }
  }

  nodes <- unlist(columns)

  # Create data frame for nodes and positions
  layoutDf <- data.frame(name = nodes, stringsAsFactors = FALSE)

  # Assign x/y layout
  colStart <- 1
  xVals <- numeric(length(nodes))
  yVals <- numeric(length(nodes))

  for (columnNodes in columns) {
    columnHeight <- length(columnNodes)

    ys <- seq((columnHeight - 1) / 2, -(columnHeight - 1) / 2)

    for (j in seq_along(columnNodes)) {
      idx <- which(layoutDf$name == columnNodes[j])
      xVals[idx] <- colStart
      yVals[idx] <- ys[j]
    }

    colStart <- colStart + 1
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

  layoutDf$x <- xVals
  layoutDf$y <- yVals
  layoutDf$fitness <- unlist(lapply(layoutDf$name, getFitness))

  # Convert edge list to data frame
  edgeDf <- do.call(rbind, edges)
  colnames(edgeDf) <- c("from", "to")
  edgeDf <- as.data.frame(edgeDf, stringsAsFactors = FALSE)

  # Create tidygraph object
  gTbl <- tbl_graph(
    nodes = layoutDf,
    edges = edgeDf,
    directed = FALSE
  )

  # Plot with ggraph
  ggraph(gTbl,
         layout = "manual",
         x = x,
         y = y
  ) +
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
#' @param numNodes Number of nodes (default is the larger number between 4 and the largest index given in the poset)
#'
#' @return Plot (gg object) visualization of CBN model
#' @export
#'
#' @examples
#' poset <- readPoset(getExamples()[1])
#' visualizeCBNModel(poset$sets)
visualizeCBNModel <- function(poset, nodeColor = "darkgreen", numNodes=max(4, max(poset))) {
  if (dim(poset)[2]<2) {
    numNodes = 4
  } else {
    edges <- as.data.frame(poset)
    colnames(edges) <- c("from", "to")
  }
  nodes <- data.frame(name = 1:numNodes)

  if (dim(poset)[2]<2) {
    gTbl <- tbl_graph(
      nodes = nodes,
      directed = TRUE
    )
  } else {
    gTbl <- tbl_graph(
      nodes = nodes,
      edges = edges,
      directed = TRUE
    )
  }

  ggraph(gTbl) +
    geom_edge_link(
      colour = "black",
      arrow = arrow(length = unit(16, "pt")),
      end_cap = circle(12, "pt")
    ) +
    geom_node_point(
      fill = nodeColor,
      shape = 21,
      size = 12,
      stroke = 0.3,
      color = "black"
    ) +
    geom_node_text(aes(label = name), color = "white", size = 5) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("CBN Model")
}

inverseFactorial <- function(n) {
  if (n < 1) {
    return(NA)
  }

  logN <- log(n)
  logFact <- 0
  k <- 1

  while (logFact <= logN) {
    k <- k + 1
    logFact <- logFact + log(k)
  }

  return(k - 1)
}

generateGgText <- function(text, bg, color = "black") {
  data.frame(label = text, x = 0, y = 0) %>%
    ggplot(aes(x, y, label = label)) +
    geom_text(parse = TRUE, family = "serif", color = color) +
    theme_void() +
    theme(panel.background = element_rect(fill = bg))
}

generateGeomNodePoint <- function(gra, fill, color, name, arrowColor = "black") {
  nodeNames <- gra %>%
    activate(nodes) %>%
    pull(name)
  if (length(fill) > 1) {
    color <- c(NA, rep(color, length(fill)), NA)
    strokes <- c(NA, rep(0.3, length(fill)), NA)
    fill <- c(NA, fill, NA)
  } else {
    strokes <- c(NA, rep(0.3, length(nodeNames) - 2), NA)
    fill <- c(NA, rep(fill, length(nodeNames) - 2), NA)
  }
  variableCapSize(gra, rep(4, length(nodeNames)), arrowColor) +
    geom_node_point(
      fill = fill,
      shape = 21,
      size = 6,
      stroke = strokes
    ) +
    geom_node_text(aes(label = name), color = color, size = 3)
}

generateGeomNodeText <- function(gra, color, name, arrowColor) {
  nodeNames <- gra %>%
    activate(nodes) %>%
    pull(name)
  if (length(color) > 1) {
    color <- c(NA, color, NA)
  }
  variableCapSize(gra, rep(max(nchar(nodeNames)) * 3, length(nodeNames)), arrowColor) +
    geom_node_text(aes(label = name), color = color, size = 3)
}

variableCapSize <- function(gTbl, nodeSizes, arrowColor) {
  gTbl <- gTbl %>%
    activate(edges) %>%
    mutate(cap_size = 1:gsize(gTbl))

  graph <- ggraph(gTbl,
                  layout = "manual",
                  x = x,
                  y = y
  ) + theme_void()

  for (i in 2:(gsize(gTbl) - 1)) {
    filterExpr <- call2("==", sym("cap_size"), i)

    graph <- graph + geom_edge_link(
      aes(filter = !!filterExpr),
      colour = arrowColor,
      arrow = arrow(length = unit(5, "pt")),
      end_cap = circle(nodeSizes[[i + 1]] + 5, "pt"),
      start_cap = circle(nodeSizes[[i]] + 5, "pt")
    )
  }

  graph
}

ptToMm <- function(pts) {
  pts / 2.83465
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
#' visualizeProbabilities(c(0.05, 0.03, 0.12, 0.04, 0.02, 0, 0.05, 0.04, 0.05, 0.06, 0.04, 0.02, 0.03, 0.02, 0.05, 0.03, 0.01, 0.09, 0.06, 0.04, 0, 0.08, 0.05, 0.02))
#'
#' visualizeProbabilities(c(0.05, 0.03, 0.12, 0.04, 0.02, 0, 0.05, 0.04, 0.05, 0.06, 0.04, 0.02, 0.03, 0.02, 0.05, 0.03, 0.01, 0.09, 0.06, 0.04, 0, 0.08, 0.05, 0.02), geneNames = c("AAAA", "BBBB", "CCCC", "DDDD"))
#'
#' mat <- matrix(c(0.1, 0.3, 0, 0.2, 0.4, 0, 0.2, 0.2, 0.1, 0, 0.2, 0.3), ncol = 2)
#' visualizeProbabilities(mat, columnTitles = TRUE)
visualizeProbabilities <- function(probabilities,
                                   outputFile = NULL,
                                   geneNames = as.character(1:inverseFactorial(length(probabilities))),
                                   geneColors = rainbow(length(geneNames), v = 0.5),
                                   columnTitles = TRUE) {
  numCol <- 1
  if (!is.matrix(probabilities)) {
    probabilities <- matrix(probabilities, ncol = 1)
  } else {
    numCol <- ncol(probabilities)
  }

  pathwayLength <- inverseFactorial(nrow(probabilities))

  if (factorial(pathwayLength) != nrow(probabilities)) {
    stop("Length of probabilities is not a factorial.")
  }

  perms <- permutations(pathwayLength, pathwayLength)
  labels <- sprintf("Pi[%d]", 1:length(probabilities))
  if (numCol == 1) {
    perms <- perms[order(probabilities, decreasing = TRUE), ]
    # labels = labels[order(probabilities, decreasing = TRUE)]
    probabilities <- matrix(sort(probabilities, decreasing = TRUE), ncol = 1)
  }

  generateRow <- function(row, padding = FALSE) {
    row <- unlist(lapply(row, function(x) {
      geneNames[[x]]
    }))
    widthMiddle <- length(row)
    if (padding) {
      row <- c(" ", row, " ")
      xCoord <- c(0, 15)
      for (i in 1:(length(row) - 3)) {
        currentX <- xCoord[[i + 1]]
        diff <- max(nchar(row)) * 3 + 10 + 20
        xCoord <- c(xCoord, currentX + diff)
      }
      xCoord <- c(xCoord, xCoord[[length(row) - 1]] + 15)
      widthMiddle <- xCoord[[length(xCoord) - 1]]
      nodes <- data.frame(
        name = row,
        x = xCoord,
        y = rep(0, length(row))
      )
    } else {
      nodes <- data.frame(
        name = row,
        x = 1:length(row),
        y = rep(0, length(row))
      )
    }
    edges <- data.frame(from = row[-length(row)], to = row[-1])

    list(tbl_graph(
      nodes = nodes,
      edges = edges,
      directed = TRUE
    ), widthMiddle)
  }

  elements <- vector("list", (factorial(pathwayLength) + 1) * (2 + numCol))

  if (columnTitles) {
    elements[[1]] <- generateGgText("pi", "lightgray")
    for (i in 1:numCol) {
      elements[[2 + i]] <- generateGgText("P(pi)", "lightgray")
    }
    elements[[2]] <- generateGgText("Pathways", "lightgray")
  }

  for (i in 1:factorial(pathwayLength)) {
    if (i %% 2 == 0) {
      bgCol <- "floralwhite"
    } else {
      bgCol <- "white"
    }

    textColor <- "black"
    if (probabilities[[i]] == 0 & numCol == 1) {
      textColor <- "lightgray"
    }
    elements[[(i) * (2 + numCol) + 1]] <- generateGgText(labels[[i]], bgCol, textColor)

    genRow <- generateRow(perms[i, ], TRUE)
    gra <- genRow[[1]]
    widthMiddle <- genRow[[2]]

    if (all(unlist(lapply(as.list(geneNames), function(x) {
      nchar(x) == 1
    })))) {
      if (numCol == 1 & probabilities[[i]] == 0) {
        gra <- generateGeomNodePoint(gra, "white", "black", name, textColor)
      } else {
        gra <- generateGeomNodePoint(gra, unlist(lapply(perms[i, ], function(x) {
          geneColors[[x]]
        })), "white", name, textColor)
      }
    } else {
      if (probabilities[[i]] == 0) {
        gra <- generateGeomNodeText(gra, "gray", name, textColor)
      } else {
        gra <- generateGeomNodeText(gra, unlist(lapply(perms[i, ], function(x) {
          geneColors[[x]]
        })), name, textColor)
      }
    }

    elements[[(i) * (2 + numCol) + 2]] <- gra +
      theme(panel.background = element_rect(fill = bgCol))

    for (colI in 1:numCol) {
      elements[[(i) * (2 + numCol) + 2 + colI]] <- generateGgText(sprintf("%.2f", probabilities[i, colI]), bgCol, textColor)
    }
  }

  if (!columnTitles) {
    elements <- elements[-(1:(2 + numCol))]
  }

  out <- wrap_plots(elements,
                    ncol = 2 + numCol,
                    widths = c(24, widthMiddle * 2, rep(32, numCol))
  )

  if (is.null(outputFile)) {
    plot(out)
  } else {
    ggsave(
      outputFile,
      out,
      width = ptToMm(24 + widthMiddle * 2 + 32 * numCol),
      height = 2 * length(perms),
      limitsize = FALSE,
      units = "mm"
    )
    return(outputFile)
  }
}
