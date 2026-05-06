#' @importFrom Matrix colSums
#' @importFrom stats lm predict
scQC <- function(
    X,
    mtThreshold = 0.1,
    minLSize = 1000,
    minCells = 25,
    qc = TRUE
) {
  if (inherits(X, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Please install 'SeuratObject' package.")
    }
    assay <- SeuratObject::DefaultAssay(X)
    X <- SeuratObject::GetAssayData(X[[assay]], "counts")
  }
  if (isFALSE(qc)) {
    return(X)
  }
  librarySize <- Matrix::colSums(X)
  X <- X[, librarySize >= minLSize, drop = FALSE]

  librarySize <- Matrix::colSums(X)
  nGenes <- Matrix::colSums(X > 0)

  genesLM <- lm(nGenes ~ librarySize)
  genesLM <- as.data.frame(predict(
    genesLM,
    data.frame(librarySize),
    interval = "prediction"
  ))

  mtGenes <- grep("^MT-", toupper(rownames(X)))
  if (isFALSE(length(mtGenes) > 0)) {
    selectedCells <- (nGenes > genesLM$lwr) &
      (nGenes < genesLM$upr) & (librarySize < 2 * mean(librarySize))
  } else {
    mtCounts <- Matrix::colSums(X[mtGenes, , drop = FALSE])
    mtProportion <- mtCounts / librarySize
    mtLM <- lm(mtCounts ~ librarySize)
    mtLM <- as.data.frame(predict(
      mtLM,
      data.frame(librarySize),
      interval = "prediction"
    ))
    selectedCells <- (mtCounts > mtLM$lwr) &
      (mtCounts < mtLM$upr) &
      (nGenes > genesLM$lwr) &
      (nGenes < genesLM$upr) &
      (mtProportion <= mtThreshold) &
      (librarySize < 2 * mean(librarySize))
  }
  X <- X[, selectedCells, drop = FALSE]
  X <- X[Matrix::rowSums(X > 0) >= minCells, , drop = FALSE]
  X
}

#' @importFrom grDevices boxplot.stats
#' @importFrom Matrix colSums rowMeans
#' @importFrom methods as
#' @importClassesFrom Matrix CsparseMatrix
scQC2 <- function(
    X,
    minLibSize = 1000,
    removeOutlierCells = TRUE,
    minPCT = 0.05,
    maxMTratio = 0.1
) {
  nCellsInit <- ncol(X)
  nGenesInit <- nrow(X)

  # Remove negative values
  X[X < 0] <- 0

  # Filter by minimum library size
  lSize <- Matrix::colSums(X)
  X <- X[, lSize > minLibSize, drop = FALSE]

  # Remove outlier cells
  if (removeOutlierCells) {
    lSize <- Matrix::colSums(X)
    X <- X[, !lSize %in% boxplot.stats(lSize)$out, drop = FALSE]
  }

  # Filter by mitochondrial ratio
  mtGenes <- grepl('^MT-', toupper(rownames(X)), ignore.case = TRUE)
  if (sum(mtGenes) > 0) {
    mtRate <- Matrix::colSums(X[mtGenes, ]) / Matrix::colSums(X)
    X <- X[, mtRate < maxMTratio, drop = FALSE]
  }

  # Filter by minimum expression percentage
  X <- X[Matrix::rowMeans(X != 0) > minPCT, , drop = FALSE]

  # Convert to sparse matrix
  gNames <- rownames(X)
  cNames <- colnames(X)
  X <- as(X, 'CsparseMatrix')
  rownames(X) <- gNames
  colnames(X) <- cNames
  X
}


