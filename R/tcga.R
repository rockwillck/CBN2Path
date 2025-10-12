#' Get Raw TCGA Data
#'
#' @param project TCGA project ID
#'
#' @return data frame of TCGA data for given project
#' @export
#'
#' @examples
#' getRawTCGAData("TCGA-BLCA")
getRawTCGAData <- function(project) {
  query <- GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    access = "open",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  GDCdownload(query)
  maf <- GDCprepare(query)
  as.data.frame(maf)
}

#' Generate TCGA Genotype Matrix
#'
#' @param rawData Raw TCGA data generated using `getRawTCGAData`
#' @param genes Genes to generate genotype matrix on
#'
#' @return A genotype matrix where each row is a patient and each column is a gene
#' @export
#'
#' @examples
#' generateTCGAMatrix()
generateTCGAMatrix <- function(rawData=suppressMessages(getRawTCGAData("TCGA-BLCA")), genes=c("TP53","ARID1A","KDM6A","PIK3CA","RB1","EP300","FGFR3","CREBBP","STAG2","ATM")) {
  mutantRows = list()
  for (gene in genes) {
    rows = rawData[rawData$Hugo_Symbol == gene,]
    rows = rows[rows$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation"),]
    mutantRows[[gene]] = rows[!duplicated(rows$Tumor_Sample_UUID),]
  }
  
  allIDs = unique(rawData$Tumor_Sample_UUID)
  
  mat = matrix(nrow = 0, ncol = length(mutantRows))
  for (id in allIDs) {
    patientV = c()
    for (rowDF in mutantRows) {
      if (id %in% rowDF$Tumor_Sample_UUID) {
        patientV = c(patientV, 1)
      } else {
        patientV = c(patientV, 0)
      }
    }
    mat = rbind(mat, patientV)
  }
  
  rownames(mat) = NULL
  return(mat)
}
