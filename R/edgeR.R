#############################################
##  do the DEG analyis use edgeR package   ##
##  Author: Erjie Zhao 2021.5.13           ##
#############################################

#' @name diff_edger_count
#' @title Do differential expression analysis on read count data using \code{edgeR}
#' @description Do differential expression analysis using \code{edgeR} package,
#' this is for the read count data in RNA-seq.
#' @usage diff_edger_count(count_mat, label_list, groups)
#' @param count_mat A matrix-like data object containing the read count data, with rows corresponding to genes and columns to samples.
#' @param label_list a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param groups expressions, or character strings which can be parsed to expressions,
#' specifying contrasts, it corresponds to \code{\link[limma]{makeContrasts}}.
#' @return A dataframe with a row for the number top genes and the following columns:
#' \item{logFC}{log2-fold change of expression between conditions being tested.}
#' \item{logCPM}{average log2-counts per million.}
#' \item{LR}{likelihood ratio statistics.}
#' \item{PValue}{p-values.}
#' \item{FDR}{false discovery rate.}
#'
#' @import edgeR
#' @author Erjie Zhao <2055469819@qq.com>
#' @export
#' @examples
#' data(OSCC)
#' res <- diff_edger_count(count_mat = OSCC$count, label_list = OSCC$group, groups = c('normal-tumor'))
#'

diff_edger_count <- function(count_mat, label_list, groups) {
  ## reference: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

  ## extract group
  group_list <- label_list[match(label_list[,1], colnames(count_mat)), 2]

  ## Construct Design Matrices
  design <- model.matrix(~ 0 + factor(group_list))
  colnames(design) <- levels(factor(group_list))
  rownames(design) <- colnames(count_mat)

  ## Construct Matrix of Custom Contrasts
  contrast_matrix <- makeContrasts(contrasts = groups, levels = design)

  ## DGEList Constructor
  dgelist <- DGEList(counts = count_mat, group = factor(group_list))
  dgelist <- estimateDisp(dgelist, design, robust = TRUE)
  fit <- glmFit(dgelist, design)
  lrt <- glmLRT(fit, contrast = contrast_matrix)
  results_edgeR <- topTags(lrt, n = nrow(count_mat), sort.by = "none")

  ## extract result
  res <- results_edgeR@.Data[[1]]
  return(res)
}
