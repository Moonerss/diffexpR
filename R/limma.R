#############################################
##  do the DEG analyis use limma package   ##
##  Author: Erjie Zhao 2021.5.12           ##
#############################################

#' @name diff_limma_array
#' @title Do differential expression analysis on array data using limma
#' @description Do differential expression analysis using \code{limma} package,
#' this is for the microarray data.
#' @usage diff_limma_array(array_mat, label_list, groups)
#' @param array_mat A matrix-like data object containing log-ratios or log-expression
#' values for a series of arrays, with rows corresponding to genes and columns to samples.
#' @param label_list a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param groups expressions, or character strings which can be parsed to expressions,
#' specifying contrasts, it corresponds to \code{\link[limma]{makeContrasts}}.
#' @return A dataframe with a row for the number top genes and the following columns:
#'   \item{genelist}{one or more columns of probe annotation, if genelist was included as input}
#'   \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast
#'   (for \code{topTableF} there may be several columns of log-fold-changes)}
#'   \item{CI.L}{left limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{CI.R}{right limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as Amean in the \code{MarrayLM} object}
#'   \item{t}{moderated t-statistic (omitted for \code{topTableF})}
#'   \item{F}{moderated F-statistic (omitted for \code{topTable} unless more than one coef is specified)}
#'   \item{P.Value}{raw p-value}
#'   \item{adj.P.Value}{adjusted p-value or q-value}
#'   \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}
#'
#' @import limma
#' @author Erjie Zhao <2055469819@qq.com>
#'
#' @export

diff_limma_array <- function(count_mat, label_list, groups) {

  ## extract group
  group_list <- label_list[match(label_list[,1], colnames(array_mat)), 2]

  ## Construct Design Matrices
  design <- model.matrix(~ 0 + factor(group_list))
  colnames(design) <- levels(factor(group_list))
  rownames(design) <- colnames(array_mat)

  ## Construct Matrix of Custom Contrasts
  contrast_matrix <- makeContrasts(groups, levels = design)

  ## differential analysis
  fit <- lmFit(array_mat, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  tempOutput <- topTable(fit2, coef = 1, number=Inf)
  limma_deg <- na.omit(tempOutput)

  return(limma_deg)
}

#' @name diff_limma_count
#' @title Do differential expression analysis on read count data using limma
#' @description Do differential expression analysis using \code{limma} package,
#' this is for the read count data in RNA-seq.
#' @usage diff_limma_count(expr_mat, label_list, groups, methods = c("voom", "limma-trend"))
#' @param count_mat A matrix-like data object containing the read count data, with rows corresponding to genes and columns to samples.
#' @param label_list a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param groups expressions, or character strings which can be parsed to expressions,
#' specifying contrasts, it corresponds to \code{\link[limma]{makeContrasts}}.
#' @param methods The method to the differential analysis. \code{voom} is powerful when
#' the library sizes are quite variable between samples, \code{limma-trend} the simplest
#' and most robust approach to differential exis when the sequencing depth is reasonably
#' consistent across the RNA samples.
#' @return A dataframe with a row for the number top genes and the following columns:
#'   \item{genelist}{one or more columns of probe annotation, if genelist was included as input}
#'   \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast
#'   (for \code{topTableF} there may be several columns of log-fold-changes)}
#'   \item{CI.L}{left limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{CI.R}{right limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as Amean in the \code{MarrayLM} object}
#'   \item{t}{moderated t-statistic (omitted for \code{topTableF})}
#'   \item{F}{moderated F-statistic (omitted for \code{topTable} unless more than one coef is specified)}
#'   \item{P.Value}{raw p-value}
#'   \item{adj.P.Value}{adjusted p-value or q-value}
#'   \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}
#'
#' @import limma
#' @import edgeR
#'
#' @author Erjie Zhao <2055469819@qq.com>
#'
#' @export
#'
diff_limma_count <- function(count_mat, label_list, groups, methods = c("voom", "limma-trend")) {

  ## extract group
  group_list <- label_list[match(label_list[,1], colnames(count_mat)), 2]

  ## Construct Design Matrices
  design <- model.matrix(~ 0 + factor(group_list))
  colnames(design) <- levels(factor(group_list))
  rownames(design) <- colnames(count_mat)

  ## DGEList Constructor
  dgelist <- DGEList(counts = count_mat, group = factor(group_list))

  ## Normalization and filtering
  keepgenes <- filterByExpr(dgelist, design)
  dgelist <- dgelist[keepgenes, , keep.lib.sizes = FALSE]
  dgelist <- calcNormFactors(dgelist)

  ## Differential expression analysis: limma-trend, voom
  method <- match.arg(methods)
  if (method == "limma-trend") {
    logCPM <- cpm(dgelist, log = TRUE, prior.count = 3)
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend = TRUE)
    tempoutput <- topTable(fit, coef  =ncol(design), number = Inf)
  } else {
    v <- voom(dgelist, design, plot = TRUE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    tempoutput <- topTable(fit, coef = ncol(design), number = Inf)
  }

  ## get results
  limma_deg <- na.omit(tempOutput)

  return(limma_deg)
}


#' @name diff_limma_fpkm
#' @title Do differential expression analysis on normalized RNA-seq data using limma
#' @description Do differential expression analysis using \code{limma} package,
#' this is for the normalized RNA-seq data, such as \code{TPM}, \code{RPKM}, etc.
#' @usage diff_limma_fpkm(expr_mat, label_list, groups)
#' @param expr_mat A matrix-like data object containing log-ratios or
#' log-expression values for normalized RNA-seq data, with rows corresponding
#' to genes and columns to samples.
#' @param a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param groups expressions, or character strings which can be parsed to expressions,
#' specifying contrasts, it corresponds to \code{\link[limma]{makeContrasts}}.
#'
#' @return A dataframe with a row for the number top genes and the following columns:
#'   \item{genelist}{one or more columns of probe annotation, if genelist was included as input}
#'   \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast
#'   (for \code{topTableF} there may be several columns of log-fold-changes)}
#'   \item{CI.L}{left limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{CI.R}{right limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as Amean in the \code{MarrayLM} object}
#'   \item{t}{moderated t-statistic (omitted for \code{topTableF})}
#'   \item{F}{moderated F-statistic (omitted for \code{topTable} unless more than one coef is specified)}
#'   \item{P.Value}{raw p-value}
#'   \item{adj.P.Value}{adjusted p-value or q-value}
#'   \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}
#'
#' @import limma
#' @author Erjie Zhao <2055469819@qq.com>
#'
#' @export
#'
#' reference：https://support.bioconductor.org/p/56275/#56299
diff_limma_fpkm <- function(expr_mat, label_list, groups) {
  ## check log2
  logs <- check_log2(expr_mat)

  ## log2 tranform
  if (logs) {
    message("The data didn't log2 tranform, do it ...")
    expr_mat[which(expr_mat <= 0)] <- 0
    expr_mat <- log2(expr_mat + 1)
  }

  ## extract group
  group_list <- label_list[match(label_list[,1], colnames(expr_mat)), 2]

  ## Construct Design Matrices
  design <- model.matrix(~ 0 + factor(group_list))
  colnames(design) <- levels(factor(group_list))
  rownames(design) <- colnames(expr_mat)

  ## Construct Matrix of Custom Contrasts
  contrast_matrix <- makeContrasts(groups, levels = design)

  ## differential analysis
  fit <- lmFit(expr_mat, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
  tempOutput <- topTable(fit2, coef = 1, number=Inf)
  limma_deg <- na.omit(tempOutput)

  return(limma_deg)
}



# check whether done log2 transformation
# logical value, TRUE not log2, FALSE, log2ed

check_log2 <- function(mat) {
  qx <- as.numeric(quantile(mat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  loged <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return(loged)
}


