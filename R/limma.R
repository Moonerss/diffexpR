#############################################
##  do the DEG analyis use limma package   ##
##  Author: Erjie Zhao 2021.5.12           ##
#############################################

#' @name diff_limma_array
#' @title Do differential expression analysis on array data using limma
#' @description Do differential expression analysis using \code{limma} package,
#' this is for the microarray data.
#' @param array_mat A matrix-like data object containing log-ratios or log-expression
#' values for a series of arrays, with rows corresponding to genes and columns to samples.
#' @param label_list a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param ref.groups a character string specifying the reference group. If specified, for a given grouping variable, each of the group levels will be compared to the reference group (i.e. control group).
#' @return A dataframe with a row for the number top genes and the following columns:
#'   \item{Gene}{the name of genes}
#'   \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast}
#'   \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as Amean in the \code{MarrayLM} object}
#'   \item{t}{moderated t-statistic (omitted for \code{topTableF})}
#'   \item{P.Value}{raw p-value}
#'   \item{adj.P.Value}{adjusted p-value or q-value}
#'   \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}
#' @import limma
#' @importFrom stats na.omit model.matrix
#' @importFrom cli cli_alert_warning cli_alert_info
#' @author Erjie Zhao <2055469819@qq.com>
#' @examples
#' data(eset)
#' res <- diff_limma_array(array_mat = eset$array, label_list = eset$group, ref.groups = "liver")
#' @export

diff_limma_array <- function(array_mat, label_list, ref.groups = NULL) {

  ## extract group
  group_list <- label_list[match(label_list[,1], colnames(array_mat)), 2]

  if (is.character(group_list)) {
    cli::cli_alert_warning('Set group information into factors')
    group_list <- as.factor(group_list)
  }

  condition_level <- levels(group_list)
  cli::cli_alert_info('factor level: {.val {paste(condition_level, collapse = \' \')}}')

  if (!is.element(ref.groups, condition_level)) {
    cli::cli_abort(c('x' = '{.var ref.groups} must in {.var {condition_level}}'))
  }
  ## Construct Design Matrices
  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(array_mat)

  ## Construct Matrix of Custom Contrasts
  vs_group <- setdiff(condition_level, ref.groups)
  groups <- paste(vs_group, ref.groups, sep = '-')
  cli::cli_alert_info('compare condition: {.val {vs_group}} vs {.val {ref.groups}}')
  contrast_matrix <- makeContrasts(contrasts = groups, levels = design)

  ## differential analysis
  fit <- lmFit(array_mat, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  tempoutput <- topTable(fit2, coef = 1, number = Inf)
  limma_deg <- na.omit(tempoutput)

  limma_deg <- limma_deg %>%
    tibble::rownames_to_column(var = 'Gene')

  return(limma_deg)
}

#' @name diff_limma_count
#' @title Do differential expression analysis on read count data using limma
#' @description Do differential expression analysis using \code{limma} package,
#' this is for the read count data in RNA-seq.
#' @usage diff_limma_count(count_mat, label_list, groups, methods = c("voom", "limma-trend"))
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
#'   \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast}
#'   \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as Amean in the \code{MarrayLM} object}
#'   \item{t}{moderated t-statistic (omitted for \code{topTableF})}
#'   \item{P.Value}{raw p-value}
#'   \item{adj.P.Value}{adjusted p-value or q-value}
#'   \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}
#'
#' @import limma
#' @import edgeR
#' @importFrom stats na.omit model.matrix
#'
#' @author Erjie Zhao <2055469819@qq.com>
#'
#' @examples
#' data(OSCC)
#' res <- diff_limma_count(count_mat = OSCC$count, label_list = OSCC$group,
#'                         groups = c("normal-tumor"), methods = "voom")
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

  ## Construct Matrix of Custom Contrasts
  contrast_matrix <- makeContrasts(contrasts = groups, levels = design)

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
    fit <- contrasts.fit(fit, contrast_matrix)
    fit <- eBayes(fit, trend = TRUE)
    tempoutput <- topTable(fit, coef = 1, number = Inf)
  } else if (method == "voom") {
    v <- voom(dgelist, design)
    fit <- lmFit(v, design)
    fit <- contrasts.fit(fit, contrast_matrix)
    fit <- eBayes(fit)
    tempoutput <- topTable(fit, coef = 1, number = Inf)
  }

  ## get results
  limma_deg <- na.omit(tempoutput)

  return(limma_deg)
}


#' @name diff_limma_normalize
#' @title Do differential expression analysis on normalized RNA-seq data using limma
#' @description Do differential expression analysis using \code{limma} package,
#' this is for the normalized RNA-seq data, such as \code{FPKM}, \code{TPM}, \code{RPKM}, etc.
#' @usage diff_limma_normalize(expr_mat, label_list, groups)
#' @param expr_mat A matrix-like data object containing log-ratios or
#' log-expression values for normalized RNA-seq data, with rows corresponding
#' to genes and columns to samples.
#' @param label_list a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param groups expressions, or character strings which can be parsed to expressions,
#' specifying contrasts, it corresponds to \code{\link[limma]{makeContrasts}}.
#'
#' @return A dataframe with a row for the number top genes and the following columns:
#'   \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast}
#'   \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as Amean in the \code{MarrayLM} object}
#'   \item{t}{moderated t-statistic (omitted for \code{topTableF})}
#'   \item{P.Value}{raw p-value}
#'   \item{adj.P.Value}{adjusted p-value or q-value}
#'   \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}
#'
#' @import limma
#' @importFrom stats na.omit quantile model.matrix
#' @author Erjie Zhao <2055469819@qq.com>
#'
#' @examples
#' data(OSCC)
#' res <- diff_limma_normalize(expr_mat = OSCC$rpkm, label_list = OSCC$group,
#'                             groups = c("normal-tumor"))
#'
#' @export
#'
diff_limma_normalize <- function(expr_mat, label_list, groups) {
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
  contrast_matrix <- makeContrasts(contrasts = groups, levels = design)

  ## differential analysis
  fit <- lmFit(expr_mat, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)

  ## the reference：https://support.bioconductor.org/p/56275/#56299
  fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
  tempoutput <- topTable(fit2, coef = 1, number=Inf)
  limma_deg <- na.omit(tempoutput)

  return(limma_deg)
}



# check whether done log2 transformation
# logical value, TRUE not log2, FALSE, log2ed

check_log2 <- function(mat) {
  qx <- as.numeric(quantile(mat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  loged <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return(loged)
}


