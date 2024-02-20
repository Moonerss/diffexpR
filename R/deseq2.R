#############################################
##  do the DEG analyis use DESeq2 package   ##
##  Author: Erjie Zhao 2021.5.13           ##
#############################################


#' @name diff_deseq2_count
#' @title Do differential expression analysis on read count data using \code{DESeq2}
#' @description Do differential expression analysis using \code{DESeq2} package,
#' this is for the read count data in RNA-seq.
#' @param count_mat A matrix-like data object containing the read count data, with rows corresponding to genes and columns to samples.
#' @param label_list a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param ref.groups a character string specifying the reference group. If specified, for a given grouping variable, each of the group levels will be compared to the reference group (i.e. control group).
#' @return A dataframe with a row for the number top genes and the following columns:
#' \item{Gene}{the name of genes}
#' \item{baseMean}{}
#' \item{log2FoldChange}{log2-fold change of expression between conditions being tested.}
#' \item{lfcSE}{the standard error of the \code{log2FoldChange}}
#' \item{stat}{Wald statistic}
#' \item{pvalue}{p-values.}
#' \item{padj}{adjusted p-values by 'BH'}
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom cli cli_alert_warning cli_alert_info cli_abort
#' @author Erjie Zhao <2055469819@qq.com>
#' @examples
#' data(OSCC)
#' res <- diff_deseq2_count(count_mat = OSCC$count, label_list = OSCC$group,
#'                          ref.groups = "normal")
#' @export
diff_deseq2_count <- function(count_mat, label_list, ref.groups = NULL) {
  ## ref: https://www.plob.org/article/11506.html
  ## extract group
  condition <- label_list[match(label_list[,1], colnames(count_mat)), 2]

  if (is.character(condition)) {
    cli::cli_alert_warning('Set group information into factors')
    condition <- as.factor(condition)
  }

  condition_level <- levels(condition)
  cli::cli_alert_info('factor level: {.val {paste(condition_level, collapse = \' \')}}')


  if (!is.element(ref.groups, condition_level)) {
    cli::cli_abort(c('x' = '{.var ref.groups} must in {.var {condition_level}}'))
  }

  ## construct DESeqDataSet object
  dds_obj <- DESeqDataSetFromMatrix(countData = count_mat,
                                    colData = data.frame(condition),
                                    design = ~ condition)
  ## differential expression analysis
  vs_group <- setdiff(levels(condition), ref.groups)
  contrast_group <- c('condition', vs_group, ref.groups)
  dds <- DESeq(dds_obj, quiet = TRUE)
  res <- results(dds, contrast = contrast_group)


  ## get result
  for (i in res@elementMetadata$description) {
    cli::cli_alert_info(i)
  }
  res <- res %>% as.data.frame() %>%
    tibble::rownames_to_column(var = 'Gene')

  return(res)
}
