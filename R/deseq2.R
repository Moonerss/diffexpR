#############################################
##  do the DEG analyis use DESeq2 package   ##
##  Author: Erjie Zhao 2021.5.13           ##
#############################################


#' @name diff_deseq2_count
#' @title Do differential expression analysis on read count data using \code{DESeq2}
#' @description Do differential expression analysis using \code{DESeq2} package,
#' this is for the read count data in RNA-seq.
#' @usage diff_deseq2_count(count_mat, label_list, groups)
#' @param count_mat A matrix-like data object containing the read count data, with rows corresponding to genes and columns to samples.
#' @param label_list a dataframe contain samples in first column and groups information
#' in second column in the form of factors.
#' @param groups the same as \code{contrast} argument in \code{\link[DESeq2]{results}}, as the name of a factor in the design formula is \code{condition}
#' @return A dataframe with a row for the number top genes and the following columns:
#' \item{baseMean}{}
#' \item{log2FoldChange}{log2-fold change of expression between conditions being tested.}
#' \item{lfcSE}{the standard error of the \code{log2FoldChange}}
#' \item{stat}{Wald statistic}
#' \item{pvalue}{p-values.}
#' \item{padj}{adjusted p-values by 'BH'}
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @author Erjie Zhao <2055469819@qq.com>
#' @examples
#' data(OSCC)
#' res <- diff_deseq2_count(count_mat = OSCC$count, label_list = OSCC$group,
#'                          groups = c("condition", "tumor", "normal"))
#' @export
diff_deseq2_count <- function(count_mat, label_list, groups) {
  ## ref: https://www.plob.org/article/11506.html

  ## extract group
  condition <- label_list[match(label_list[,1], colnames(count_mat)), 2]

  ## construct DESeqDataSet object
  dds_obj <- DESeqDataSetFromMatrix(countData = count_mat,
                                    colData = data.frame(condition),
                                    design = ~ condition)
  ## differential expression analysis
  dds <- DESeq(dds_obj)
  res <- results(dds, contrast = groups)

  ## get result
  res <- as.data.frame(res)
  return(res)
}
