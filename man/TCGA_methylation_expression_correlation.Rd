% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGA_methylation_expression_correlation.R
\name{TCGA_methylation_expression_correlation}
\alias{TCGA_methylation_expression_correlation}
\title{Methylation-Expression correlation of any genes in TCGA
samples}
\usage{
TCGA_methylation_expression_correlation(
  tumor = "all",
  gene = NULL,
  nt_up = 1000,
  nt_down = 200,
  min_probe_number = 3,
  include_normal_tissues = FALSE,
  values_only = FALSE
)
}
\arguments{
\item{tumor}{\code{character} defining the TCGA tumor type. Can be one
of "SKCM", "LUAD", "LUSC", "COAD", "ESCA", "BRCA", "HNSC", or
"all" (default).}

\item{gene}{\code{character} selected gene.}

\item{nt_up}{\code{numeric(1)} indicating the number of nucleotides
upstream the TSS to define the promoter region (1000 by
default)}

\item{nt_down}{\code{numeric(1)} indicating the number of nucleotides
downstream the TSS to define the promoter region (200 by
default)}

\item{min_probe_number}{\code{numeric(1)} indicating the minimum
number of probes (with methylation values) within the selected region
to calculate its mean methylation level. Default is 3.}

\item{include_normal_tissues}{\code{logical(1)}. If \code{TRUE},
the function will include normal peritumoral tissues in addition to
tumoral samples. Default is \code{FALSE}.}

\item{values_only}{\code{logical(1)}. If \code{TRUE}, the function will return the
methylation and expression values in TCGA samples instead of the
heatmap. Default is \code{FALSE}.}
}
\value{
A scatter plot representing for each TCGA sample, gene expression
and mean methylation values of probe(s) located in its promoter region
(defined as 1000 nucleotides upstream TSS and 200 nucleotides downstream
TSS by default). If \code{values_only = TRUE}, methylation and expression values
are returned in a tibble instead.
}
\description{
Plots the correlation between methylation and expression values of
a gene in TCGA samples.
}
\details{
The coefficient of correlation is set to \code{NA} if no probes
are found in promoter regions or if less than 1\% of tumors are
positive (TPM >= 1) for the gene.
}
\examples{
\dontrun{
TCGA_methylation_expression_correlation("LUAD", gene = "TDRD1")
}
}
