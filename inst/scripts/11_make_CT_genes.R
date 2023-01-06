## Code to prepare `CT_genes` dataset goes here

library("tidyverse")
library("SummarizedExperiment")
library("biomaRt")
library("CTexploreR")

load("../extdata/CT_list.rda")
load(file = "../../data/TCGA_CT_methylation.rda")
load(file = "../../data/CT_mean_methylation_in_tissues.rda")
load(file = "../../data/CT_methylation_in_tissues.rda")

################################################################################
## Add CpG densities and promoter methylation analysis in normal tissues
################################################################################
CT_genes <- CT_genes %>%
  left_join(as_tibble(rowData(CT_mean_methylation_in_tissues)))

################################################################################
## Based on DAC induction and on methylation levels in normal tissues
## (if available).For some genes, methylation analysis was not possible due
## to multimapping issues. In this case, genes are still considered as regulated
## by methylation if they show a strong activation in cells treated with
## 5-Aza-2′-Deoxycytidine.
################################################################################

CT_genes <- CT_genes %>%
  mutate(regulation = case_when(
    DAC == "induced" &
      (is.na(methylation_in_tissues) |
         methylation_in_tissues == "methylated_in_somatic_unmethylated_in_germline")
    ~ "methylation"))

CT_genes[is.na(CT_genes$regulation), "regulation"] <- "not_methylation"

################################################################################
## Add correlation value between methylation and expression from TCGA data
################################################################################

met_exp_corr_TCGA <- sapply(CT_genes$external_gene_name,
                            TCGA_methylation_expression_correlation,
                            expression_database = TCGA_TPM,
                            methylation_database = TCGA_CT_methylation,
                            tumor = c("SKCM", "LUAD", "LUSC", "COAD", "ESCA"),
                            corr_coeff = TRUE)
met_exp_corr_TCGA <- enframe(met_exp_corr_TCGA) %>%
  dplyr::rename(external_gene_name = name, met_exp_corr_TCGA = value)
CT_genes <- CT_genes %>%
  left_join(met_exp_corr_TCGA)

CT_genes <- CT_genes %>%
  mutate(X_linked = ifelse(chromosome_name == "X", "chrX", "not_chrX"))

################################################################################
## Flag genes as "oncogenic" or "tumor suppressor" using
## [Cancermine](http://bionlp.bcgsc.ca/cancermine/), a literature-mined database
## of drivers, oncogenes and tumor suppressors in cancer.
################################################################################
if (!file.exists("../../../CTdata/inst/extdata/cancermine_collated.tsv")) {
  download.file(url = "http://bionlp.bcgsc.ca/cancermine/session/d64941b57313d6774447c97940158a37/download/gene_download_collated_all?w=",
                destfile = "../../../CTdata/inst/extdata/cancermine_collated.tsv")
}
cancermine <- read_tsv("../../../CTdata/inst/extdata/cancermine_collated.tsv")
oncogenes <- cancermine %>%
  filter(role == "Oncogene") %>%
  pull(gene_normalized)
TS <- cancermine %>%
  filter(role == "Tumor_Suppressor") %>%
  pull(gene_normalized)
CT_genes <- CT_genes %>%
  mutate(oncogene = case_when(external_gene_name %in% oncogenes ~ "oncogene")) %>%
  mutate(tumor_suppressor = case_when(external_gene_name %in% TS ~ "tumor_suppressor"))

################################################################################
# Reorder the final table
################################################################################
CT_genes <- CT_genes %>%
  dplyr::select("ensembl_gene_id", "external_gene_name", "family", "X_linked",
                "GTEX_category", "TPM_testis", "max_TPM_somatic",
                "q75_TPM_somatic", "multimapping_analysis", "testis_specificity",
                "percent_of_positive_CCLE_cell_lines",
                "percent_of_negative_CCLE_cell_lines", "max_TPM_in_CCLE",
                "CCLE_category", "percent_pos_tum", "percent_neg_tum",
                "max_TPM_in_TCGA", "TCGA_category","DAC", "methylation_in_tissues",
                "regulation", "met_exp_corr_TCGA", "CpG_density","CpG_promoter",
                "external_transcript_name", "ensembl_transcript_id",
                "chromosome_name", "strand", "transcription_start_site",
                "transcript_length", "transcript_biotype", "oncogene",
                "tumor_suppressor")

usethis::use_data(CT_genes, overwrite = TRUE)



