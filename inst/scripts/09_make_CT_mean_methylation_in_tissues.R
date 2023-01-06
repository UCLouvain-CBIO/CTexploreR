## Code to prepare `CT_mean_methylation_in_tissues` dataset goes here

library(GenomicRanges)
library(tidyverse)
library(SummarizedExperiment)

load("../extdata/CT_list.rda")
load("../../data/CT_methylation_in_tissues.rda")

## Promoter region is defined as `nt_up` nucleotides upstream TSS
## and `nt_down` nucleotides downstream TSS
nt_up <- 1000
nt_down <- 200

## Calculate mean methylation of each promoter in tissues
## and store CpG number by promoter
prom_mean_met_in_tissues <- tibble(tissue =
                                     c(colnames(CT_methylation_in_tissues),
                                       "CpG_number"))

for (gene in CT_list$external_gene_name) {

  TSS <- CT_list %>%
    filter(external_gene_name == gene) %>%
    pull(transcription_start_site)

  chr <- CT_list %>%
    filter(external_gene_name == gene) %>%
    pull(chromosome_name)

  strand <- CT_list %>%
    filter(external_gene_name == gene) %>%
    pull(strand)

  TSS <- CT_list %>%
    filter(external_gene_name == gene) %>%
    pull(transcription_start_site)

  if (strand == 1) { # Analyse region at +/- nt_up and nt_down around TSS
    promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                           strand = '+',
                           ranges = IRanges(start = TSS - nt_up,
                                            end = TSS + nt_down))
    promoter_gr$TSS <- TSS
  }

  if (strand == -1) { # Analyse region at +/- nt_up and nt_down around TSS
    promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                           strand = '-',
                           ranges = IRanges(start = TSS - nt_down,
                                            end = TSS + nt_up))
    promoter_gr$TSS <- TSS
  }

  promoter_methylation <- subsetByOverlaps(CT_methylation_in_tissues,
                                           promoter_gr)
  tmp <- enframe(colMeans(assay(promoter_methylation), na.rm = TRUE),
                 name = "tissue", value = gene)
  # Store the number of CpG in the promoter region
  tmp <- rbind(tmp, c("CpG_number", dim(promoter_methylation)[1]))
  prom_mean_met_in_tissues <- left_join(prom_mean_met_in_tissues, tmp)
}

## Store CpG number by promoter
CT_CpG_number <- prom_mean_met_in_tissues %>%
  filter(tissue == "CpG_number") %>%
  pivot_longer(names_to = "external_gene_name", values_to = "CpG_number",
               -tissue) %>%
  dplyr::select(-tissue) %>%
  mutate(CpG_number = as.integer(CpG_number))

## Store mean methylation level by tissue by promoter
prom_mean_met_in_tissues <- prom_mean_met_in_tissues %>%
  filter(tissue != "CpG_number") %>%
  pivot_longer(names_to = "external_gene_name", values_to = "mean_methylation",
               -tissue) %>%
  mutate(mean_methylation = as.numeric(mean_methylation)) %>%
  pivot_wider(names_from = tissue, values_from = mean_methylation)
mat <- prom_mean_met_in_tissues %>%
  dplyr::select(-external_gene_name) %>%
  as.matrix()
rownames(mat) <- prom_mean_met_in_tissues$external_gene_name

## Calculate CpG densities and ratios of methylation in somatic tissues vs sperm
methylation_analysis <- tibble(
  external_gene_name = prom_mean_met_in_tissues$external_gene_name,
  somatic_met = prom_mean_met_in_tissues %>%
    dplyr::select(-c(external_gene_name, placenta, testis, sperm)) %>%
    rowwise %>%
    rowMeans(na.rm = TRUE),
  sperm_met = prom_mean_met_in_tissues %>%
    dplyr::select(sperm) %>%
    pull(sperm))

methylation_analysis <- methylation_analysis %>%
  left_join(CT_list %>%
              dplyr::select(external_gene_name, ensembl_gene_id)) %>%
  mutate(ratio_somatic_sperm = somatic_met / sperm_met) %>%
  left_join(CT_CpG_number) %>%
  mutate(CpG_density = CpG_number / (nt_up + nt_down) * 100) %>%
  mutate(CpG_promoter = case_when(CpG_density < 2 ~ "CpG_low",
                                 CpG_density >= 2 &
                                   CpG_density < 4 ~ "CpG_intermediate",
                                 CpG_density >= 4 ~ "CpG_high")) %>%
  mutate(methylation_in_tissues =
           case_when(somatic_met < 50 ~ "unmethylated_in_somatic",
                     somatic_met >= 50 &
                       ratio_somatic_sperm > 2 ~ "methylated_in_somatic_unmethylated_in_germline",
                     somatic_met >= 50 &
                       ratio_somatic_sperm <= 2 ~ "methylated_in_somatic_and_germline")) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, CpG_density,
                CpG_promoter, methylation_in_tissues)


CT_mean_methylation_in_tissues <- SummarizedExperiment(assays = mat,
                                                       rowData = methylation_analysis)

usethis::use_data(CT_mean_methylation_in_tissues, overwrite = TRUE)
