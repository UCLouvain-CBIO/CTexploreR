#' Check spelling of entered variables
#'
#' @description
#'
#' Checks the spelling of a vector of entered variable(s) comparing it
#' to a vector of valid names, and removes the ones that are absent
#' from the vector of valid names.
#'
#' @param variable `character()` containing the names of variables to
#'     check.
#'
#' @param valid_vector `character()` with valid variable names.
#'
#' @return A character with valid variables.
#'
#' @examples
#' CTexploreR:::check_names(variable = c("Ovarian", "leukemia", "wrong_name"),
#'                          valid_vector = c("ovarian", "leukemia"))
check_names <- function(variable, valid_vector) {
    in_valid <- variable %in% valid_vector
    if (!all(in_valid)) {
        msg <- paste0(sum(!in_valid), " out of ",
                     length(in_valid), " names invalid: ",
                     paste(variable[!in_valid], collapse = ", "),
                     ".")
        warning(paste(strwrap(msg), collapse = "\n"),
                "\nSee the manual page for valid types.",
                call. = FALSE)
    }
    variable[in_valid]
}

## Color vectors for heatmaps

legend_colors <- c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
                   "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F",
                   "#9E0142")


CCLE_colors <- c("lung" = "seagreen3", "skin" = "red3",
                 "bile_duct" = "mediumpurple1", "bladder" = "mistyrose2",
                 "colorectal" = "plum", "lymphoma" = "steelblue1",
                 "uterine" = "darkorange4", "myeloma" = "turquoise3", 
                 "kidney" = "thistle4",
                 "pancreatic" = "darkmagenta", "brain" = "palegreen2",
                 "gastric" = "wheat3", "breast" = "midnightblue",
                 "bone" = "sienna1", "head_and_neck" = "deeppink2",
                 "ovarian" = "tan3", "sarcoma" = "lightcoral",
                 "leukemia" = "steelblue4", "esophageal"= "khaki",
                 "neuroblastoma" = "olivedrab1")


DAC_colors <- c("B2-1" = "olivedrab2", "HCT116" = "lightcoral",
                "HEK293T" = "seagreen3", "HMLER" = "mediumpurple1",
                "IMR5-75" = "deeppink2", "NCH1681" = "steelblue2",
                "NCH612" = "red3", "TS603" = "darkmagenta")


TCGA_colors <- c("BRCA" = "midnightblue", "COAD" = "darkorchid2",
                 "ESCA" = "gold", "HNSC" = "deeppink2",
                 "LUAD" = "seagreen", "LUSC" = "seagreen3",
                 "SKCM" = "red3")

testis_colors <- c("SSC" = "floralwhite", "Spermatogonia" = "moccasin",
                   "Early_spermatocyte" = "gold", 
                   "Late_spermatocyte" = "orange",
                   "Round_spermatid" = "red2", 
                   "Elongated_spermatid" = "darkred",
                   "Sperm1" = "violet", "Sperm2" = "purple", 
                   "Sertoli" = "gray", 
                   "Leydig" = "cadetblue2", "Myoid" = "springgreen3", 
                   "Macrophage" = "gray10",
                   "Endothelial" = "steelblue")

cell_type_colors <- c("Germ_cells" = "mediumvioletred", 
                      "Trophoblast_cells" = "brown4",
                      "Adipocyte_cells" = "lightgoldenrod1", 
                      "Blood_immune_cells" = "ivory2",
                      "Endocrine_cells" = "lightslateblue", 
                      "Endothelial_cells" = "plum2",
                      "Glandular_cells" = "deepskyblue", 
                      "Glial_cells" = "sienna1", 
                      "Mesenchymal_cells" = "peachpuff", 
                      "Pigment_cells" = "palegreen2", "Myoid" = "gray54", 
                      "Specialized_epithelial_cells" = "indianred1")

#' CT genes description table
#'
#' @description
#'
#' Cancer-Testis (CT) genes description, imported from `CTdata`
#'
#' @format
#'
#' A `tibble` object with 308 rows and 37 columns.
#'
#' - Rows correspond to CT genes
#'
#' - Columns give CT genes characteristics
#'
#' @details
#'
#' See `CTdata::CT_genes` documentation for details
#' 
#' @source
#'
#' See `scripts/make_CT_genes.R` in `CTdata` for details on how this list of
#' curated CT genes was created.
#' 
#' @export
#' 
#' @name CT_genes
#'
#' @docType data
load("inst/extdata/CT_genes.rda")

# CT_genes <- CTdata::CT_genes()


