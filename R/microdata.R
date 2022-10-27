#' Microbiome Data from Article Supplement
#'
#' This file was created from the file `microbiome_data.xlsx` in the supplementary information for the log-ratio paper.
#' It two extra measurements than the `microbiome_data.xlsx` file in the supplementary information for JASA 2022 paper.
#' These two extra measurements are outliers and have `IndividualID` of `2079` and `2280`.
#'
#' The file was created using the below
#'
#' ```
#' microdata <- readxl::read_excel("microbiome_data.xlsx")
#' microdata <- as.data.frame(microdata)
#' save(microdata, file = "./data/microdata.rda")
#' ```
#'
#' @format A `data.frame` with 94 rows and 31 columns.
#' @source is there a source?
"microdata"
