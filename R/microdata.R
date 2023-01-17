#' Microbiome Data from Article Supplement
#'
#' This file was created from the file `microbiome_data.xlsx` in the supplementary information of the log-ratio draft manuscript.
#' It contains two more measurements than the `microbiome_data.xlsx` file in the supplementary information of \insertCite{scealy2022sc}{scorecompdir},
#' which were deemed outliers by Scealy and Wood, and have `IndividualID` of `2079` and `2280`.
#' 
#' This data is currently private, pending permission from the original creators of the data.
#'
#' The file was created `microbiome_data.xlsx` using the below code
#'
#' ```
#' microdata <- readxl::read_excel("microbiome_data.xlsx")
#' microdata <- as.data.frame(microdata)
#' save(microdata, file = "./data/microdata.rda")
#' ```
#'
#' @format A `data.frame` with 94 rows and 31 columns.
#' @source To be confirmed.
#' @references
#' \insertAllCited{}
"microdata"
