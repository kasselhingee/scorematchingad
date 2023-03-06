#' Microbiome Data for Soil-Transmitted Helminths
#'
#' The data `microbiome` contains paired DNA samples before and at 21 months after treatment for helminth infections \insertCite{Martin2019mi}{scorecompdir}.
#' This data was analysed by Martin et al \insertCite{Martin2019mi}{scorecompdir} and a further subset was studied by Scealy et al \insertCite{Scealy2022sc}{scorecompdir}.
#' The data come from a study into the effect of helminth infections on the course of malaria infections (ImmunoSPIN-Malaria) in the Nangapanda subdistrict, Indonesia \insertCite{Wiria2010do}{scorecompdir}.
#' As part of the study, some participants were given 400mg of albendazole every three months for 1.5 years,
#' remaining participants were given a placebo
#'  \insertCite{Wiria2010do}{scorecompdir}.
#'
#' The measurements in the data come from stool samples before and after treatment.
#' Gut microbiome prevalence was measured using 16s rRNA 454 sequencing \insertCite{Martin2019mi}{scorecompdir}.
#' Helminth infections were detected by PCR or microscopy \insertCite{Martin2019mi}{scorecompdir}.
#' 
#' The subset studied by Scealy and Wood  \insertCite{Scealy2022sc}{scorecompdir} contained only the measurements from before treatment, and only those individuals with a helminth infection.
#' These measurements can be obtained by running...
#' Two further individuals (`IndividualID` of `2079` and `2280`) deemed outliers by Scealy and Wood \insertCite{Scealy2022sc}{scorecompdir}, and removed in their analyses.
#'
#' @format 
#' A dataframe with ?? rows and ?? columns:
#' \describe{
#'   \item{IndividualID}{An integer uniquely specifying the individual.}
#'   \item{Year}{The collection year for the sample.}
#'   \item{Sex}{`1` if female, `0` otherwise.}
#'   \item{Treatment}{`1` if individual given 400mg of albendazole every three months for 1.5 years, `0` otherwise.}
#'   \item{Age}{Age at first sample.}
#'   \item{ct_Al}{The qPCR cycle threshold (CT) for *Ascaris lumbricoides* (large roundworm). *Ascaris lumbricoides* can be considered present if the value is 30 or less.}
#'   \item{ct_Na}{The qPCR cycle threshold (CT) for *Necator americanus* (a hookworm). *Necator americanus* can be considered present if the value is 30 or less.}
#'   \item{ct_Ad}{The qPCR cycle threshold (CT) for *Ancylostoma duodenale* (a hookworm). *Ancylostoma duodenale* can be considered present if the value is 30 or less.}
#'   \item{micr_Tt}{The presence of *Trichuris trichiura* as determined by microscopy. A value of `1` means *Trichuris trichiura* was detected.}
#'   \item{helminth}{If any of the above helminths were detected then `1`, otherwise `0`.}
#' }
#' 
#' 16s measurements of stool samples pre and post sample
#' 16s rRNA gene was processed with the 454 pyrosequencing technique
#' phyla
#' treatment = 400mg of albendazole every three months for 1.5 years
#' Nangapanda subdistrict, Indonesia
#' data collected to assess the effect of the treatment on soil-transmitted helminth infections
#' results of first (baseline) measurements of individuals only
#' also ONLY for individuals infected by helminth according to Scealy 2022 (Martin et al also had 94 individuals with helmiths at baseline - 47 in treatment and 47 without treatment)
#' 
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
