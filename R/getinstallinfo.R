#' @title Helper function to get install information about this package
#' @return A data frame of information

getinstallinfo <- function(){
cdabyppiinfo <- packageDescription("cdabyppi")
pkgimportant <- paste(cdabyppiinfo[c("Suggests", "Imports", "Depends")], collapse = ", ")
pkgimportant <- gsub("[^,ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz]", "", pkgimportant) %>%
  strsplit(",") %>%
  unlist()
pkgimportant <- c(pkgimportant, "cdabyppi")
pkginfo <- as.data.frame(installed.packages(field = c("Built", "Packaged")))
pkginfo <- pkginfo[pkginfo$Package %in% pkgimportant, 
                   c("Package", "LibPath", "Version",
                     "Built", "Packaged")]
return(pkginfo)
}

