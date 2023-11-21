# @title Helper function to get install information about this package
# @return A data frame of information

getinstallinfo <- function(){
scorecompdirinfo <- utils::packageDescription("scorecompdir")
pkgimportant <- paste(scorecompdirinfo[c("Suggests", "Imports", "Depends")], collapse = ", ")
pkgimportant <- unlist(strsplit(gsub("[^,ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz]", "", pkgimportant), ","))
pkgimportant <- c(pkgimportant, "scorecompdir")
pkginfo <- as.data.frame(utils::installed.packages(fields = c("Built", "Packaged")))
pkginfo <- pkginfo[pkginfo$Package %in% pkgimportant, 
                   c("Package", "LibPath", "Version",
                     "Built", "Packaged")]
return(pkginfo)
}

