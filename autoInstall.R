
if(!require("dada2",quietly=TRUE,warn.conflicts =FALSE)){
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager",repos="https://cloud.r-project.org", version = "3.11")
  }
  if (!requireNamespace("dada2", quietly = TRUE))
    BiocManager::install("dada2", version = "3.11")
} else {
  print("dada2 is already installed.")
}

if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){
  print("phyloseq")
  source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local =TRUE)
} else {
  print("phyloseq is already installed.")
}

my_packages <- c("phyloseq", "dada2") 

not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]

if (length(not_installed)>0) {
  cat(paste("Package", not_installed, "could not be installed. Please install it manually in your R environment",sep="\t"))
}
