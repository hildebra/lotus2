
if(!require("dada2",quietly=TRUE,warn.conflicts =FALSE)){
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager",repos="https://cloud.r-project.org", version = "3.11")
  }
  if (!requireNamespace("dada2", quietly = TRUE))
    BiocManager::install("dada2")
} else {
  #print("dada2 is already installed.")
}

if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){
  #print("phyloseq")
  source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local =TRUE)
} else {
  #print("phyloseq is already installed.")
}


if(!require("dplyr",quietly=TRUE,warn.conflicts =FALSE)){
	install.packages("dplyr",repos="https://cloud.r-project.org")
	#require("dplyr",warn.conflicts =FALSE)
}

if(0 && !require("lulu",quietly=TRUE,warn.conflicts =FALSE)){ #no longer needed.. dplyr now
	if(!require("devtools",quietly=TRUE,warn.conflicts =FALSE)){
		install.packages("devtools",repos="https://cloud.r-project.org",quiet=TRUE);require(devtools)
	}
	library(devtools)
	install_github("tobiasgf/lulu")  
}


my_packages <- c("phyloseq", "dada2", "dplyr") 

not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]

if (length(not_installed)>0) {
  cat(paste("Package", not_installed, "could not be installed. Please install it manually in your R environment",sep="\t"))
}
