
#path_TABLE=$outdir/OTU.txt
#Classification_method   #phyloseq option needs to have a flag specifying  the output of the classification to use (for ex BLAST)
#path_TAX=$outdir/*$Classification_method*.txt   
#path_SD=$map

#RScript lotus2-phyloseq.R $outdir/OTU.tx $outdir/*$Classification_method*.txt $map




args = commandArgs(trailingOnly=TRUE)

path_TABLE=args[1]
path_TAX=args[2]
path_SD=args[3]



# Test if taxonomy table is produced: 

if(!(file.exists(path_TAX))) {
	stop(paste("Please run the taxonomic classification first:",path_TAX), call.=FALSE)
}



# Require the pyloseq package:
if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local=TRUE);require("phyloseq")}

library("phyloseq")
packageVersion("phyloseq")
# ‘1.30.0'


cat("\nPhyloseq object is being created...\n");
  
  
# Read the OTU or ASV table:
otu=as.matrix(read.table(path_TABLE,row.names=1,header=TRUE,sep="\t"))

# Read the metadata:
#sd=read.table(text = gsub(",", "\t", readLines(path_SD)))
sdA=scan(file=path_SD,nlines =1,sep ="\t",what="character")
sd=read.table(path_SD,sep="\t",row.names=1,header=FALSE,comment.char="#",as.is=TRUE)
colnames(sd) = sdA[-1]

# Read the taxonomy table:
tax= as.matrix(read.table(path_TAX, row.names=1,header=TRUE,sep="\t") )
if (dim(tax)[1] < dim(otu)[1]){
	taxA = matrix("?", nrow(otu)-nrow(tax), ncol(tax))
	idx = !dimnames(otu)[[1]] %in% dimnames(tax)[[1]] 
	dimnames(taxA)[[1]] = dimnames(otu)[[1]][idx]
	tax=rbind (tax, taxA)
}
#actual conversion
sam1 <- sample_data(sd) 
otu1 <- otu_table(otu, taxa_are_rows=TRUE)
tax1 <- tax_table(tax)

# Create the phyloseq object
 physeq <- phyloseq(otu1,tax1,sam1)
 
# Save the phyloseq object as an R object:
saveRDS(physeq,file=paste0(strsplit(path_TABLE,"O")[[1]][1],"phyloseq.Rdata"))


  cat("Phyloseq object is created: phyloseq.Rdata\n");