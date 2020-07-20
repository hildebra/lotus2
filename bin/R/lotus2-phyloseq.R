
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

  stop("Please run the taxonomic classification first", call.=FALSE)
  
  }



# Require the pyloseq package:
if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local=TRUE);require("phyloseq")}


#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local=TRUE)

library("phyloseq")
packageVersion("phyloseq")
# ‘1.30.0'


cat("\n\nPhyloseq object is being created...\n\n");
  
  
# Read the OTU or ASV table:
otu=read.table(path_TABLE,row.names=1,head=T)
otu=as.matrix(otu)
otu1 <- otu_table(otu, taxa_are_rows=FALSE)

# Read the metadata:
sd=read.table(text = gsub(",", "\t", readLines(path_SD)))
sam1 <- sample_data(sd) 

# Read the taxonomy table:
tax= read.table(path_TAX, row.names=1,head=T,sep="\t")
tax=as.matrix(tax)
tax1 <- tax_table(tax)

# Create the phyloseq object
 physeq <- phyloseq(otu1,tax1,sam1)
 
# Save the phyloseq object as an R object:
saveRDS(physeq,file=paste0(strsplit(path_TABLE,"O")[[1]][1],"phyloseq_object.R"))


  cat("\n\nPhyloseq object is created: phyloseq_object.R\n\n");