
# Require the pyloseq package:
if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){
	source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local=TRUE);
	require("phyloseq")
}

library("phyloseq")
#packageVersion("phyloseq")
#path_TABLE=$outdir/OTU.txt
#Classification_method   #phyloseq option needs to have a flag specifying  the output of the classification to use (for ex BLAST)
#path_TAX=$outdir/*$Classification_method*.txt   
#path_SD=$map
#RScript lotus2-phyloseq.R $outdir/OTU.tx $outdir/*$Classification_method*.txt $map

args = commandArgs(trailingOnly=TRUE)

path_TABLE=args[1]
path_TAX=args[2]
path_SD=args[3]
path_TREE=args[4]


# Test if taxonomy table is produced: 
RDPhier=FALSE
if(!(file.exists(path_TAX))) {
	#stop(paste("Please run the taxonomic classification first:",path_TAX), call.=FALSE)
	path_TAX = gsub("hiera_BLAST.txt","hiera_RDP.txt",path_TAX)
	RDPhier=TRUE
	if(!(file.exists(path_TAX))) { #still doesn't exist
		cat("Could not find taxonomy file.. aborting l2phyloseq\n")
		q("no")
	}
}

# Test if phylogenetic tree is produced:
# If the phylogenetic tree exists, require the "ape" package and read the tree:
if((file.exists(path_TREE))) {
	if(!require("ape",quietly=TRUE,warn.conflicts =FALSE)){
		install.packages("ape",repos="https://cloud.r-project.org",quiet=TRUE);require(ape)
	}
	tree=read.tree(path_TREE)
}


cat("\nPhyloseq object is being created...\n");
  
  
# Read the OTU or ASV table:
#otu=as.matrix(read.table(path_TABLE,row.names=1,header=TRUE,sep="\t"))

otu <- read.delim(path_TABLE,row.names = 1,head=TRUE,sep="\t")

# Read the metadata:
#sd=read.table(text = gsub(",", "\t", readLines(path_SD)))
sdA=scan(file=path_SD,nlines =1,sep ="\t",what="character")
sd=read.table(path_SD,sep="\t",row.names=1,header=FALSE,comment.char="#",as.is=TRUE)
colnames(sd) = sdA[-1]

# Read the taxonomy table:
tax= as.matrix(read.delim(path_TAX, row.names=NULL,header=TRUE,sep="\t") )
if (RDPhier){
	rownames(tax) = tax[,dim(tax)[2],drop=FALSE]
	tax = tax[,1:(dim(tax)[2]-1),drop=FALSE]
} else {
	rownames(tax) = tax[,1]
	tax = tax[,2:(dim(tax)[2]),drop=FALSE]
}

if (dim(tax)[1] < dim(otu)[1]){
	taxA = matrix("?", nrow(otu), ncol(tax))
	rownames(taxA) = dimnames(otu)[[1]]
	idx = dimnames(otu)[[1]] %in% dimnames(tax)[[1]] 
	taxA[idx,] = tax[dimnames(otu)[[1]][idx],,drop=FALSE]
	tax = taxA
	#rownames(taxA) = dimnames(otu)[[1]][idx]
	#tax=rbind (tax, taxA)
}
#actual conversion
colnames(tax)=c("Domain", "Phylum","Class","Order","Family","Genus","Species")
sam1 <- sample_data(sd) 
otu1 <- otu_table(otu, taxa_are_rows=TRUE)
tax1 <- tax_table(tax)


# Create the phyloseq object
if (length(args) == 3) {
  physeq = phyloseq(otu1,tax1,sam1)}else if (length(args) == 4)
      {physeq = phyloseq(otu1,tax1,sam1,tree)}
 
physeq=subset_taxa(physeq, Domain != "?")
# Save the phyloseq object as an R object:
save(physeq,file=paste0(strsplit(path_TABLE,"O")[[1]][1],"phyloseq.Rdata"))


  cat("Phyloseq object is created: phyloseq.Rdata\n")
