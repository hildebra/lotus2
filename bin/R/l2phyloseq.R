
# Require the pyloseq package:
if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){
	source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local=TRUE);
	require("phyloseq")
}
if(!require("ape",quietly=TRUE,warn.conflicts =FALSE)){
	install.packages("ape",repos="https://cloud.r-project.org",quiet=TRUE);
}

if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){#no use..
	cat("Could not find phyloseq R package, possibly a problem installing this.\nWill not write phyloseq package.\n")
	q("n")
}

library("phyloseq")
library("ape")
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
keepUnclass=args[4] #0 or anything else than 0
path_TREE=args[5]


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


cat("\nCreating Phyloseq object ...\n");
  
  
# Read the OTU or ASV table:
#otu=as.matrix(read.table(path_TABLE,row.names=1,header=TRUE,sep="\t"))

otu <- read.delim(path_TABLE,row.names = 1,head=TRUE,sep="\t",check.names=FALSE)

tree=NULL
if((file.exists(path_TREE))) {
	tree=read.tree(path_TREE)
	if (sum(!tree$tip.label %in% rownames(otu)) >0){
		rmS=tree$tip.label[tree$tip.label %in% rownames(otu)]
		tree = drop.tip(tree,names(rmS))
	}
}

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

#First let's check whether any of the samples names starts with a number. If so, phyloseq does not include these samples in the phyloseq object
#If so, we will add "S_" in the beginning of the sample name:
sample_names=rownames(sd)
start_nr <- grep("^\\d", sample_names)
new_names=replace(sample_names, start_nr, paste0("S_", sample_names[start_nr]))
rownames(sd)=new_names

if(length(start_nr)!=0){sample_names_otu=colnames(otu); start_nr_otu <- grep("^\\d", sample_names_otu)
new_names_otu=replace(sample_names_otu, start_nr_otu, paste0("S_", sample_names_otu[start_nr_otu]))
colnames(otu)=new_names_otu}


if(!(all(sample_names==new_names))) {cat("WARNING: Phyloseq doesn't recognize sample names starting with a number. Some of the sample names start with a number..\"S_\" is added into beginning of these sample names.\n")}

#actual conversion
colnames(tax)=c("Domain", "Phylum","Class","Order","Family","Genus","Species")
sam1 <- sample_data(sd) 
otu1 <- otu_table(otu, taxa_are_rows=TRUE)
tax1 <- tax_table(tax)


# Create the phyloseq object
physeq=NULL
if (is.null(tree)) {
	physeq = phyloseq(otu1,tax1,sam1)
} else {
	physeq = phyloseq(otu1,tax1,sam1,tree)
}

if(keepUnclass==0 &&!is.null(physeq)){
	physeq=subset_taxa(physeq, Domain != "?")
}
# Save the phyloseq object as an R object:
save(physeq,file=paste0(dirname(path_TABLE),"/phyloseq.Rdata"))


cat("Phyloseq object is created: phyloseq.Rdata\n")
