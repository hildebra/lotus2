
#How a matchlist is obtained (matchid 0.84): vsearch --usearch_global OTU.fna --db OTU.fna --self --id .97 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

#We have five/six input files:
args = commandArgs(trailingOnly=TRUE)

path_TABLE=args[1]
path_TAX=args[2]
path_matchlist=args[3]
path_SD=args[4]
keepUnclass=args[5] #0 or anything else than 0
path_TREE=args[6]


# Test if taxonomy table is produced: 
RDPhier=FALSE
if(!(file.exists(path_TAX))) {
	  #stop(paste("Please run the taxonomic classification first:",path_TAX), call.=FALSE)
	  path_TAX = gsub("hiera_BLAST.txt","hiera_RDP.txt",path_TAX)
	  RDPhier=TRUE
	  if(!(file.exists(path_TAX))) { #still doesn't exist
		cat("Could not find taxonomy file.. aborting ASV_clustering.R\n")
		q("no")
	  }
}


otutable <- read.delim(path_TABLE) #rownames are ASV names
hiera_BLAST <- read.delim(path_TAX) #First columns is "OTU" names


#check if file is empty
if (!file.exists(path_matchlist) || file.info(path_matchlist)$size == 0){
	#no matches.. probably nothing to do..
	cat("No sequence matches found, nothing to do for ASV_clustering.R\n")
	q("no")
}
matchlist <- read.delim(path_matchlist, header=FALSE)
colnames(matchlist)=c("OTUid","hit","match")

otutab=merge(otutable,hiera_BLAST,by="OTU")
rownames(otutab)=otutab[,"OTU"]
tab=otutab[,2:dim(otutable)[2]]
tax=otutab[,(dim(otutab)[2]-6):dim(otutab)[2]]
tax[,"OTU"]=rownames(tax)


#Make correlation matrix of the otu counts:
corM=cor(t(tab))
diag(corM)=0


#Clustering criteria: ASVs having at least 10 reads in total, co-occurence (overlap) & co-abundance (0.70) & no dramatic diff btw counts of two ASVs (<2) & matching tax ("Family"):
out_list=as.list(1:dim(tab)[1])
for (line in 1:dim(tab)[1]){

	potential_ids=as.vector(matchlist[matchlist[,"OTUid"]==rownames(tab)[line],"hit"])
	potential_parents=match(potential_ids,rownames(tab))
	daughter_samples <- tab[line, ]

	# Only include ASVs showing >0.70 corr in abundances:
	shouldbemerged=which(corM[line,]>0.70)


	for (i in potential_parents){
		# this co-occurrence did not take into account the samples of the second vector? [ "/ sum(..)"
		#rel_cooccur=sum((daughter_samples[tab[i, ] > 0]) > 0)/sum(daughter_samples > 0)
		rel_cooccur = sum(tab[line, ]>0 & tab[i, ]>0)   /  sum(tab[line, ] + tab[i, ]>0 )
		
		relcoabsence = sum(tab[line, ]==0 & tab[i, ]==0)   /  sum(tab[line, ] + tab[i, ]==0 )
		
		#suppressWarnings : prob makes this slower, you could just vote to ignore warnings?  ##Do you mean to suppress them globally?
		relative_abund=suppressWarnings(min(tab[i, ][daughter_samples > 0]/daughter_samples[daughter_samples > 0]))
		inv_relative_abund=suppressWarnings(min(daughter_samples[daughter_samples > 0]/tab[i, ][daughter_samples > 0]))
		highest_rank=min((length(which(tax[line,1:7]!="?"))),(length(which(tax[i,1:7]!="?"))))
		zero_line=as.numeric(which(tab[line,]==0))
		zero_i=which(tab[i,]==0)


		#preselection based on some arbitrary criteria
		if(rel_cooccur>=0.40
		   && sum(tab[line,])>10  #needed?? if there are ASVs occurring in very low abundance, they are usually identified as being copies; however they are probably just noise.
		   
		   #this term could be too hard on actually abundant samples.. I inserted a value for this above: relcoabsence >= 0.8
		   && relcoabsence >= 0.8
		   #&& length(intersect(zero_line,zero_i))/dim(tab)[2]>=0.80 #check if zero abundance happens in most of the same samples
		   && (relative_abund<2) && (inv_relative_abund<2)
		   &&  highest_rank>=5 #Match the tax at least at the Family level
		   &&  tax[which(tax[,"OTU"]==row.names(tab)[line]),highest_rank]==tax[which(tax[,"OTU"]==row.names(tab)[i]),highest_rank]
		   &&  length(shouldbemerged)!=0
		   &&  any(names(corM[line,shouldbemerged])%in%rownames(tab)[i]))#if the potential parent has >0.70 correlation
		{
			out_list[[line]]=c(out_list[[line]],i) 
		}

	}
}
names(out_list)=rownames(tab)[1:length(out_list)]


#Remove ASVs that did not match with any other ASV:
out_list2=out_list[lengths(out_list) != 1]

#Find the duplicate rows:
a=combn(1:length(out_list2),m=2)
mat=data.frame(matrix(ncol=2,,));colnames(mat)=c("v1","v2");for (i in (1:dim(a)[2])){if(setequal(out_list2[[a[1,i]]],out_list2[[a[2,i]]])){match=cbind(a[1,i],a[2,i]);colnames(match)=c("v1","v2");mat=rbind(mat,match)}}
mat=mat[-1,]

#Choose one of the duplicates randomly and drop it:
random=c();for (i in (1:dim(mat)[1])){random=c(random,sample(mat[i,],1))}
random=as.vector(unlist(random))
out_list_nodup=out_list2[-random]

#Single linkage clustering:
union_index=list()
for (i in 1:length(out_list_nodup)){
  for (j in 1:length(out_list_nodup)){
    if(any(out_list_nodup[[i]]%in%out_list_nodup[[j]])){
      out_list_nodup[[i]]=union(out_list_nodup[[i]],out_list_nodup[[j]]); union_index[[i]]=c(union_index,i,j)}
    else {next}
  }}

#Removing duplicate rows randomly:
out_list_nodup2=out_list_nodup[!(duplicated(lapply(out_list_nodup,sort)))]


#Summed up the abundances and assigned them to the first index, which is the parent ASV (ASVs summed up in the tab).
for (i in out_list_nodup2){tab[i[1],]=colSums(tab[i,])}
#Assigned 0 for the daughter ASVs in the tab:
for (i in out_list_nodup2){tab[i[-1],]=0} 
#And removed them both from the tab and tax table:
tab_denoised=tab[which(rowSums(tab)!=0),]
tax_new=tax[(tax[,"OTU"]%in%rownames(tab_denoised)),]

merged_table=matrix(ncol=2,nrow=length(out_list_nodup2));for(i in 1:length(out_list_nodup2)){
  merged_table[i,1]=rownames(tab)[out_list_nodup2[[i]]][1]
  merged_table[i,2]=paste0(rownames(tab)[out_list_nodup2[[i]]],collapse = ",")
}
colnames(merged_table)=c("Parent_OTU","Merged_OTUs")
colnames(tax_new)[8]=c("Parent_OTU")
#The taxa of the parent ASV is assigned to the merged ASVs:
merged_table=merge(merged_table,tax_new,by="Parent_OTU")
colnames(tax_new)[8]=c("OTU")

#Write the output files:
write.table(merged_table,file="clustered_ASVs.txt",row.names=FALSE,sep="\t")
write.table(tab_denoised,file="OTU.txt",sep="\t")
write.table(tax_new[,-8],file="hiera_BLAST.txt",row.names=TRUE,sep="\t")

#### Create a phyloseq object ###:
cat("\nPhyloseq object is being created...\n");


# Require the pyloseq package:
if(!require("phyloseq",quietly=TRUE,warn.conflicts =FALSE)){
  source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local=TRUE);
  require("phyloseq")
}
if(!require("ape",quietly=TRUE,warn.conflicts =FALSE)){
  install.packages("ape",repos="https://cloud.r-project.org",quiet=TRUE);require(ape)
}

library("phyloseq")

# Test if phylogenetic tree is produced:
# If the phylogenetic tree exists, require the "ape" package and read the tree:
tree=NULL
if((file.exists(path_TREE))) {
  tree=read.tree(path_TREE)
}

#Read the sample data:
sdA=scan(file=path_SD,nlines =1,sep ="\t",what="character")
sd=read.table(path_SD,sep="\t",row.names=1,header=FALSE,comment.char="#",as.is=TRUE)
colnames(sd) = sdA[-1]

otu=tab_denoised
tax_phylo=tax_new[,-8]
if (dim(tax_phylo)[1] < dim(otu)[1]){
  taxA = matrix("?", nrow(otu), ncol(tax_phylo))
  rownames(taxA) = dimnames(otu)[[1]]
  idx = dimnames(otu)[[1]] %in% dimnames(tax)[[1]] 
  taxA[idx,] = tax_phylo[dimnames(otu)[[1]][idx],,drop=FALSE]
  tax_phylo = taxA
  #rownames(taxA) = dimnames(otu)[[1]][idx]
  #tax=rbind (tax, taxA)
}


tax=as.matrix(tax_phylo)


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

#Create the phyloseq object
if (length(args) == 5) {
physeq = phyloseq(otu1,tax1,sam1)}else if (length(args) == 6)
{physeq = phyloseq(otu1,tax1,sam1,tree)}

if(keepUnclass==0){
  physeq=subset_taxa(physeq, Domain != "?")}
# Save the phyloseq object as an R object:
save(physeq,file=paste0(strsplit(path_TABLE,"O")[[1]][1],"phyloseq.Rdata"))
