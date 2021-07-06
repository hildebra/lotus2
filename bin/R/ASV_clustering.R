
#How a matchlist is obtained (matchid 0.84): vsearch --usearch_global OTU.fna --db OTU.fna --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10


#We have three input files: matchlist, otutable and hiera_blast:
args = commandArgs(trailingOnly=TRUE)

path_TABLE=args[1]
path_TAX=args[2]
path_matchlist=args[3]

otutable <- read.delim(path_TABLE) #rownames are ASV names
hiera_BLAST <- read.delim(path_TAX) #First columns is "OTU" names
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
    rel_cooccur=sum((daughter_samples[tab[i, ] > 0]) > 0)/sum(daughter_samples > 0)
    relative_abund=suppressWarnings(min(tab[i, ][daughter_samples > 0]/daughter_samples[daughter_samples > 0]))
    inv_relative_abund=suppressWarnings(min(daughter_samples[daughter_samples > 0]/tab[i, ][daughter_samples > 0]))
    highest_rank=min((length(which(tax[line,1:7]!="?"))),(length(which(tax[i,1:7]!="?"))))
    zero_line=as.numeric(which(tab[line,]==0))
    zero_i=which(tab[i,]==0)
    
    if(rel_cooccur>=0.80
       && sum(tab[line,])>10 
       && length(intersect(zero_line,zero_i))/dim(tab)[2]>=0.80 #check if zero abundance happens in most of the same samples
       && (relative_abund<2) && (inv_relative_abund<2)
       &&  highest_rank>=5 #Match the tax at least at the Family level
       &&  tax[which(tax[,"OTU"]==row.names(tab)[line]),highest_rank]==tax[which(tax[,"OTU"]==row.names(tab)[i]),highest_rank]
       &&  length(shouldbemerged)!=0
       &&  any(names(corM[line,shouldbemerged])%in%rownames(tab)[i]))#if the potential parent has >0.70 correlation
    {out_list[[line]]=c(out_list[[line]],i) 
    
    }
    else {next}
  }}
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

#tab_denoised and tax_new can be passed into the phyloseq object as an argument.
#However, phylogenetic tree should be generated after this script, since some ASVs will be removed and this will make problem while mergeing the phyloseq object otherwise.
#That's why I did not bind this script and l2phyloseq.R together.
