
#How a matchlist is obtained (matchid 0.99): vsearch --usearch_global otus.fa --db otus.fa --self --id .99 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10

#We have three input files: matchlist, otutable and hiera_blast:
args = commandArgs(trailingOnly=TRUE)

path_TABLE=args[1]
path_TAX=args[2]
path_matchlist=args[3]

otutable <- read.delim(path_TABLE) #rownames are ASV names
hiera_BLAST <- read.delim(path_TAX) #First columns is "OTU" names
matchlist <- read.delim(path_matchlist, header=FALSE)
colnames(matchlist)=c("OTUid","hit","match")


#Make correlation matrix of the otu counts:
corM=cor(t(otutable))
diag(corM)=0


#Clustering criteria: co-occurence (total overlap) & co-abundance (0.95) & no dramatic diff btw counts of two ASVs (1.5) & matching tax ("Genus"):
out_list=as.list(1:dim(otutable)[1])
for (line in 1:dim(otutable)[1]){
  
  #print(paste0("progress: ",
  #round(((line/nrow(statistics_table)) * 100), 0), "%"))
  potential_ids=as.vector(matchlist[matchlist[,"OTUid"]==rownames(otutable)[line],"hit"])
  potential_parents=match(potential_ids,rownames(otutable))
  daughter_samples <- otutable[line, ]
  
  # Only include ASVs having >0.95 match in vsearch:
  shouldbemerged=which(corM[line,]>0.95)
 
  
  for (i in potential_parents){
    rel_cooccur=sum((daughter_samples[otutable[i, ] > 0]) > 0)/sum(daughter_samples > 0)
    relative_abund=min(otutable[i, ][daughter_samples > 0]/daughter_samples[daughter_samples > 0])
    inv_relative_abund=min(daughter_samples[daughter_samples > 0]/otutable[i, ][daughter_samples > 0])
    highest_rank=min((length(which(hiera_BLAST[line,]!="?"))-1),(length(which(hiera_BLAST[i,]!="?"))-1))
   
     if(
      rel_cooccur==1
       && sum(as.numeric(which(otutable[line,]==0)!=which(otutable[i,]==0)))==0 #check if zero abundance happens in the same samples
       && (relative_abund<1.5) && (inv_relative_abund<1.5)
       &&  highest_rank>=6 #Match the tax at least at the Genus level
       &&  hiera_BLAST[which(hiera_BLAST[,"OTU"]==row.names(otutable)[line]),highest_rank]==hiera_BLAST[which(hiera_BLAST[,"OTU"]==row.names(otutable)[i]),highest_rank]
       &&  length(shouldbemerged)!=0
       &&  any(names(corM[line,shouldbemerged])%in%rownames(otutable)[i]))#if the potential parent has >0.95 correlation
       {out_list[[line]]=c(out_list[[line]],i) 
       
   }
    else {next}
  }}
names(out_list)=rownames(otutable)[1:length(out_list)]


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
for (i in 1:length(out_list_nodup)){
  for (j in 1:length(out_list_nodup)){
if(any(out_list_nodup[[i]]%in%out_list_nodup[[j]])){
  out_list_nodup[[i]]=union(out_list_nodup[[i]],out_list_nodup[[j]])
    }}}


#Removing duplicate rows randomly:
dup_combn=combn(1:length(out_list_nodup),m=2)
mat=data.frame(matrix(ncol=2,,));colnames(mat)=c("v1","v2");for (i in (1:dim(dup_combn)[2])){if(setequal(out_list_nodup[[dup_combn[1,i]]],out_list_nodup[[dup_combn[2,i]]])){match=cbind(dup_combn[1,i],dup_combn[2,i]);colnames(match)=c("v1","v2");mat=rbind(mat,match)}}
mat=mat[-1,]

#Choose one of the duplicates randomly and drop it:
random=c();for (i in (1:dim(mat)[1])){random=c(random,sample(mat[i,],1))}
random=as.vector(unlist(random))
out_list_nodup2=out_list_nodup[-random]


#Summed up the abundances and assigned them to the first index, which is the parent ASV (ASVs summed up in the otutable).
for (i in out_list_nodup2){otutable[i[1],]=colSums(otutable[i,])}
#Assigned 0 for the daughter ASVs in the otutable:
for (i in out_list_nodup2){otutable[i[-1],]=0} 
#And removed them both from the otutable and tax table:
otutable_denoised=otutable[which(rowSums(otutable)!=0),]
hiera_BLAST_new=hiera_BLAST[(hiera_BLAST[,"OTU"]%in%rownames(otutable_denoised)),]

merged_table=matrix(ncol=2,nrow=length(out_list_nodup2));for(i in 1:length(out_list_nodup2)){
  merged_table[i,1]=rownames(otutable)[out_list_nodup2[[i]]][1]
  merged_table[i,2]=paste0(rownames(otutable)[out_list_nodup2[[i]]],collapse = ",")
}
colnames(merged_table)=c("Parent_OTU","Merged_OTUs")
colnames(hiera_BLAST_new)[1]=c("Parent_OTU")
#The taxa of the parent ASV is assigned to the merged ASVs:
merged_table=merge(merged_table,hiera_BLAST_new,by="Parent OTU")
colnames(hiera_BLAST_new)[1]=c("OTU")

#Write the output files:
write.table(merged_table,file="clustered_ASVs.txt",row.names=FALSE,sep="\t")
write.table(otutable_denoised,file="OTU.txt",sep="\t")
write.table(hiera_BLAST_new,file="hiera_BLAST.txt",row.names=FALSE,sep="\t")

#otutable_denoised and hiera_BLAST_new can be passed into the phyloseq object as an argument.
#However, phylogenetic tree should be generated after this script, since some ASVs will be removed and this will make problem while mergeing the phyloseq object otherwise.
#That's I did not bind this script and l2phyloseq.R together.

