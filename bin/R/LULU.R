#uses LULU to clean OTU tables 

if(!require("lulu",quietly=TRUE,warn.conflicts =FALSE)){
	if(!require("devtools",quietly=TRUE,warn.conflicts =FALSE)){
		install.packages("devtools",repos="https://cloud.r-project.org",quiet=TRUE);require(devtools)
	}
	library(devtools)
	install_github("tobiasgf/lulu")  
	require("lulu")
}
library(lulu)

#args=c("/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//tmpFiles//lulu_match_list.txt","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//OTU.txt","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S/")
args = commandArgs(trailingOnly=TRUE)
matchL = args[1]
otuF = args[2]
logD = args[3]

matchs = read.table(matchL,header=TRUE,as.is=TRUE)
otuM = read.table(otuF,header=TRUE,as.is=TRUE,row.names=1)

lulu <- lulu(otuM, matchs)
write.table(lulu$curated_table,quote =FALSE,sep="\t",file=otuF)
write.table(lulu$discarded_otus,quote =FALSE,sep="\t",file=paste0(matchL,".rm"),col.names =FALSE,row.names =FALSE)
#save the lulu object for log purposes
save(lulu,file=paste0(logD,"LULU.Rdata"));