
if(!require("dada2",quietly=TRUE,warn.conflicts =FALSE)){
	if (!requireNamespace("BiocManager", quietly = TRUE)){
		install.packages("BiocManager",repos="https://cloud.r-project.org")
	}
	BiocManager::install("dada2")
}
#packageVersion("dada2")

#ARGS parsing
#args=c("/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//demultiplexed/","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//demultiplexed/","0","12","/hpc-home/hildebra/dev/lotus/maps/AngeTest1.16S.sm.map")
args = commandArgs(trailingOnly=TRUE)
# test if all the arguments are there: 
if (length(args) <= 4) {
	stop("All the arguments must be supplied: pathF path_output seed_num ncores map_path.\n", call.=FALSE)
}


pathF=args[1]
path_output = args[2]
seed_num = as.integer(args[3])
ncores = as.integer(args[4])
map_path=args[5]


mapping=as.matrix(read.delim(map_path,check.names = FALSE,as.is=TRUE,header=TRUE,sep="\t"))


subset_map=array("A",dim(mapping)[1])
if ("SequencingRun" %in% colnames(mapping)){ #explicitly defined groups of sequencing runs
	subset_map = mapping[,"SequencingRun"]
#	for (i in 1:length(unique(mapping[,"SequencingRun"]))){
#		subset_map[[i]]=subset(mapping,SequencingRun==row.names(table(mapping[,"SequencingRun"]))[i])
#	}
} else if ("fastqFile" %in% colnames(mapping) ){ #1)check if different BCs 2)check if different folders to impute seq runs
	tfqFs = table(mapping[,"fastqFile"])
	hasSlash=grepl("/",names(tfqFs))
	if (any(tfqFs) > 1){ #demultiplexing included
		subset_map = mapping[,"fastqFile"]
		cat(paste("Found ",length(tfqFs)," unique fastqFile's, assumming each represents a sequencing run\n\n",sep=""))
	} else if (any(hasSlash)) {
		leastS=function (x){lx=length(x);if (lx==1){""} else {paste(x[1:(lx-1)],collapse="/")}}
		fastqP = unlist(lapply(strsplit(mapping[,"fastqFile"],"/"),leastS))
		tfqFs=table(fastqP)
		if (length(tfqFs) > 1){
			cat(paste("Found ",length(tfqFs)," unique path's to fastq's, assumming each path represents a sequencing run\n\n",sep=""))
			subset_map = fastqP
		}
	}
} else {
	cat("No sample run information, samples are considered as sequenced in the same run\n\n")
}



# File parsing
listF=list();listR=list()
tSuSe = table(subset_map)
for (i in names(tSuSe)){
	#forward read error pattern
	idx = subset_map == i
	pattern1 = paste0(mapping[idx,"#SampleID"],".1.fq[.gz]?",sep="")
	listF[[i]] = list.files(pathF, pattern=paste0(pattern1, collapse="|"),full.names = TRUE)
	#reverse read
	pattern1 = paste0(mapping[idx,"#SampleID"],".2.fq[.gz]?",sep="")
	listR[[i]] = list.files(pathF, pattern=paste0(pattern1, collapse="|"),full.names = T)
}   

cat("\n\nStarting DADA2 ASV clustering... please be patient\n\n");

list_seqtabs=list()
list_seqtabs.nochim=list()
i=1

for (i in 1:length(listR)){
	sampleNames <- sapply(strsplit(basename(listF[[i]]), ".1.fq"),`[`, 1) # 
	sampleNamesR <- sapply(strsplit(basename(listR[[i]]), ".2.fq"),`[`, 1) #
	if(!identical(sampleNames, sampleNamesR)) stop("Forward and reverse files do not match:",sampleNames,sampleNamesR)
	cat(paste0("Detected ",length(sampleNames)," samples .. Computing error profiles\n"));
	#start of dada2 learning
	set.seed(seed_num)
	filtFs <- file.path(path_output, "matched_IDs", paste0(sampleNames, "_F_matched.fastq.gz"))
	filtRs <- file.path(path_output, "matched_IDs", paste0(sampleNames, "_R_matched.fastq.gz"))
	#filtFs <- file.path(listF[[i]])
	names(filtFs) <- sampleNames
	names(filtRs) <- sampleNames 

	#Not doing any filtering but just to match the IDs of reverse and forward reads:
	filterAndTrim(fwd=listF[[i]], filt=filtFs, rev=listR[[i]], filt.rev=filtRs,matchIDs=TRUE,multithread=ncores)

	# Learn forward error rates
	cat(paste0("Learning error profiles for the forward reads:\n"));
	
	try( #dada2 is too instable
		errF <- learnErrors(filtFs, multithread=ncores)
	)
	if (is.null(errF)) {
		errF <- learnErrors(filtFs, multithread=1)
	}
	
	

	# Learn reverse error rates
	cat(paste0("Learning error profiles for the reverse reads:\n"));	
	#dada2 is too instable
	try( 
		errR <- learnErrors(filtRs, multithread=ncores)
	)
	if (is.null(errF)) {
		errR <- learnErrors(filtRs, multithread=1)
	}

	# Save the plots of error profiles, for a sanity check:
	pdf(paste0(path_output,"/dada2_p",i,"_errF.pdf"),useDingbats = FALSE)
	plotErrors(errF, nominalQ=TRUE)
	dev.off()
	pdf(paste0(path_output,"/dada2_p",i,"_errR.pdf"),useDingbats = FALSE)
	plotErrors(errR, nominalQ=TRUE)
	dev.off()


	# Sample inference and merger of paired-end reads
	mergers <- vector("list", length(sampleNames))
	names(mergers) <- sampleNames
	for(sam in sampleNames) {
		cat("Processing:", sam, "\n")
		derepF <- derepFastq(filtFs[[sam]])
		
		try ( ddF <- dada(derepF, err=errF, multithread=ncores) )
		if (is.null(ddF)){ddF <- dada(derepF, err=errF, multithread=1)}
		derepR <- derepFastq(filtRs[[sam]])
		
		try(ddR <- dada(derepR, err=errR, multithread=ncores))
		if (is.null(ddR)){ddR <- dada(derepR, err=errR, multithread=1)}
		merger <- mergePairs(ddF, derepF, ddR, derepR)
		mergers[[sam]] <- merger
	}
	rm(derepF); rm(derepR)

	# Remove the matched_IDs folder:
	unlink(paste0(path_output,"matched_IDs"),recursive = TRUE)


	# Construct sequence table and remove chimeras
	list_seqtabs[[i]] <- makeSequenceTable(mergers)

	list_seqtabs.nochim[[i]] <- removeBimeraDenovo(list_seqtabs[[i]], method="consensus", multithread=ncores, verbose=FALSE)
} 


##Merge the tables:
if (length(list_seqtabs)>1){
	mergetab <- mergeSequenceTables(tables=list_seqtabs)
	mergetab.nochim <- mergeSequenceTables(tables=list_seqtabs.nochim)
} else {
	mergetab = list_seqtabs[[1]]
	mergetab.nochim = list_seqtabs.nochim[[1]]
}

cat("%",(1-sum(mergetab.nochim)/sum(mergetab))*100,"of the reads were assigned to be chimeric and removed (DADA2)")
uniquesToFasta(getUniques(mergetab.nochim), fout=paste0(path_output,"/uniqueSeqs.fna"), ids=paste0("ASV", seq(length(getUniques(mergetab.nochim)))))
#saveRDS(seqtab, path_output) 
cat("\nDADA2 clustering finished\n")

# If you work with unfiltered data (involving Ns), it will give this error:			  
# Error in dada(drps, err = NULL, errorEstimationFunction = errorEstimationFunction,  : 
# Invalid derep$uniques vector. Sequences must be made up only of A/C/G/T.
