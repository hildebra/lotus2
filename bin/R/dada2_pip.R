
if(!require("dada2",quietly=TRUE,warn.conflicts =FALSE)){
	if (!requireNamespace("BiocManager", quietly = TRUE)){
		install.packages("BiocManager",repos="https://cloud.r-project.org")
	}
	BiocManager::install("dada2")
}
#packageVersion("dada2")

derepFastqRead= function (fls, n = 1e+06, verbose = FALSE, qualityType = "Auto") 
{
	suppressMessages(library(ShortRead,quietly = TRUE))
    if (!is.character(fls)) {
        stop("File paths must be provided in character format.")
    }
    if (!all(file.exists(fls))) {
        stop("Not all provided files exist.")
    }
    rval <- list()
	fl <- fls[[1]]#was loop before..
	if (verbose) {
		message("Dereplicating sequence entries in Fastq file: ", fl, appendLF = TRUE)
	}
	rds = readFastq(fl,qualityType=qualityType)
	
	IDspl = lapply(lapply(as.character(id(rds)),strsplit,";size="),"[[",1)
	IDs = as.character(id(rds)) #currently I need the full names
	#unlist(lapply(IDspl,"[[",1))
	derepCounts = as.numeric(unlist(lapply(lapply(IDspl,"[[",2),strsplit,";")))
	derepQuals = as(quality(rds), "matrix")
	names(derepCounts) = as.character(sread(rds))
	rownames(derepQuals) = names(derepCounts)
	
	ord <- order(derepCounts, decreasing = TRUE)
	derepCounts <- derepCounts[ord]
	IDs = IDs[ord]
	derepQuals <- derepQuals[ord, , drop = FALSE]
	derepMap=NULL
	#pseudo val, completely useless but fullfilling dada2 data requirements
	#derepMap=seq(sum(derepCounts))
	#derepMap <- match(derepMap, ord)
	if (verbose) {
		message("Encountered ", length(derepCounts), 
			" unique sequences from ", sum(derepCounts), 
			" total sequences read.")
	}
	derepO <- list(uniques = derepCounts, quals = derepQuals, map = derepMap, IDs=IDs)
	derepO <- as(derepO, "derep")
    return(derepO)
}




#ARGS parsing
#args=c("/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S/tmpFiles/demultiplexed/","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//tmpFiles//","0","12","/hpc-home/hildebra/dev/lotus/maps/AngeTest1.16S.sm.map","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S/tmpFiles/derep.fas")
#args=c("/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//demultiplexed/","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//demultiplexed/","0","12","/hpc-home/hildebra/dev/lotus/maps/Lucas_16S__map2.txt")
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


#-------------- prep the input (demultiplexed) file paths --------------

subset_map=array("A",dim(mapping)[1])
if ("SequencingRun" %in% colnames(mapping)){ #explicitly defined groups of sequencing runs
	subset_map = mapping[,"SequencingRun"]
#	for (i in 1:length(unique(mapping[,"SequencingRun"]))){
#		subset_map[[i]]=subset(mapping,SequencingRun==row.names(table(mapping[,"SequencingRun"]))[i])
#	}
} else if ("fastqFile" %in% colnames(mapping) ){ #1)check if different BCs 2)check if different folders to impute seq runs
	tfqFs = table(mapping[,"fastqFile"])
	hasSlash=grepl("/",names(tfqFs))
	if (any(tfqFs > 1)){ #demultiplexing included
		subset_map = mapping[,"fastqFile"]
		write(paste("Found ",length(tfqFs)," unique fastqFile's, assumming each represents a sequencing run\n\n",sep=""),stderr())
	} else if (any(hasSlash)) {
		leastS=function (x){lx=length(x);if (lx==1){""} else {paste(x[1:(lx-1)],collapse="/")}}
		fastqP = unlist(lapply(strsplit(mapping[,"fastqFile"],"/"),leastS))
		tfqFs=table(fastqP)
		if (length(tfqFs) > 1){
			write(paste("Found ",length(tfqFs)," unique path's to fastq's, assumming each path represents a sequencing run\n\n",sep=""),stderr())
			subset_map = fastqP
		}
	} else {
		write("No sample run information, samples are considered as sequenced in the same run\n\n",stderr())
	}
} else {
	write("No sample run information, samples are considered as sequenced in the same run\n\n",stderr())
}



# File parsing
listF=list();listR=list()
maxReg=500
tSuSe = table(subset_map)
for (i in names(tSuSe)){
	#forward read error pattern
	idx = which(subset_map == i)
	listF[[i]] = paste0(pathF, paste0(mapping[idx,"#SampleID"],".1.fq",sep="") )
	listR[[i]] = paste0(pathF, paste0(mapping[idx,"#SampleID"],".2.fq",sep="") )
	for (XX in 1:length(listF[[i]])){
		if (!file.exists(listF[[i]][XX])){
			listF[[i]][XX]= paste0(listF[[i]][XX],".gz")
			listR[[i]][XX]= paste0(listR[[i]][XX],".gz")
		}
		if (!file.exists(listF[[i]][XX])){
			cat(paste("Could not find expected file", listF[[i]][XX] ,"\n"))
			listF[[i]][XX] = ""
			listR[[i]][XX] = ""
		}
	}
	#next;
}
#listF[[names(tSuSe)[1]]][1] = ""


#-------------- pooled clustering --------------
for (i in 1:length(listF)){
	idx = which(subset_map == names(tSuSe)[i])
	sampleNames = mapping[idx,"#SampleID"]
	cat(paste0("Detected ",length(sampleNames)," samples in batch ",i,"/",length(listF), " .. Computing error profiles\n"));
	#start of dada2 learning
	set.seed(seed_num)
	filtFs <- file.path(listF[[i]]);	names(filtFs) <- sampleNames;	
	filtFs = filtFs[file.exists(filtFs)]
	# Learn forward error rates
	cat(paste0("Learning error profiles for the forward reads:\n"));
	#forward read error rates
	
	try( errF <- learnErrors(filtFs, multithread=ncores))#dada2 is too instable
	if (is.null(errF)) {errF <- learnErrors(filtFs, multithread=1)}
	# Learn reverse error rates

	# Save the plots of error profiles, for a sanity check:
	pdf(paste0(path_output,"/dada2_p",i,"_errF.pdf"),useDingbats = FALSE)
	plotErrors(errF, nominalQ=TRUE);	dev.off()
	
	#cat(paste0("Learning error profiles for the reverse reads:\n"));	
	filtRs <- file.path(listR[[i]]);names(filtRs) <- sampleNames 
	if (0){
		# second read
		try( errR <- learnErrors(filtRs, multithread=ncores))
		if (is.null(errR)) {errR <- learnErrors(filtRs, multithread=1)}
		pdf(paste0(path_output,"/dada2_p",i,"_errR.pdf"),useDingbats = FALSE)
		plotErrors(errR, nominalQ=TRUE);	dev.off()
	}
	break;#no reason currently to go further..
}

##########################################################
#new implementation relying on sdm derep
##########################################################
if (length(args)>5){ 
	derepFile=args[6]
	tDerep = derepFastqRead(derepFile)
	tdada=dada(tDerep,err=errF,multithread=ncores)
	ASVseq=as.character(names(tdada$denoised))
	ASVab=tdada$clustering$abundance

	mapA=tdada$map #links to tDerep
	tars = table(mapA)
	#nTs = names(tars)
	Dids = tDerep$IDs
	ASVname="otu"
	#create .uc file from dada2
	fileUC=file(paste0(path_output,"dada2.uc"),open ="wt")
	fileFNA=file(paste0(path_output,"uniqueSeqs.fna"),open ="wt")
	for (i in 1:length(ASVseq)){
		idx = which(mapA == i)
		for (j in 1:length(idx)){
			if (j==1){
				cat(paste0(Dids[idx[j]],"\t",ASVname,i,"\t*\n"),file=fileUC)
				cat(paste0(">",Dids[idx[j]],"\n",ASVseq[[i]],"\n"),file=fileFNA) #ASV
			} else {
				cat(paste0(Dids[idx[j]],"\tmatch\tdqt=1;top=",Dids[idx[1]],"(99%);\n"),file=fileUC)
			}
		}
	}
	close(fileUC); close (fileFNA);
	#write unique sequences
	cat(paste0("Found ",length(ASVseq)," ASVs, summing to ",sum(ASVab)," reads (dada2)\n"));
	q("no");
}

cat("Wrong branch\n")
q('no')

cat("\n\nStarting DADA2 ASV clustering... please be patient\n\n");

list_seqtabs=list()
list_seqtabs.nochim=list()
mergers = list()

i=1
for (i in 1:length(listF)){
	idx = which(subset_map == names(tSuSe)[i])
	sampleNames = mapping[idx,"#SampleID"]
	cat(paste0("Detected ",length(sampleNames)," samples in batch ",i,"/",length(listF), " .. Computing error profiles\n"));
	#start of dada2 learning
	set.seed(seed_num)
	filtFs <- file.path(listF[[i]]);	filtRs <- file.path(listR[[i]])
	names(filtFs) <- sampleNames;	names(filtRs) <- sampleNames 
	# Learn forward error rates
	cat(paste0("Learning error profiles for the forward reads:\n"));
	#forward read error rates
	try( errF <- learnErrors(filtFs, multithread=ncores))#dada2 is too instable
	if (is.null(errF)) {errF <- learnErrors(filtFs, multithread=1)}
	# Learn reverse error rates
	cat(paste0("Learning error profiles for the reverse reads:\n"));	
	# second read
	try( errR <- learnErrors(filtRs, multithread=ncores))
	if (is.null(errR)) {errR <- learnErrors(filtRs, multithread=1)}

	# Save the plots of error profiles, for a sanity check:
	pdf(paste0(path_output,"/dada2_p",i,"_errF.pdf"),useDingbats = FALSE)
	plotErrors(errF, nominalQ=TRUE);	dev.off()
	pdf(paste0(path_output,"/dada2_p",i,"_errR.pdf"),useDingbats = FALSE)
	plotErrors(errR, nominalQ=TRUE);	dev.off()


	# Sample inference and merger of paired-end reads
	#mergers = vector("list",length(sampleNames))
	#names(mergers) <- sampleNames
	
	for(sam in sampleNames) {
		cat("Processing:", sam, "\n")
		derepF <- derepFastq(filtFs[[sam]],n=5e5,verbose=FALSE)
		try ( ddF <- dada(derepF, err=errF, multithread=ncores) )
		if (is.null(ddF)){ddF <- dada(derepF, err=errF, multithread=1)}
		#dada(..., HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32) #IT, 454

		
		derepR <- derepFastq(filtRs[[sam]],n=5e5,verbose=FALSE)
		try(ddR <- dada(derepR, err=errR, multithread=ncores))
		if (is.null(ddR)){ddR <- dada(derepR, err=errR, multithread=1)}
		merger <- mergePairs(ddF, derepF, ddR, derepR,trimOverhang=TRUE,minOverlap=6)
		mergers[[sam]] <- merger
		rm(derepF); rm(derepR)
	}

	# Remove the matched_IDs folder:
	#unlink(paste0(path_output,"matched_IDs"),recursive = TRUE)


	# Construct sequence table and remove chimeras
#	list_seqtabs[[i]] <- makeSequenceTable(mergers)
#	list_seqtabs.nochim[[i]] <- removeBimeraDenovo(list_seqtabs[[i]], method="consensus", multithread=ncores, verbose=FALSE)
}

##Merge the tables:
#if (length(list_seqtabs)>1){
#	mergetab <- mergeSequenceTables(tables=list_seqtabs)
#	mergetab.nochim <- mergeSequenceTables(tables=list_seqtabs.nochim)
#} else {
#	mergetab = list_seqtabs[[1]]
#	mergetab.nochim = list_seqtabs.nochim[[1]]
#}
mergetab = makeSequenceTable(mergers)
mergetab.nochim <- removeBimeraDenovo(mergetab, method="consensus", multithread=ncores, verbose=FALSE)



cat(paste0(round((1-sum(mergetab.nochim)/sum(mergetab))*100,2),"% of reads (",dim(mergetab)[2]-dim(mergetab.nochim)[2]," of ",dim(mergetab)[2]," ASVs) were chimeric and will be removed (DADA2)"))
uniquesToFasta(getUniques(mergetab.nochim), fout=paste0(path_output,"/uniqueSeqs.fna"), ids=paste0("ASV", seq(length(getUniques(mergetab.nochim)))))
#saveRDS(seqtab, path_output) 
cat("\nDADA2 clustering finished\n")

# If you work with unfiltered data (involving Ns), it will give this error:			  
# Error in dada(drps, err = NULL, errorEstimationFunction = errorEstimationFunction,  : 
# Invalid derep$uniques vector. Sequences must be made up only of A/C/G/T.
