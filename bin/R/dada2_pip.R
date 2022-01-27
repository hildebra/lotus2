
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
	existFiles=file.exists(fls)
    if (!all(existFiles)) {
        stop(paste("Not all provided files exist:",fls[!existFiles]))
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
combineDada2 = function (samples, orderBy = "abundance") 
{
	if (class(samples) %in% c("dada", "derep", "data.frame")) {
		samples <- list(samples)
	}
	if (!is.list(samples)) {
		stop("Requires a list of samples.")
	}
	unqs <- lapply(samples, getUniques) #gets $denoised seqs from dada2
	unqsqs <- unique(do.call(c, lapply(unqs, names))) #compares them
	sums = 0;
	for (x in names(samples)){sums = sums + sum((samples[[x]]$clustering$abundance))}
	rval <- matrix(0L, nrow = length(unqs), ncol = length(unqsqs))
	colnames(rval) <- unqsqs
	for (i in seq_along(unqs)) {
		rval[i, match(names(unqs[[i]]), colnames(rval))] <- unqs[[i]]
	}
	if (!is.null(names(unqs))) {
		rownames(rval) <- names(unqs)
	}
	if (!is.null(orderBy)) {
		if (orderBy == "abundance") {
			rval <- rval[, order(colSums(rval), decreasing = TRUE), 
				drop = FALSE]
		}
		else if (orderBy == "nsamples") {
			rval <- rval[, order(colSums(rval > 0), decreasing = TRUE), 
				drop = FALSE]
		}
	}
	return(list(rval=rval,sums=sums))
}

isBimeraDenovo2 = function (unqs, minFoldParentOverAbundance = 2, minParentAbundance = 8, 
    allowOneOff = FALSE, minOneOffParentDistance = 4, maxShift = 16, 
    multithread = FALSE, verbose = FALSE) 
{
	stop("outdate isBimeraDenovo2")
    unqs.int <- getUniques(unqs, silence = TRUE)
    abunds <- unname(unqs.int)
    seqs <- names(unqs.int)
    seqs.input <- getSequences(unqs)
    rm(unqs)
    gc(verbose = FALSE)
    if (is.logical(multithread)) {
        if (multithread == TRUE) {
            mc.cores <- getOption("mc.cores", detectCores())
        }
    }
    else if (is.numeric(multithread)) {
        mc.cores <- multithread
        multithread <- TRUE
    }
    else {
        warning("Invalid multithread parameter. Running as a single thread.")
        multithread <- FALSE
    }
    loopFun <- function(i, unqs.loop, minFoldParentOverAbundance, 
        minParentAbundance, allowOneOff, minOneOffParentDistance, 
        maxShift) {
        sq <- names(unqs.loop)[[i]]
        abund <- unqs.loop[[i]]
        pars <- names(unqs.loop)[(unqs.loop > (minFoldParentOverAbundance * 
            abund) & unqs.loop > minParentAbundance)]
        if (length(pars) < 2) {
            return(FALSE)
        }
        else {
            isBimera(sq, pars, allowOneOff = allowOneOff, minOneOffParentDistance = minOneOffParentDistance, 
                maxShift = maxShift)
        }
    }
    if (multithread) {
        mc.indices <- sample(seq_along(unqs.int), length(unqs.int))
        bims <- mclapply(mc.indices, loopFun, unqs.loop = unqs.int, 
            allowOneOff = allowOneOff, minFoldParentOverAbundance = minFoldParentOverAbundance, 
            minParentAbundance = minParentAbundance, minOneOffParentDistance = minOneOffParentDistance, 
            maxShift = maxShift, mc.cores = mc.cores)
        bims <- bims[order(mc.indices)]
    }
    else {
        bims <- lapply(seq_along(unqs.int), loopFun, unqs.loop = unqs.int, 
            allowOneOff = allowOneOff, minFoldParentOverAbundance = minFoldParentOverAbundance, 
            minParentAbundance = minParentAbundance, minOneOffParentDistance = minOneOffParentDistance, 
            maxShift = maxShift)
    }
    bims <- unlist(bims)
    bims.out <- seqs.input %in% seqs[bims]
    names(bims.out) <- seqs.input
    if (verbose) 
        message("Identified ", sum(bims.out), " bimeras out of ", 
            length(bims.out), " input sequences.")
    return(bims.out)
}






#ARGS parsing
#args=c("/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//tmpFiles//demultiplexed/","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//tmpFiles/","0","2","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//primary/in.map","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//tmpFiles//derep.fas")
#args=c("/hpc-home/hildebra/grp/data/results/lotus/Anh//tmpFiles//demultiplexed/","/hpc-home/hildebra/grp/data/results/lotus/Anh//tmpFiles/","0","14","/hpc-home/hildebra/grp/data/results/lotus/Anh//primary/in.map","/hpc-home/hildebra/grp/data/results/lotus/Anh//tmpFiles//derep.fas")
#args=c("/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/data/test_data/lotus2/Arabi_output/Arabidopsis_ITSd2/tmpFiles//demultiplexed/","/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/data/test_data/lotus2/Arabi_output/Arabidopsis_ITSd2/tmpFiles/","0","12","/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/data/test_data/lotus2/Arabi_output/Arabidopsis_ITSd2/primary/in.map","/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/data/test_data/lotus2/Arabi_output/Arabidopsis_ITSd2/tmpFiles//derep.fas")

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
bp4error = 5e7

#read map to get file locations, sampleRuns etc
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
listF=list();listR=list();listM=list();
maxReg=500
tSuSe = table(subset_map)
for (i in names(tSuSe)){
	#forward read error pattern
	idx = which(subset_map == i)
	listF[[i]] <- file.path(pathF, paste0(mapping[idx, "#SampleID"], ".1.fq"))
	listR[[i]] <- file.path(pathF, paste0(mapping[idx, "#SampleID"], ".2.fq"))
	listM[[i]] <- file.path(pathF, paste0(mapping[idx, "#SampleID"], ".merg.fq"))

	for (XX in 1:length(listF[[i]])){
		if (!file.exists(listF[[i]][XX])){
			listF[[i]][XX]= paste0(listF[[i]][XX],".gz")
			listR[[i]][XX]= paste0(listR[[i]][XX],".gz")
		}
		if (!file.exists(listF[[i]][XX])){
			listF[[i]][XX]= paste0(listF[[i]][XX],".gz")
			listR[[i]][XX]= paste0(listR[[i]][XX],".gz")
			if (!file.exists(listF[[i]][XX])){
				#cat(paste("Could not find expected file", listF[[i]][XX] ,"\n"))
				listF[[i]][XX] = ""
				listR[[i]][XX] = ""
			}
		}
		if (!file.exists(listM[[i]][XX])){
			listM[[i]][XX]= paste0(listM[[i]][XX],".gz")
			if (!file.exists(listM[[i]][XX])){
				listM[[i]][XX] = ""
			}
		}
	}
	#last resort..
	if (all(!file.exists(listF[[i]]))){#single file case
		listF[[i]] <- file.path(pathF, paste0(mapping[idx, "#SampleID"], ".fq"))
		for (XX in 1:length(listF[[i]])){
			if (!file.exists(listF[[i]][XX])){
				listF[[i]][XX]= paste0(listF[[i]][XX],".gz")
			}
		}
	}
	#next;
}

#listF[[names(tSuSe)[1]]][1] = ""

#do we deal with merged data? needs to be clear early on
mergedData=FALSE
if (length(args)>5){ 
	derepFile1=args[6]
	xd=strsplit(derepFile1,"\\.")[[1]]
	if (xd[length(xd)-1] == "merg"){#working on merged data..
		mergedData=TRUE
	}
}


#double check all relevant files
for (i in names(tSuSe)){
	if (!mergedData){
		if (all(!file.exists(listF[[i]]))) {#length(listM[[i]]) == 0){
			stop(paste0("Can't find files for block ",i," expected files such as \n",listF[[i]][0],"\n\nAborting dada2 run\n"))
		}
	} else {
		if (all(!file.exists(listM[[i]]))) {#length(listM[[i]]) == 0){
			stop(paste0("Can't find files for block ",i," expected files such as \n",listM[[i]][0],"\n\nAborting dada2 run\n"))
		}
	}
}

errorsF = errorsM = list()
#-------------- pooled clustering --------------
cnt = 0
i = sort(names(tSuSe))[1]
for (i in sort(names(tSuSe))){
	idx = which(subset_map == i)
	cnt = cnt +1
	sampleNames = mapping[idx,"#SampleID"]
	cat(paste0("Detected ",length(sampleNames)," samples in batch ",i," (",cnt,"/",length(listF), ") .. Computing error profiles\n"));
	#start of dada2 learning
	set.seed(seed_num)
	
	# Learn forward error rates
	#forward read error rates
	if (mergedData){
		defMergeFile = listM[[i]][0]
		filtMs <- file.path(listM[[i]]);	names(filtMs) <- sampleNames;filtMs = filtMs[file.exists(filtMs)]
		if (length(filtMs) == 0){
			stop(paste0("Can't find derep for block ",i," expected\n",defMergeFile,"\n\nAborting dada2 run\n"))
		}
		cat(paste0("Learning error profiles for merged reads\n"));
		#DEBUG - to find problematic input files
		#for (j in 1:length(filtMs)){learnErrors(filtMs[j], nbases = bp4error, multithread=ncores)}
		try( errM <- learnErrors(filtMs, nbases = bp4error, multithread=ncores))#dada2 is too instable
		if (is.null(errM)) {errM <- learnErrors(filtMs, nbases = bp4error, multithread=1)}
		# Save the plots of error profiles, for a sanity check:
		pdf(file.path(path_output, paste0("dada2_p", i, "_errM.pdf")), useDingbats = FALSE)
		suppressWarnings(print(plotErrors(errM, nominalQ=TRUE)));
		dev.off()
		errorsM[[i]] = errM
	} 
	if (!mergedData || 1){
		filtFs <- file.path(listF[[i]]);	names(filtFs) <- sampleNames;filtFs = filtFs[file.exists(filtFs)]
		cat(paste0("Learning error profiles for the forward reads:\n"));
		try( errF <- learnErrors(filtFs, nbases = bp4error, multithread=ncores))#dada2 is too instable
		if (is.null(errF)) {errF <- learnErrors(filtFs, nbases = bp4error, multithread=1)}
		# Save the plots of error profiles, for a sanity check:
		pdf(file.path(path_output, paste0("dada2_p", i, "_errF.pdf")), useDingbats = FALSE)
		suppressWarnings(print(plotErrors(errF, nominalQ=TRUE)))
		dev.off() 
		errorsF[[i]] = errF
	}
	
	#cat(paste0("Learning error profiles for the reverse reads:\n"));	
	filtRs <- file.path(listR[[i]]);names(filtRs) <- sampleNames 
	if (0){
		# second read
		try( errR <- learnErrors(filtRs, nbases = bp4error, multithread=ncores))
		if (is.null(errR)) {errR <- learnErrors(filtRs, nbases = bp4error, multithread=1)}
		pdf(paste0(path_output,"/dada2_p",i,"_errR.pdf"),useDingbats = FALSE)
		suppressWarnings(print(plotErrors(errR, nominalQ=TRUE)));	dev.off()
	}
	#break;#no reason currently to go further..
}

##########################################################
#new implementation relying on sdm derep
##########################################################


if (length(args)>5){ 
	derepFile1=args[6]
	xd=strsplit(derepFile1,"\\.")[[1]]
	if (mergedData){
		derePref = paste(xd[1:(length(xd)-2)],collapse=".")
		derePost = paste(xd[(length(xd)-1):length(xd)],collapse=".")
	} else {
		derePref = paste(xd[1:(length(xd)-1)],collapse=".")
		derePost = xd[length(xd)]
	}
	ldada=list();tDerep=list()
	cnt=1
	for (i in sort(names(tSuSe))){
		derepFile = paste0(derePref,".",i,".",derePost)
		cat("Running dada on derep fastq ",cnt,"/",length(tSuSe),"(",derepFile,")\n")
		tDerep[[i]] = derepFastqRead(derepFile)
		if (length(tDerep[[i]]$IDs) == 0){
			ldada[[i]]=NULL
		} else {
		if (mergedData){
			locErr = errorsM[[i]];
		} else {
			locErr = errorsF[[i]];
		}
			ldada[[i]]=dada(tDerep[[i]],err=locErr,multithread=ncores)
		}
		cnt = cnt+1
	}
	for (kk in 1:length(ldada)){
		#remove chimeras in subsets
		num_prev = length(ldada[[kk]]$denoised)
		sum_prev = sum(ldada[[kk]]$clustering$abundance)
		#pooled  consensus  
		#isBimeraDenovo2
		isBimera = isBimeraDenovo(ldada[[kk]],multithread=ncores,verbose=TRUE) 
		#ldada[[kk]]$denoised = ldada[[kk]]$denoised 
		
		sum_aft = sum(ldada[[kk]]$clustering$abundance[unname(isBimera)])
		num_aft = sum(isBimera)
		cat(paste0("Removed ",num_prev-num_aft," chimeric ASVs(",sum_prev-sum_aft," read counts)\n"))
		
		ldada[[kk]]$clustering = ldada[[kk]]$clustering[!isBimera,]
		ldada[[kk]]$denoised  = ldada[[kk]]$denoised [!isBimera]
		
		
		tmp =names(ldada[[kk]]$denoised)
		ldada[[kk]]$denoised = seq(length(ldada[[kk]]$denoised))
		names(ldada[[kk]]$denoised) = tmp

	}
	#tdada=ldada[[1]]
	cat("Finished initial per-subset dada2 clustering\n")
	cat("Merging dada2 clusters between subsets..\n")
	tdada2 = combineDada2(ldada)
	ASVseq=as.character(colnames(tdada2$rval))
	#names(ASVseqFnd) = ASVseq
	ASVab = tdada2$sums
	
	#old implementation relying on single dada2 object
	#ASVseq=as.character(names(tdada$denoised))
	#ASVab=tdada$clustering$abundance

	#mapA=tdada$map #links to tDerep
	#tars = table(mapA)
	#nTs = names(tars)
	#Dids = tDerep$IDs
	ASVname="otu"
	tothits=0
	ASVseqFnd=array("",length(ASVseq))
	#create .uc file from dada2
	cat("Converting dada2 clustering to deterministic read mapping..\n")
	fileUC <- file(file.path(path_output, "dada2.uc"), open = "wt")
	fileFNA <- file(file.path(path_output, "uniqueSeqs.fna"), open = "wt")
	for (i in 1:length(ASVseq)){#OTUs & their original sequencing need to be written out in a block ??
		for (kk in 1:length(ldada)){
			nASV = ldada[[kk]]$denoised[ASVseq[i]]
			if (is.na(nASV)){next}
			idx = which(ldada[[kk]]$map == nASV)
			Dids = tDerep[[kk]]$IDs
			for (j in 1:length(idx)){
				if (ASVseqFnd[i] == ""){
					cat(paste0(Dids[idx[j]],"\t",ASVname,i,"\t*\n"),file=fileUC)
					cat(paste0(">",Dids[idx[j]],"\n",ASVseq[[i]],"\n"),file=fileFNA) #ASV
					ASVseqFnd[i]=Dids[idx[j]]
				} else {
					cat(paste0(Dids[idx[j]],"\tmatch\tdqt=1;top=",ASVseqFnd[i],"(99%);\n"),file=fileUC)
				}
				tothits=tothits+1
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
	try( errF <- learnErrors(filtFs,nbases = bp4error,  multithread=ncores))#dada2 is too instable
	if (is.null(errF)) {errF <- learnErrors(filtFs, nbases = bp4error, multithread=1)}
	# Learn reverse error rates
	cat(paste0("Learning error profiles for the reverse reads:\n"));	
	# second read
	try( errR <- learnErrors(filtRs, multithread=ncores))
	if (is.null(errR)) {errR <- learnErrors(filtRs, multithread=1)}

	# Save the plots of error profiles, for a sanity check:
	pdf(file.path(path_output, paste0("dada2_p", i, "_errF.pdf")), useDingbats = FALSE)
	plotErrors(errF, nominalQ=TRUE);	dev.off()
	pdf(file.path(path_output, paste0("dada2_p", i, "_errR.pdf")), useDingbats = FALSE)
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
	#unlink(file.path(path_output, "matched_IDs"), recursive = TRUE)


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
uniquesToFasta(getUniques(mergetab.nochim), fout=file.path(path_output, "uniqueSeqs.fna"), ids=paste0("ASV", seq(length(getUniques(mergetab.nochim)))))
#saveRDS(seqtab, path_output) 
cat("\nDADA2 clustering finished\n")

# If you work with unfiltered data (involving Ns), it will give this error:			  
# Error in dada(drps, err = NULL, errorEstimationFunction = errorEstimationFunction,  : 
# Invalid derep$uniques vector. Sequences must be made up only of A/C/G/T.
