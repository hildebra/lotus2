
if(!require("dada2",quietly=TRUE,warn.conflicts =FALSE)){
	if (!requireNamespace("BiocManager", quietly = TRUE)){
		install.packages("BiocManager",repos="https://cloud.r-project.org")
	}
	BiocManager::install("dada2", version = "3.11")
}
#packageVersion("dada2")

#ARGS parsing
args = commandArgs(trailingOnly=TRUE)
# test if all the arguments are there: 
if (length(args) <= 3) {
	stop("All the arguments must be supplied: pathF path_output seed_num ncores.\n", call.=FALSE)
}
pathF=args[1]
path_output=args[2]
seed_num=args[3]
ncores =args[4]


cat("\n\nStarting DADA2 ASV clustering... please be patient\n\n");

# File parsing

#grep for demultiplexed files on HDD
fileFs <- list.files(pathF, pattern=".1.fq", full.names = TRUE)
fileRs <- list.files(pathF, pattern=".2.fq", full.names = TRUE)
#get sample names
sampleNames <- sapply(strsplit(basename(fileFs), ".1.fq"),`[`, 1) # 
sampleNamesR <- sapply(strsplit(basename(fileRs), ".2.fq"),`[`, 1) #
if(!identical(sampleNames, sampleNamesR)) stop("Forward and reverse files do not match:",sampleNames,sampleNamesR)
names(fileFs) <- sampleNames
names(fileRs) <- sampleNames


print(sampleNames)

cat(paste0("Detected ",length(sampleNames)," samples\nDetecting error profiles\n"));

#start of dada2 learning
set.seed(seed_num)

filtFs <- file.path(path_output, "matched_IDs", paste0(sampleNames, "_F_matched.fastq.gz"))
filtRs <- file.path(path_output, "matched_IDs", paste0(sampleNames, "_R_matched.fastq.gz"))
names(filtFs) <- sampleNames
names(filtRs) <- sampleNames

#Not doing any filtering but just to match the IDs of reverse and forward reads:
filterAndTrim(fwd=fileFs, filt=filtFs, rev=fileRs, filt.rev=filtRs,matchIDs=TRUE)

# Learn forward error rates
cat(paste0("Learning error profiles for the forward reads:\n"));
errF <- learnErrors(filtFs, nbases=1e4, multithread=ncores)

cat(paste0("Learning error profiles for the reverse reads:\n"));
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e4, multithread=ncores)

# Save the plots of error profiles, for a sanity check:
pdf(paste0(path_output,"dada2_p_errF.pdf"),useDingbats = FALSE)
plotErrors(errF, nominalQ=TRUE)
dev.off()
pdf(paste0(path_output,"dada2_p_errR.pdf"),useDingbats = FALSE)
plotErrors(errR, nominalQ=TRUE)
dev.off()
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sampleNames))
names(mergers) <- sampleNames
for(sam in sampleNames) {
	cat("Processing:", sam, "\n")
	derepF <- derepFastq(filtFs[[sam]])
	ddF <- dada(derepF, err=errF, multithread=ncores)
	derepR <- derepFastq(filtRs[[sam]])
	ddR <- dada(derepR, err=errR, multithread=ncores)
	merger <- mergePairs(ddF, derepF, ddR, derepR)
	mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Remove the matched_IDs folder:
unlink(paste0(path_output,"matched_IDs"),recursive = TRUE)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=ncores, verbose=TRUE)
cat("%",(1-sum(seqtab.nochim)/sum(seqtab))*100,"of the reads were assigned to be chimeric and removed (DADA2)")
uniquesToFasta(getUniques(seqtab.nochim), fout=paste0(path_output,"/uniqueSeqs.fna"), ids=paste0("ASV", seq(length(getUniques(seqtab.nochim)))))
#saveRDS(seqtab, path_output) 
cat("\nDADA2 clustering finished\n")

# If you work with unfiltered data (involving Ns), it will give this error:			  
# Error in dada(drps, err = NULL, errorEstimationFunction = errorEstimationFunction,  : 
# Invalid derep$uniques vector. Sequences must be made up only of A/C/G/T.

