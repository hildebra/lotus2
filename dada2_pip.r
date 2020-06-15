
if(!require("dada2",quietly=TRUE,warn.conflicts =FALSE)){
	if (!requireNamespace("BiocManager", quietly = TRUE)){
		install.packages("BiocManager",repos="https://cloud.r-project.org")
	}
	BiocManager::install("dada2", version = "3.11")
}
packageVersion("dada2")

args = commandArgs(trailingOnly=TRUE)


# test if all the arguments are there: 

if (length(args) <= 3) {
  
  stop("All the arguments must be supplied: path path_output seed_num.\n", call.=FALSE)
  
}


# File parsing
pathF=args[1]
path_output=args[2]
seed_num=args[3]
ncores =args[4]


  cat("\n\nStarting DADA2 ASV clustering... please be patient\n\n");


#grep for demultiplexed files on HDD
  fileFs <- list.files(pathF, pattern=".1.fq", full.names = TRUE)
  fileRs <- list.files(pathF, pattern=".2.fq", full.names = TRUE)
#get sample names
  sample.names <- sapply(strsplit(basename(fileFs), "_"), `[`, 1) # 
  sample.namesR <- sapply(strsplit(basename(fileRs), "_"), `[`, 1) #
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(fileFs) <- sample.names
  names(fileRs) <- sample.names

  cat(paste0("Detected ",length(sample.names)," samples\nDetecting error profiles\n"));
  
#start of dada2 learning
  set.seed(seed_num)
# Learn forward error rates
  errF <- learnErrors(fileFs, nbases=1e4, multithread=ncores)
# Learn reverse error rates
  errR <- learnErrors(fileRs, nbases=1e4, multithread=ncores)

# Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(fileFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=ncores)
    derepR <- derepFastq(fileRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=ncores)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
  }
  rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  cat("%",(1-sum(seqtab.nochim)/sum(seqtab))*100,"of the reads were assigned to be chimeric and removed (DADA2)")
  uniquesToFasta(getUniques(seqtab.nochim), fout=paste0(path_output,"/uniqueSeqs.fasta"), ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))
  #saveRDS(seqtab, path_output) 


exit(0);

# If you work with unfiltered data (involving Ns), it will give this error:			  
# Error in dada(drps, err = NULL, errorEstimationFunction = errorEstimationFunction,  : 
# Invalid derep$uniques vector. Sequences must be made up only of A/C/G/T.

#After integrating dada2 into Lotus, you may add a script similar in dada2 (Track reads through the pipeline) to create a table summarizing filtering processes.



#system.time when processing samples individually:
#Timing stopped at: 631.9 57 274

# system.time when polling:
#Timing stopped at: 640 61.24 284.8

