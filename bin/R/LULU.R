#uses LULU to clean OTU tables 

if(!require("dplyr",quietly=TRUE,warn.conflicts =FALSE)){
	install.packages("dplyr")
	require("dplyr",warn.conflicts =FALSE)
}
#library(lulu)


#lulu too complicated to install via devtools, use this one instead
lulu = function (otutable, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95) 
{
    suppressWarnings(require(dplyr,warn.conflicts =FALSE))
    start.time <- Sys.time()
    colnames(matchlist) <- c("OTUid", "hit", "match")
    matchlist = matchlist[which(matchlist$hit != "*"), ]
    matchlist = matchlist[which(matchlist$hit != matchlist$OTUid), ]
	otutable=otutable[rowSums(otutable)>0,,drop=FALSE]
    statistics_table <- otutable[, 0]
    statistics_table$total <- rowSums(otutable)
    statistics_table$spread <- rowSums(otutable > 0)
    statistics_table <- statistics_table[with(statistics_table, order(spread, total, decreasing = TRUE)), ]
    otutable <- otutable[match(row.names(statistics_table), row.names(otutable)),]
    statistics_table$parent_id <- "NA"
    log_con <- file(file.path(logD, paste0("lulu.log_", format(start.time, "%Y%m%d_%H%M%S"))), open = "a")
	tarProg=0.1
    for (line in seq(1:nrow(statistics_table))) {
		if (line/nrow(statistics_table)>tarProg){
			print(paste0("progress: ", round(((line/nrow(statistics_table)) * 100), 0), "%"))
			tarProg = tarProg+0.1
		}
        potential_parent_id <- row.names(otutable)[line]
        cat(paste0("\n", "####processing: ", potential_parent_id, " #####"), file = log_con)
        daughter_samples <- otutable[line, ]
        hits <- matchlist[which(matchlist$OTUid == potential_parent_id & matchlist$match > minimum_match), "hit"]
        cat(paste0("\n", "---hits: ", hits), file = log_con)
        last_relevant_entry <- sum(statistics_table$spread >= statistics_table$spread[line])
        potential_parents <- which(row.names(otutable)[1:last_relevant_entry] %in% hits)
        cat(paste0("\n", "---potential parent: ", row.names(statistics_table)[potential_parents]), file = log_con)
        success <- FALSE
        if (length(potential_parents) > 0) {
            for (line2 in potential_parents) {
                cat(paste0("\n", "------checking: ", row.names(statistics_table)[line2]), file = log_con)
                if (!success) {
                  relative_cooccurence <- sum((daughter_samples[otutable[line2, ] > 0]) > 0)/sum(daughter_samples > 0)
                  cat(paste0("\n", "------relative cooccurence: ", relative_cooccurence), file = log_con)
                  if (relative_cooccurence >= minimum_relative_cooccurence) {
                    cat(paste0(" which is sufficient!"), file = log_con)
                    if (minimum_ratio_type == "avg") {
                      relative_abundance <- mean(otutable[line2, ][daughter_samples > 0]/daughter_samples[daughter_samples > 
                        0])
                      cat(paste0("\n", "------mean avg abundance: ", relative_abundance), file = log_con)
                    } else {
                      relative_abundance <- min(otutable[line2, ][daughter_samples > 0]/daughter_samples[daughter_samples > 
                        0])
                      cat(paste0("\n", "------min avg abundance: ", relative_abundance), file = log_con)
                    }
                    if (relative_abundance > minimum_ratio) {
                      cat(paste0(" which is OK!"), file = log_con)
                      if (line2 < line) {
                        statistics_table$parent_id[line] <- statistics_table[row.names(otutable)[line2], "parent_id"]
                        cat(paste0("\n", "SETTING ", potential_parent_id, " to be an ERROR of ", (statistics_table[row.names(otutable)[line2], "parent_id"]), "\n"), file = log_con)
                      }
                      else {
                        statistics_table$parent_id[line] <- row.names(otutable)[line2]
                        cat(paste0("\n", "SETTING ", potential_parent_id, " to be an ERROR of ", (row.names(otutable)[line2]), "\n"), file = log_con)
                      }
                      success <- TRUE
                    }
                  }
                }
            }
        }
        if (!success) {
            statistics_table$parent_id[line] <- row.names(statistics_table)[line]
            cat(paste0("\n", "No parent found!", "\n"), file = log_con)
        }
    }
    close(log_con)
    total_abundances <- rowSums(otutable)
    curation_table <- cbind(nOTUid = statistics_table$parent_id, otutable)
    statistics_table$curated <- "merged"
    curate_index <- row.names(statistics_table) == statistics_table$parent_id
    statistics_table$curated[curate_index] <- "parent"
    statistics_table <- transform(statistics_table, rank = ave(total, FUN = function(x) rank(-x, ties.method = "first")))
    curation_table <- as.data.frame(curation_table %>% group_by(nOTUid) %>% summarise_all(funs(sum)))
    row.names(curation_table) <- as.character(curation_table$nOTUid)
    curation_table <- curation_table[, -1]
    curated_otus <- names(table(statistics_table$parent_id))
    curated_count <- length(curated_otus)
    discarded_otus <- setdiff(row.names(statistics_table), curated_otus)
    discarded_count <- length(discarded_otus)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    result <- list(curated_table = curation_table, curated_count = curated_count, curated_otus = curated_otus, discarded_count = discarded_count, discarded_otus = discarded_otus, runtime = time.taken, minimum_match = minimum_match, minimum_relative_cooccurence = minimum_relative_cooccurence, otu_map = statistics_table, original_table = otutable)
    return(result)
}
library(compiler)
lulu=cmpfun(lulu, options = NULL)


#args=c("lulu_match_list.txt","OTU.txt","logs")
#args=c("/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//tmpFiles//lulu_match_list.txt","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S//OTU.txt","/hpc-home/hildebra/grp/data/results/lotus/Angela/Test1.16S/")
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3){stop("Not enough commandline args!")}

matchL = args[1]
otuF = args[2]
logD = args[3]

info = file.info(matchL)
if (info$size == 0){
	q("no");
}

matchs = read.table(matchL,header=FALSE,as.is=TRUE)
otuM = read.table(otuF,header=TRUE,as.is=TRUE,row.names=1)
if (dim(otuM)[1] <= 1 || dim(otuM)[2] <= 1){
	q("no");
}

lulu <- lulu(otuM, matchs)
write.table(lulu$curated_table,quote =FALSE,sep="\t",file=otuF,col.names = NA)
write.table(lulu$discarded_otus,quote =FALSE,sep="\t",file=paste0(matchL,".rm"),col.names =FALSE,row.names =FALSE)
#save the lulu object for log purposes
save(lulu,file=paste0(logD,"LULU.Rdata"));
