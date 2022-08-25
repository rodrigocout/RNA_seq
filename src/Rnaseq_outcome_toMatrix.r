#!/usr/bin/RScript
# htseq-combine_all.R
# copy this code to RStudio and adapt file locations to match yours
 #from https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4
# Take 'all' htseq-count results and melt them in to one big dataframe
 
# where are we?
basedir <- "directory where your data are"
setwd(basedir)
 
sampleFiles<-grep('New',list.files(basedir),value=TRUE)
 
cntdir <- paste(basedir, sep="/")
#pat <- "_all_counts.txt"
tophat.all <- list.files(path = cntdir,
	pattern = pat,
	all.files = TRUE,
	recursive = FALSE,
	ignore.case = FALSE,
	include.dirs = FALSE)
 
# we choose the 'all' series
myfiles <- sampleFiles
DT <- list()
 
# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
	infile = paste(cntdir, myfiles[i], sep = "/")
	DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
	cnts <- gsub("(.*)_all_counts.txt", "\\1", myfiles[i])
	colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}
 
# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]
 
# inspect
head(data)
 
# ID SRR1039508
# 1 ENSG00000000003        667
# 2 ENSG00000000005          0
# 3 ENSG00000000419        430
# 4 ENSG00000000457        256
# 5 ENSG00000000460         56
# 6 ENSG00000000938          0
 
# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
	y <- DT[[myfiles[i]]]
	z <- merge(data, y, by = c("ID"))
	data <- z
}
 
 #remove duplicates
 df <- data[!duplicated(data),]
  # ID column becomes rownames
rownames(df) <- df$ID
df <- df[,-1]
 
## add total counts per sample
data <- rbind(df, tot.counts=colSums(df))
 
# inspect and look at the top row names!
head(data)
 
# SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512 SRR1039513
# __alignment_not_unique    6029848    5535220    6917679    5736110    7174795    3829823
# __ambiguous                609790     567580     563265     582869     711348     449909
# __no_feature              2048439    1871251    2186546    2036064    2564408    1387796
# __not_aligned                   0          0          0          0          0          0
# __too_low_aQual                 0          0          0          0          0          0
# ENSG00000000003               667        434        721        444        862        401
# SRR1039514 SRR1039515 SRR1039516 SRR1039517 SRR1039518 SRR1039519
# __alignment_not_unique    8865980    5533574    7349997    7829016    7357178    5435207
# __ambiguous               1021469     607559     688966     890454     809693     575376
# __no_feature              3464421    2230297    2530061    3020635    2678110    1990138
# __not_aligned                   0          0          0          0          0          0
# __too_low_aQual                 0          0          0          0          0          0
# ENSG00000000003               958        853       1133       1050       1087       1104
# SRR1039520 SRR1039521 SRR1039522 SRR1039523
# __alignment_not_unique    5019316    5159442    6941580    7078894
# __ambiguous                554414     613151     709211     818199
# __no_feature              1939908    1931496    2641504    2825142
# __not_aligned                   0          0          0          0
# __too_low_aQual                 0          0          0          0
# ENSG00000000003               750        562       1075        826
 
tail(data)
 
# SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512 SRR1039513 SRR1039514
# ENSG00000273487          3          4          4          6          0          9          5
# ENSG00000273488          5          3          4          3          5          3         10
# ENSG00000273489          0          0          1          0          0          2          0
# ENSG00000273492          0          0          0          0          1          0          2
# ENSG00000273493          0          0          0          0          0          0          0
# tot.counts        26792115   24519985   27416297   25783685   33081391   19381676   45865162
# SRR1039515 SRR1039516 SRR1039517 SRR1039518 SRR1039519 SRR1039520 SRR1039521
# ENSG00000273487          4          6          4          5          3          3          8
# ENSG00000273488          4          3         10          6          6          6         11
# ENSG00000273489          2          1          0          2          1          0          0
# ENSG00000273492          1          0          0          1          0          0          0
# ENSG00000273493          0          0          0          0          0          0          0
# tot.counts        27862053   32316085   39563928   35371127   25763223   24653719   26874854
# SRR1039522 SRR1039523
# ENSG00000273487          6          9
# ENSG00000273488          7          6
# ENSG00000273489          3          1
# ENSG00000273492          1          0
# ENSG00000273493          0          1
# tot.counts        32888531   35869342
 
####################################
# take summary rows to a new table
# ( not starting with ENS with invert=TRUE )
 
# transpose table for readability
data.all.summary <- data[grep("^ENS", rownames(data), perl=TRUE, invert=TRUE), ]
 
# review
data.all.summary
 
# SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512 SRR1039513
# __alignment_not_unique    6029848    5535220    6917679    5736110    7174795    3829823
# __ambiguous                609790     567580     563265     582869     711348     449909
# __no_feature              2048439    1871251    2186546    2036064    2564408    1387796
# __not_aligned                   0          0          0          0          0          0
# __too_low_aQual                 0          0          0          0          0          0
# tot.counts               26792115   24519985   27416297   25783685   33081391   19381676
# SRR1039514 SRR1039515 SRR1039516 SRR1039517 SRR1039518 SRR1039519
# __alignment_not_unique    8865980    5533574    7349997    7829016    7357178    5435207
# __ambiguous               1021469     607559     688966     890454     809693     575376
# __no_feature              3464421    2230297    2530061    3020635    2678110    1990138
# __not_aligned                   0          0          0          0          0          0
# __too_low_aQual                 0          0          0          0          0          0
# tot.counts               45865162   27862053   32316085   39563928   35371127   25763223
# SRR1039520 SRR1039521 SRR1039522 SRR1039523
# __alignment_not_unique    5019316    5159442    6941580    7078894
# __ambiguous                554414     613151     709211     818199
# __no_feature              1939908    1931496    2641504    2825142
# __not_aligned                   0          0          0          0
# __too_low_aQual                 0          0          0          0
# tot.counts               24653719   26874854   32888531   35869342
 
# transpose table
t(data.all.summary)
 
# write summary to file
write.csv(data.all.summary, file = "chr22-htseq_counts_all-summary.csv")
 
####################################
# take all data rows to a new table
 
data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]
 
# final merged table
head(data.all, 3)
 
# write data to file
write.csv(data.all, file = "chr22-htseq_counts_all.csv")
 
# cleanup intermediate objects
rm(y, z, i, DT)