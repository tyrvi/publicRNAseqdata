gtex <- read.delim("GTEx-3-Tissues-N.txt", sep="\t")
gtexcolumns <- colnames(gtex)
# pick a random sample from samples of bladder, prostate and thyroid
random.bladder <- gtexcolumns[sample(grep("bladder", gtexcolumns), size=1)]
random.prostate <- gtexcolumns[sample(grep("prostate", gtexcolumns), size=1)]
random.thyroid <- gtexcolumns[sample(grep("bladder", gtexcolumns), size=1)]

# remove dot in ENSGID
gtex.names <- gtex[,"Name"]
temp_list <- strsplit(as.character(gtex.names), split="\\.")
gtex.names.nodot <- unlist(temp_list)[2*(1:length(gtex.names))-1]
# get columns from samples chosen and place in gtex.fpkms
gtex.fpkms <- data.frame(ENSG_ID=gtex.names.nodot, GTEx_bladder=gtex[,c(random.bladder)], GTEx_prostate=gtex[,c(random.prostate)], GTEx_thyroid=gtex[,c(random.thyroid)])

# write to gtex file
write.table(gtex.fpkms, file="gtexfpkms.txt", quote=FALSE, sep="\t")