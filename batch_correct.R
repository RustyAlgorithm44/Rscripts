# This is the function, the function call is below

ComBat_data <- function(expr, batch.info, batch = "Batch", NameString = "")
{
  print ("===========================Batch Effects Adjustment using ComBat=====================")
  #matching IDs for ComBat
  if (is.character(batch.info[,1])){
    match.id <- match (as.character(colnames(expr)), batch.info[,1])
    batch.id <- batch.info[match.id, 2]
  } else if (is.numeric(batch.info[,1])){
    match.id <- match (as.numeric(colnames(expr)), batch.info[,1])
    batch.id <- batch.info[match.id, 2]
  }
  print("Performing batch correction using ComBat...")
  batch_corrected <- sva::ComBat(dat=expr,
                                 batch = batch.id,
                                 mod=NULL,
                                 par.prior=TRUE,
                                 prior.plots=FALSE)
  date <- as.character(format(Sys.Date(), "%Y%m%d"))
  if(NameString==""){
    outfile  <- outFile <- paste0(date, "_data_", "batch_corrected_", batch, ".txt")
  } else {
    outFile <- paste0(date, "_data_", NameString, "_batch_corrected_", batch, ".txt")
  }
  write.table(batch_corrected, outFile, quote=FALSE, sep="\t")
  print(paste("Batch corrected data written to file... ", outFile, sep=""))
  return (batch_corrected)
}

#m1 is the expression file, with samples in the column and genes in the row
#m2 is the batches file, with Samples in 1st column and Batch in 2nd column (containing values like 1,2,3)
m1 <- read.csv("merged_exp.csv")
rownames(m1) <- m1[,1]
m1 <- m1[,-1]
m2 <- read.csv("batches.csv")
ComBat_data(m1,m2,batch = "Batch")