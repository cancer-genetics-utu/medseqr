fastqc.summary <- function(files, outfolder = NULL, out.id=NULL, top=5, pdf.width=12, text.size=9, title=NULL){ ## !! add new variables to body, e.g., allow for adding a keyword to the title
  dat <- read.table(files[1], header = F, sep = "\t", stringsAsFactors = F)
  names <- dat[,2]
  names <- as.vector(sapply(names, gsub, pattern = " ", replacement = "_"))
  res0 <- matrix(ncol=length(names)+1)
  colnames(res0) <- c('File',names)
  res0[1,] <- c(dat[1,3], dat[,1])
  for(i in 2:length(files)){
    dat <- read.table(files[i], header = F, sep = "\t", stringsAsFactors = F)
    res0 <- rbind(res0, c(dat[1,3], dat[,1]))
  }
  if (!is.null(outfolder)) write.csv(res0, file.path(outfolder, paste0(out.id, "allParamPerSamp.csv")), quote=F, row.names=F)
  vlong <- split(as.data.frame(res0), as.data.frame(res0)$File)
  vlong <- lapply(vlong, function(x) {
			  y <- as.data.frame(t(x), stringsAsFactors=F)
			  y$Module <- rownames(y)
			  y$Sample <- sub(".fq.gz|.fastq.gz|.fq|.fastq", "", y[1,1])
			  y <- y[-1,]
			  names(y)[1] <- "Value"
			  y})
  vlong <- do.call("rbind", vlong)
  cbPalette <- c(red="#EE2405", green="#00B868", blue="#005BE4")
  pdf(file.path(outfolder, paste0(out.id, "samples_per_module_barplot.pdf")))
  p <- qplot(data=vlong, x=Module, geom="bar", fill=Value, position="dodge")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=text.size))
  p <- p + ggtitle(paste("Number of samples per QC module", title, sep="\n"))
  p <- p + scale_fill_manual(values=as.character(cbPalette))
  print(p)
  dev.off()
  pdf(file.path(outfolder, paste0(out.id, "module_over_sample_binplot.pdf")), width=pdf.width, height=9)
  p <- qplot(data=vlong, x=Sample, y=Module, geom="bin2d", fill=Value)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=text.size))
  p <- p + ggtitle(paste("FastQC module outcome per sample", title, sep="\n"))
  p <- p + scale_fill_manual(values=as.character(cbPalette))
  print(p)
  dev.off()
  res <- res0
  fpw <- apply(res[, 2:ncol(res)], 2, table)
  for(i in 2:ncol(res)){
    res[which(res[,i] == 'PASS'),i] <- 0
    res[which(res[,i] == 'WARN'),i] <- 1
    res[which(res[,i] == 'FAIL'),i] <- 2
  }
  res <- as.data.frame(res, stringsAsFactors = F)
  for(i in 2:ncol(res)){
    res[,i] <- as.numeric(res[,i])
  }

 sums <- apply(res[,-1], 1, sum)
 tots <- rep(2*length(names), length(sums))
 res <- cbind(res, sums)
 res <- cbind(res, of=tots)
 res <- res[order(res$sums, decreasing = T),]
 max <- max(sums)
 inds <- 1:top


 res1 <- cbind(res[inds,1],
          Total=res[inds,'sums'],
          Of=res[inds,'of'],
          Frac=res[inds,'sums']/res[inds,'of'],
          SD=apply(res[inds,2:13], 1, sd),
          Median=apply(res[inds,-c(1,14,15)], 1, median),
          Mean=apply(res[inds,-c(1,14,15)], 1, mean))
 # Write top worst samples in terms of overall quality
 if (!is.null(outfolder)){
 write.table(res1,
             file = file.path(outfolder, paste0(out.id, 'top_worst_samples.csv')),
             sep = ",",
             row.names = F,
             col.names = T)
 }
 # Top worst variables in terms of quality
 sums <- apply(res[,-c(1, 14,15)], 2, sum)
 max <- max(sums)

 tots <- rep(dim(res)[1]*2, length(sums))

 #resmat <- matrix(nrow=length(sums))

 resdat <- data.frame(Variable = names(sums),
                 Total=sums,
                 Of=tots,
                 Frac=(sums/tots),
                 SD=apply(res[,-c(1,14,15)], 2, sd),
                 Median=apply(res[,-c(1,14,15)], 2, median),
                 Mean=apply(res[,-c(1,14,15)], 2, mean), stringsAsFactors=FALSE)
 resdat <- cbind(resdat, data.frame(PASS=0, WARN=0, FAIL=0))
 for (i in rownames(resdat)) {
	 if(length(grep("PASS", names(fpw[[i]])))>0) resdat[i, "PASS"] <- fpw[[i]][["PASS"]]
	 if(length(grep("WARN", names(fpw[[i]])))>0) resdat[i, "WARN"] <- fpw[[i]][["WARN"]]
	 if(length(grep("FAIL", names(fpw[[i]])))>0) resdat[i, "FAIL"] <- fpw[[i]][["FAIL"]]
 }
 resdat <- resdat[order(resdat[,'Total'], decreasing = T),]

 rownames(resdat) <- NULL
 # violin plot of results per module
 ## first we need to reshape the data frame to be suitable for ggplot2
 vlong <- split(resdat[, c("Variable", "FAIL", "PASS", "WARN")], resdat$Variable)
 vlong <- lapply(vlong, function(x) {
			 y <- as.data.frame(t(x), stringsAsFactors=F)
			 y$Value <- rownames(y)
			 y$Module <- y[1,1]
			 y <- y[-1,]
			 names(y)[1] <- "SampleCount"
             y$SampleCount <- as.numeric(as.character(y$SampleCount))
			 y})
 vlong <- do.call("rbind", vlong)
 ## then plot
 pdf(file.path(outfolder, paste0(out.id, "modules_outcome_violinplot.pdf")), width=pdf.width)
 p <- qplot(data=vlong, x=Value, y=SampleCount, geom="violin", fill=Value)
 p <- p + geom_jitter(height = 0, aes(col=Value))
 p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=text.size))
 p <- p + ggtitle(paste("FastQC module outcome distribution", title, sep="\n"))
 p <- p + scale_fill_manual(values=as.character(cbPalette))
 print(p)
 dev.off()


 # Write top worst samples in terms of overall quality
 if(!is.null(outfolder)){
 write.table(resdat,
             file = file.path(outfolder, paste0(out.id, 'top_worst_variables.csv')),
             sep = ",",
             row.names = F,
             col.names = T)
 }

 return(list(per_sample = res1, per_variable = resdat, all_param_per_sample=res0))
}

