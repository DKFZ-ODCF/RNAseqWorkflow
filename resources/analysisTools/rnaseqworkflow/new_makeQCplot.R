## new_makeQCplot.R
# author: Simon Steiger
# modified by Jeongbin Park

require("RColorBrewer") #installed
require("gplots") #installed
require("gtools")
require("reshape")
require("ggplot2")

makeQC_plot <- function(counts, sample="SingleCell", png_path=getwd(), chrom_counts=NULL) {
    print(sample)

    output.dir <- file.path(png_path,sample)
    dir.create(output.dir)

    counts <- counts
    if(length(grep("ERCC",rownames(chrom_counts)))==0){ercc <- FALSE} else {ercc<-TRUE}

    if(is.null(chrom_counts)){
        chrom_counts <- matrix(ncol=ncol(counts),nrow=25,0)
        colnames(chrom_counts) <- colnames(counts)
    }

    ## ERCC Calculation
    if((sum(colnames(counts) %in% colnames(chrom_counts)) - ncol(counts))>0){
        print(paste("Removed",sum(colnames(counts) %in% colnames(chrom_counts)) - ncol(counts),"cells missing in the chromosome_reads"))
    } else {
        print(paste(ncol(counts),"cells in QC for plotting"))
    }

    counts <- counts[,colnames(counts) %in% colnames(chrom_counts)]
    counts_all <- counts

    # Removing Ensembl_ID decimal (xxx.4)
    rownames(counts) <- gsub("\\.\\d+", "\\1", rownames(counts))

    require(org.Hs.eg.db)
    ids <- unlist(lapply( mget( rownames(counts),
                                org.Hs.egENSEMBL2EG,
                                ifnotfound = NA),
                          function(x) x[1] ))

    ids <- na.omit(ids)

    genechr <- unlist(lapply( mget( ids,
                                    org.Hs.egCHR,
                                    ifnotfound = NA),
                              function(x) x[1]))


    Gensymbs <- unlist(lapply( mget( ids,
                                     org.Hs.egSYMBOL,
                                     ifnotfound = NA),
                               function(x) x[1]))

    rids <- Gensymbs; names(rids) <- names(ids)

    # Switching EnsemblID to GeneSymbol

    # Linking available gene symbol to rowname / ensemblID
    nasymbs <- unlist(lapply( rownames(counts),
                              function(x) rids[x]))

    # rownames(counts) <- nasymbs

    counts <- counts[!is.na(rownames(counts)),]

    print("Finished Gene_Symb conversion")

    ColSS <- as.data.frame(colnames(counts))
    colnames(ColSS) <- "id"

    if(ercc){
        ColSS$ercc <- as.data.frame(colSums(chrom_counts[grep("ERCC",rownames(chrom_counts)),]),row.names=NULL)
        ColSS$percc <- colSums(chrom_counts[grep("ERCC",rownames(chrom_counts)),])/(colSums(chrom_counts[grep("ERCC",rownames(chrom_counts)),])+colSums(counts_all))
    }


    ColSS$Mtreads <- colSums(counts[genechr=="MT",])/colSums(counts_all)
    ColSS$Rsum <- colSums(counts_all)
    ColSS$GenCount <- colSums(counts>0)
    ColSS$ngenicReads <- (colSums(counts_all)-colSums(counts))/colSums(counts_all)
    ColSS <- ColSS[with(ColSS,mixedorder(id)),]

    write.table(ColSS,file.path(output.dir,paste0(sample,"_QCMetrics.csv")),dec=",",quote=FALSE,sep=";")

    melted <- melt(ColSS, id.vars=c("id"))
    #groupsS <- unlist(lapply(strsplit(melted$id,"_"),function(x) paste(x[c(1,2,3)],collapse="_")))

    #samplecol <- brewer.pal(length(levels(factor(groupsS))),"Set3")

    melted$id <- as.factor(melted$id)

    melted$variable <- factor(melted$variable)

    capitalize <- function(string) {
    substr(string, 1, 1) <- toupper(substr(string, 1, 1))
    string
    }
    if(ercc){
    conservation_status <- c(
    erccSum = "ERCC Counts",
    percc = "ERCC %",
    Mtreads = "MitoReads %",
    Rsum = "ReadSum",
    ngenicReads = "Non-genic Reads %",
    GenCount = "GeneCount"
    )
    } else {
    conservation_status <- c(
    Mtreads = "MitoReads %",
    Rsum = "ReadSum",
    ngenicReads = "Non-genic Reads %",
    GenCount = "GeneCount"
    )}

    print("Plotting")


    p1 <- ggplot(melted, aes(x=id, y=value, fill = variable)) +
          geom_bar(stat="identity") +
          facet_grid(variable ~ ., scales="free_y",labeller = labeller(.default=capitalize,variable=conservation_status)) +
          ggtitle(paste("Quality Control for Sample",sample)) +
          scale_x_discrete(limits=levels(melted$id)) +
          theme(legend.position="bottom"
                ,axis.title.y = element_blank()
                ,text = element_text(size=9)
                ,axis.text.y = element_text(vjust=1,size=12),
                axis.title.x = element_blank()
                ,axis.text.x = element_text(angle=-90,vjust=1,size=4))+
          guides(fill=FALSE)

    png(file=file.path(output.dir,paste(sample,"QCplot.png",sep="_")),
      width=7*300,
      height=7*300,
      res=300)
    print(p1)
    dev.off()
}