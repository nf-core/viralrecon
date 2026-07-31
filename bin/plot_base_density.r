#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(scales)
library(reshape2)
library(Biostrings)

################################################
################################################
## VALIDATE COMMAND-LINE PARAMETERS           ##
################################################
################################################

option_list <- list(make_option(c("-i", "--fasta_files"), type="character", default=NULL, help="Comma-separated list of fasta files", metavar="fasta_files"),
                    make_option(c("-s", "--prefixes"), type="character", default=NULL, help="Comma-separated list of prefixes associated with fasta files to add to plots. Must be unique and in same order as fasta file input.", metavar="prefixes"),
                    make_option(c("-o", "--output_dir"), type="character", default='./', help="Output directory", metavar="path"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## Check input files
INPUT_FILES <- unique(unlist(strsplit(opt$fasta_files,",")))
if (length(INPUT_FILES) == 0) {
    print_help(opt_parser)
    stop("At least one input file must be supplied", call.=FALSE)
}
if (!all(file.exists(INPUT_FILES))) {
    stop(paste("The following input files don't exist:",paste(INPUT_FILES[!file.exists(INPUT_FILES)], sep='', collapse=' '), sep=' '), call.=FALSE)
}

## Check prefixes for input files
PREFIXES <- basename(INPUT_FILES)
if (!is.null(opt$prefixes)){
    PREFIXES <- unique(unlist(strsplit(opt$prefixes,",")))
    if (length(INPUT_FILES) != length(PREFIXES)) {
        print_help(opt_parser)
        stop("Please provide a unique prefix for each fasta file.", call.=FALSE)
    }
}

## Check the output directory has a trailing slash, if not add one
OUTDIR <- opt$output_dir
if (tail(strsplit(OUTDIR,"")[[1]],1)!="/") {
    OUTDIR <- paste(OUTDIR,"/",sep='')
}
## Create the directory if it doesn't already exist.
if (!file.exists(OUTDIR)) {
    dir.create(OUTDIR, recursive=TRUE)
}

################################################
################################################
## READ IN DATA                               ##
################################################
################################################

dat <- NULL
for (input_file in INPUT_FILES) {
    dat <- c(dat,readDNAStringSet(input_file)[1])
}

################################################
################################################
## PLOTS                                      ##
################################################
################################################

bases_std <-  c("A","C","T","G")
base_cols <-  c("A" = "#009E73",
                "C" = "#0072B2",
                "T" = "#D55E00",
                "G" = "#000000",
                "N" = "#E69F00",
                "X" = "#999999")

for (idx in 1:length(dat)) {

    ## Table of base counts
    base_seq <- strsplit(toString(dat[[idx]]), "")[[1]]
    base_tab <- data.frame(table(base_seq), stringsAsFactors=FALSE)
    colnames(base_tab) <- c("base","freq")
    rownames(base_tab) <- base_tab$base
    for (base in 1:length(bases_std)) {
        if (!any(base_tab$base %in% bases_std[base])) {
            base_tab <- rbind(base_tab,c(bases_std[base],0))
        }
    }
    base_tab$perc <- 100 *base_tab$freq / sum(base_tab$freq)
    base_tab <- base_tab[order(base_tab$base, decreasing=FALSE),]
    base_tab <- rbind(base_tab[c(bases_std, "N"),], base_tab[!rownames(base_tab) %in% c(bases_std, "N"),])
    base_tab$base <- factor(base_tab$base, levels=rownames(base_tab))
    outfile <- paste(OUTDIR, PREFIXES[idx], ".base_counts.tsv", sep='')
    write.table(base_tab, file=outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

    ## Barplot of base frequencies
    barplot <-  ggplot(base_tab, aes(x=base,y=perc)) +
                geom_bar(stat="identity") +
                theme_classic() +
                scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100)) +
                ylab("% Observed") +
                xlab("Base") +
                ggtitle(PREFIXES[idx])
    outfile <- paste(OUTDIR, PREFIXES[idx], ".base_counts.pdf", sep='')
    ggsave(file=outfile, barplot, width=12, height=10, units="cm")

    ## Create a data frame of base coverage
    bases <- unique(c(bases_std,"N",unique(base_seq)))
    base_dat <- data.frame(sample=names(dat[[idx]])[1], position=1:length(base_seq), stringsAsFactors=FALSE)
    for (base in 1:length(bases)) {
        base_dat[,bases[base]] <- as.numeric(base_seq==bases[base])
    }

    ## Stretches of N's
    N_rle <- Rle(base_dat[,"N"])
    N_dat <- data.frame(start=cumsum(runLength(N_rle))[runValue(N_rle)==1], width=runLength(N_rle)[runValue(N_rle)==1])
    outfile <- paste(OUTDIR, PREFIXES[idx], ".N_run.tsv", sep='')
    write.table(N_dat, file=outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

    ## Running mean of bp density for standard bases
    run_k <- 1001
    run_dat <- base_dat[,c("sample", "position", bases_std)]
    for (base in bases_std) {
        run_dat[,base] <- as.numeric(runmean(Rle(base_dat[,base]), k=run_k, endrule="constant"))
    }
    run_dat <- melt(run_dat, c(1,2))
    colnames(run_dat)[3] <- "base"
    run_dat$position <- run_dat$position/1000
    lineplot <- ggplot(run_dat,aes(x=position, y=value, colour=base)) +
                geom_line() +
                theme_classic() +
                theme(panel.border=element_rect(colour="black", fill=NA, size=1)) +
                scale_y_continuous(breaks=c(0,0.25,0.50,0.75,1)) +
                xlab("Position (Kb)") +
                ylab(paste("Base density (running mean k=",run_k,")", sep='')) +
                ggtitle(PREFIXES[idx]) +
                scale_colour_manual(values=base_cols)
    outfile <- paste(OUTDIR, PREFIXES[idx], ".ACTG_density.pdf", sep='')
    ggsave(file=outfile, lineplot, width=18, height=10, units="cm")

    ## Single base density plots, nucleotide resolution.
    bases_other <- bases[!bases %in% bases_std]
    for (obase in bases_other) {
        plot_dat  <- base_dat[,c("sample", "position", obase)]
        colnames(plot_dat)[3] <- "base"
        plot_col <- ifelse(obase=="N", base_cols[["N"]], base_cols[["X"]])
        lineplot <- ggplot(plot_dat, aes(x=position/1000, y=base)) +
                    geom_line(colour=plot_col) +
                    theme_classic() +
                    theme(legend.position="none", panel.border=element_rect(colour="black", fill=NA, size=1)) +
                    scale_y_continuous(breaks=c(0,1), labels=c(0,1)) +
                    xlab("Position (Kb)") +
                    ylab(paste(obase,"density", sep=' ')) +
                    ggtitle(PREFIXES[idx])
        outfile <- paste(OUTDIR, PREFIXES[idx], ".", obase, "_density.pdf", sep='')
        ggsave(file=outfile, lineplot, width=18, height=10, units="cm")
    }
}
