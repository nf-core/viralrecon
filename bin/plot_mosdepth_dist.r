#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(scales)

################################################
################################################
## VALIDATE COMMAND-LINE PARAMETERS           ##
################################################
################################################

option_list <- list(make_option(c("-i", "--input_files"), type="character", default=NULL, help="Comma-separated list of mosdepth regions output file (typically end in *.mosdepth.global.dist.txt)", metavar="input_files"),
                    make_option(c("-s", "--input_suffix"), type="character", default='.mosdepth.global.dist.txt', help="Portion of filename after sample name to trim for plot title e.g. '.mosdepth.global.dist.txt' if 'SAMPLE1.mosdepth.global.dist.txt'", metavar="input_suffix"),
                    make_option(c("-o", "--output_dir"), type="character", default='./', help="Output directory", metavar="path"),
                    make_option(c("-p", "--output_suffix"), type="character", default='global.dist', help="Output suffix", metavar="output_suffix"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## Check input files
INPUT_FILES <- unique(unlist(strsplit(opt$input_files,",")))
if (length(INPUT_FILES) == 0) {
    print_help(opt_parser)
    stop("At least one input file must be supplied", call.=FALSE)
}
if (!all(file.exists(INPUT_FILES))) {
    stop(paste("The following input files don't exist:",paste(INPUT_FILES[!file.exists(INPUT_FILES)], sep='', collapse=' '), sep=' '), call.=FALSE)
}

## Check the output directory has a trailing slash, if not add one
OUTDIR <- opt$output_dir
if (tail(strsplit(OUTDIR,"")[[1]],1)!="/") {
    OUTDIR <- paste(OUTDIR,"/",sep='')
}
## Create the directory if it doesn't already exist.
if (!file.exists(OUTDIR)) {
  dir.create(OUTDIR,recursive=TRUE)
}

OUTSUFFIX <- trimws(opt$output_suffix, "both", whitespace = "\\.")

################################################
################################################
## READ IN DATA                               ##
################################################
################################################

## Read in data
dat <- NULL
for (input_file in INPUT_FILES) {
    sample = gsub(opt$input_suffix,'',basename(input_file))
    dat <- rbind(dat, cbind(read.delim(input_file, header=FALSE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)[,-4], sample, stringsAsFactors=F))
}
colnames(dat) <- c('chrom', 'coverage', 'frequency', 'sample')
dat <- dat[which(dat$chrom == 'total'),][,2:ncol(dat)]

################################################
################################################
## PER-SAMPLE COVERAGE PLOTS                  ##
################################################
################################################

for (sample in unique(dat$sample)) {
    sample_dat <- dat[dat$sample == sample,]
    plot <- ggplot(sample_dat,aes(x=coverage,y=frequency)) +
            geom_line(stat="identity") +
            theme_bw() +
            scale_x_continuous(expand=c(0, 0)) +
            scale_y_continuous(limits=c(0,1),
                               breaks=seq(0,1,0.2),
                               labels=seq(0,1,0.2),
                               expand=c(0, 0)) +
            ylab('Proportion of genome at coverage') +
            xlab('Coverage') +
            ggtitle(paste(sample,' genome coverage'))

      outfile <- paste(OUTDIR,sample,".",OUTSUFFIX,".coverage.pdf", sep='')
      ggsave(file=outfile, plot, height=4, width=8, units="in")
}

################################################
################################################
## COVERAGE PLOT ACROSS ALL SAMPLES           ##
################################################
################################################

if (length(INPUT_FILES) > 1) {
    plot <- ggplot(dat,aes(x=coverage,y=frequency,colour=sample)) +
            geom_line(stat="identity") +
            theme_bw() +
            scale_x_continuous(expand=c(0, 0)) +
            scale_y_continuous(limits=c(0,1),
                               breaks=seq(0,1,0.2),
                               labels=seq(0,1,0.2),
                               expand=c(0, 0)) +
            ylab('Proportion of genome at coverage') +
            xlab('Coverage') +
            ggtitle(paste('All samples genome coverage'))

    outfile <- paste(OUTDIR,"all_samples.",OUTSUFFIX,".coverage.pdf", sep='')
    ggsave(file=outfile, plot, height=6, width=12, units="in")
}

################################################
################################################
################################################
################################################
