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
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--region_file"), type="character", default=NULL, help="mosdepth regions output file (typically ends with *.regions.bed.gz)", metavar="path"),
                    make_option(c("-s", "--sample_name"), type="character", default=NULL, help="Sample name for plot title. If not provided will be extracted from --region_file", metavar="string"),
                    make_option(c("-o", "--out_file"), type="character", default=NULL, help="Full path to pdf output file. If not provide will be extracted from --region_file", metavar="path"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$region_file)){
    print_help(opt_parser)
    stop("At least one mosdepth region file must be supplied", call.=FALSE)
}

SAMPLE_NAME = opt$sample_name
if (is.null(opt$sample_name)){
    SAMPLE_NAME = gsub('.regions.bed.gz','',basename(opt$region_file))
}

OUT_FILE = opt$out_file
if (is.null(opt$out_file)){
    OUT_FILE = gsub('.gz','.pdf',opt$region_file)
}
if (file.exists(dirname(OUT_FILE)) == FALSE) {
    dir.create(dirname(OUT_FILE),recursive=TRUE)
}

################################################
################################################
## PLOT COVERAGE                              ##
################################################
################################################

## Read in data
dat <- read.csv(gzfile(opt$region_file,'r'),sep="\t", header=FALSE)
colnames(dat) <- c('chrom', 'start','end', 'coverage')
dat$coverage <- dat$coverage + 1

## Coverage plot
plot <- ggplot(dat,aes(x=end,y=coverage)) +
        geom_ribbon(aes(ymin = 0, ymax = coverage), data =) +
        theme_bw() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(trans = log10_trans(),
                           breaks = trans_breaks("log10", function(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x)),
                           expand = c(0, 0)) +
        ylab(bquote('log'[10]~'(Coverage+1)')) +
        xlab("Position (bp)") +
        ggtitle(paste(SAMPLE_NAME,"coverage"))

## Export plot to file
pdf(file=OUT_FILE,height=6,width=12)
print(plot)
dev.off()

################################################
################################################
################################################
################################################
