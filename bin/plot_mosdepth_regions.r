#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(viridis)
library(tidyverse)

################################################
################################################
## VALIDATE COMMAND-LINE PARAMETERS           ##
################################################
################################################

option_list <- list(make_option(c("-i", "--input_files"), type="character", default=NULL, help="Comma-separated list of mosdepth regions output file (typically end in *.regions.bed.gz)", metavar="input_files"),
                    make_option(c("-s", "--input_suffix"), type="character", default='.regions.bed.gz', help="Portion of filename after sample name to trim for plot title e.g. '.regions.bed.gz' if 'SAMPLE1.regions.bed.gz'", metavar="input_suffix"),
                    make_option(c("-o", "--output_dir"), type="character", default='./', help="Output directory", metavar="path"),
                    make_option(c("-p", "--output_suffix"), type="character", default='regions', help="Output suffix", metavar="output_suffix"))

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
    dat <- rbind(dat, cbind(read.delim(input_file, header=FALSE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)[,-6], sample, stringsAsFactors=F))
}

## Reformat table
if (ncol(dat) == 6) {
    colnames(dat) <- c('chrom', 'start','end', 'region', 'coverage', 'sample')
    dat$region <- factor(dat$region, levels=unique(dat$region[order(dat$start)]))
} else {
    colnames(dat) <- c('chrom', 'start','end', 'coverage', 'sample')
}
dat$sample <- factor(dat$sample, levels=sort(unique(dat$sample)))

## Write merged coverage data for all samples to file
outfile <- paste(OUTDIR,"all_samples.",OUTSUFFIX,".coverage.tsv", sep='')
write.table(dat, file=outfile, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

################################################
################################################
## PER-SAMPLE COVERAGE PLOTS                  ##
################################################
################################################

for (sample in unique(dat$sample)) {
    sample_dat <- dat[dat$sample == sample,]
    outfile <- paste(OUTDIR,sample,".",OUTSUFFIX,".coverage.tsv", sep='')
    write.table(sample_dat,file=outfile, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
    sample_dat$coverage <- sample_dat$coverage + 1

    if (ncol(sample_dat) == 6) {
        plot <- ggplot(sample_dat,aes(x=region,y=coverage)) +
                geom_bar(stat="identity", fill="#D55E00", width=0.6) +
                theme_bw() +
                theme(plot.title=element_text(size=10),
                      axis.text.x=element_text(size=10),
                      axis.text.y=element_text(size=8)) +
                coord_flip() +
                scale_x_discrete(expand=c(0, 0)) +
                scale_y_continuous(trans=log10_trans(),
                                   breaks=trans_breaks('log10', function(x) 10^x),
                                   labels=trans_format('log10', math_format(10^.x)),
                                   expand=c(0, 0)) +
                ylab(bquote('log'[10]~'(Coverage+1)')) +
                xlab('Amplicon') +
                ggtitle(paste(sample,'per amplicon coverage'))

          outfile <- paste(OUTDIR,sample,".",OUTSUFFIX,".coverage.pdf", sep='')
          ggsave(file=outfile, plot, width=16, height=3+(0.3*length(unique(sample_dat$region))), units="cm")
    } else {
        plot <- ggplot(sample_dat,aes(x=end,y=coverage)) +
                geom_ribbon(aes(ymin=0, ymax=coverage), fill="#D55E00", data=) +
                theme_bw() +
                scale_x_continuous(expand=c(0, 0)) +
                scale_y_continuous(trans=log10_trans(),
                                   breaks=trans_breaks('log10', function(x) 10^x),
                                   labels=trans_format('log10', math_format(10^.x)),
                                   expand=c(0, 0)) +
                ylab(bquote('log'[10]~'(Coverage+1)')) +
                xlab('Position (bp)') +
                ggtitle(paste(sample,'coverage'))

        outfile <- paste(OUTDIR,sample,".",OUTSUFFIX,".coverage.pdf", sep='')
        ggsave(file=outfile, plot, width=12, height=6, units="in")
    }
}

################################################
################################################
## REGION-BASED HEATMAP ACROSS ALL SAMPLES    ##
################################################
################################################

if (ncol(dat) == 6 && length(INPUT_FILES) > 1) {
    mat <- spread(dat[,c("sample", "region", "coverage")], sample, coverage, fill=NA, convert=FALSE)
    rownames(mat) <- mat[,1]
    mat <- as.matrix(log10(mat[,-1] + 1))
    heatmap <- Heatmap(mat,
                       name                = "log10(Coverage+1)",
                       cluster_rows        = FALSE,
                       cluster_columns     = FALSE,
                       show_row_names      = TRUE,
                       show_column_names   = TRUE,
                       column_names_side   = "bottom",
                       rect_gp             = gpar(col="white", lwd=1),
                       show_heatmap_legend = TRUE,
                       row_names_gp        = gpar(fontsize=6),
                       column_names_gp     = gpar(fontsize=6),
                       height              = unit(5, "mm")*nrow(mat),
                       width               = unit(5, "mm")*ncol(mat),
                       col                 = viridis(50))

    ## Size of heatmaps scaled based on matrix dimensions: https://jokergoo.github.io/ComplexHeatmap-reference/book/other-tricks.html#set-the-same-cell-size-for-different-heatmaps-with-different-dimensions
    width = 0.1969*ncol(mat) + (2*1.3150)
    height = 0.1969*nrow(mat) + 1.3150
    outfile <- paste(OUTDIR,"all_samples.",OUTSUFFIX,".heatmap.pdf", sep='')
    pdf(file=outfile, width=width, height=height)
    draw(heatmap)
    dev.off()
}

################################################
################################################
################################################
################################################