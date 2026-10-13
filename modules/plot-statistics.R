#!/usr/bin/env Rscript
#
# plot-statistics.R - graphs of the final results of Prokfinder (called by summarize-final-results.sh)
# Uses base R only, so no extra package needs to be installed.
#
# Usage: Rscript plot-statistics.R <genome_metadata.tsv> <tables_dir> <graphs_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot-statistics.R <genome_metadata.tsv> <tables_dir> <graphs_dir>")
}
metadata_file <- args[1]
tables_dir <- args[2]
graphs_dir <- args[3]

meta <- read.delim(metadata_file, sep = "\t", quote = "", na.strings = c("NULL", "NA", ""),
                   check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")

# Opens a PNG in the graphs folder
open_png <- function(name) {
  png(file.path(graphs_dir, paste0(name, ".png")), width = 1800, height = 1200, res = 200)
}

# 1. Selection funnel
funnel_file <- file.path(tables_dir, "selection-funnel.tsv")
if (file.exists(funnel_file)) {
  funnel <- read.delim(funnel_file, stringsAsFactors = FALSE)
  funnel <- funnel[!is.na(suppressWarnings(as.numeric(funnel$GENOMES))), ]
  if (nrow(funnel) > 0) {
    open_png("selection-funnel")
    par(mar = c(5, 18, 3, 2))
    counts <- rev(as.numeric(funnel$GENOMES))
    bp <- barplot(counts, names.arg = rev(sub("^[0-9]+_", "", funnel$STEP)), horiz = TRUE, las = 1,
                  col = "#4C78A8", border = NA, xlab = "Genomes", main = "Genomes left after each step",
                  xlim = c(0, max(counts) * 1.15))
    text(counts, bp, labels = counts, pos = 4, cex = 0.8)
    invisible(dev.off())
  }
}

# 2. Distribution of the quality metrics (all assessed genomes; the QC thresholds are drawn as lines)
metrics <- list(
  list(col = "CHECKM2_COMPLETENESS", label = "CheckM2 completeness (%)", line = as.numeric(Sys.getenv("QC_MIN_COMPLETENESS", "75"))),
  list(col = "CHECKM2_CONTAMINATION", label = "CheckM2 contamination (%)", line = as.numeric(Sys.getenv("QC_MAX_CONTAMINATION", "5"))),
  list(col = "CLEAN_N50", label = "N50 (bp, log10)", line = as.numeric(Sys.getenv("QC_MIN_N50", "5000")), log = TRUE),
  list(col = "CLEAN_NUM_CONTIGS", label = "Number of contigs (log10)", line = as.numeric(Sys.getenv("QC_MAX_CONTIGS", "1000")), log = TRUE),
  list(col = "CLEAN_GENOME_SIZE_BP", label = "Genome size (Mb)", scale = 1e6),
  list(col = "CLEAN_GC_PERCENT", label = "GC content (%)")
)
available <- Filter(function(m) m$col %in% names(meta) && sum(!is.na(suppressWarnings(as.numeric(meta[[m$col]])))) > 1, metrics)
if (length(available) > 0) {
  open_png("quality-metrics")
  par(mfrow = c(2, ceiling(length(available) / 2)), mar = c(4, 4, 2, 1))
  for (m in available) {
    values <- suppressWarnings(as.numeric(meta[[m$col]]))
    values <- values[!is.na(values)]
    threshold <- if (is.null(m$line)) NA else m$line
    if (isTRUE(m$log)) {
      values <- log10(values[values > 0])
      threshold <- log10(threshold)
    }
    if (!is.null(m$scale)) values <- values / m$scale
    hist(values, breaks = 40, col = "#72B7B2", border = "white", main = m$label, xlab = "", ylab = "Genomes")
    if (!is.na(threshold)) abline(v = threshold, col = "#E45756", lwd = 2, lty = 2)
  }
  invisible(dev.off())
}

# 3. Genomes per GTDB phylum and genus (top 25)
for (rank in c("phylum", "genus")) {
  table_file <- file.path(tables_dir, paste0("genomes-per-gtdb-", rank, ".tsv"))
  if (!file.exists(table_file)) next
  tab <- read.delim(table_file, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
  if (nrow(tab) == 0) next
  tab <- head(tab[order(-tab[[2]]), ], 25)
  open_png(paste0("genomes-per-gtdb-", rank))
  par(mar = c(5, 14, 3, 2))
  counts <- t(as.matrix(tab[rev(seq_len(nrow(tab))), c(3, 4)]))
  barplot(counts, beside = TRUE, horiz = TRUE, las = 1, names.arg = rev(tab[[1]]), cex.names = 0.7,
          col = c("#4C78A8", "#F58518"), border = NA, xlab = "Genomes",
          main = paste("High-quality and representative genomes per GTDB", rank, "(top 25)"))
  legend("bottomright", legend = c("High quality", "Representative"), fill = c("#4C78A8", "#F58518"), bty = "n")
  invisible(dev.off())
}

cat("Graphs saved to:", graphs_dir, "\n")
