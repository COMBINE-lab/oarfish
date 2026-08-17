#!/usr/bin/env Rscript
# Production DE over oarfish-converted salmon dirs: swish (fishpond) and
# edgeR v4 with inferential-uncertainty moderation (catchSalmon), plus naive
# edgeR without moderation as the comparator.
# Usage: run_de.R <dir-with-sample-subdirs A1..A4,B1..B4> <out-prefix>
suppressMessages({
  library(tximport)
  library(fishpond)
  library(edgeR)
  library(SummarizedExperiment)
})
args <- commandArgs(trailingOnly = TRUE)
base <- args[1]
outp <- args[2]
samples <- c(paste0("A", 1:4), paste0("B", 1:4))
cond <- factor(c(rep("A", 4), rep("B", 4)))
paths <- file.path(base, samples)
files <- file.path(paths, "quant.sf")
stopifnot(all(file.exists(files)))

## ---- swish -------------------------------------------------------------
txi <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = FALSE)
infList <- txi$infReps
nrep <- ncol(infList[[1]])
assays <- list(counts = txi$counts, length = txi$length)
for (b in seq_len(nrep)) {
  assays[[paste0("infRep", b)]] <-
    sapply(infList, function(m) m[, b])
}
se <- SummarizedExperiment(assays = assays,
                           colData = S4Vectors::DataFrame(condition = cond))
rownames(se) <- rownames(txi$counts)
set.seed(1)
se <- scaleInfReps(se, quiet = TRUE)
se <- labelKeep(se, minCount = 3, minN = 3)
se <- se[SummarizedExperiment::rowData(se)$keep, ]
se <- swish(se, x = "condition", quiet = TRUE)
sw <- data.frame(transcript = rownames(se),
                 log2FC = SummarizedExperiment::rowData(se)$log2FC,
                 stat = SummarizedExperiment::rowData(se)$stat,
                 qvalue = SummarizedExperiment::rowData(se)$qvalue)
write.table(sw, paste0(outp, ".swish.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("swish tested:", nrow(sw), "\n")

## ---- edgeR v4 with catchSalmon overdispersion --------------------------
catch <- catchSalmon(paths, verbose = FALSE)
run_edger <- function(counts, tag) {
  y <- DGEList(counts = counts, group = cond)
  keep <- filterByExpr(y, group = cond)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  design <- model.matrix(~cond)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- glmQLFTest(fit, coef = 2)
  tt <- topTags(res, n = Inf)$table
  out <- data.frame(transcript = rownames(tt),
                    log2FC = tt$logFC, PValue = tt$PValue, FDR = tt$FDR)
  write.table(out, paste0(outp, ".", tag, ".tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  cat(tag, "tested:", nrow(out), "\n")
}
scaled <- catch$counts / catch$annotation$Overdispersion
colnames(scaled) <- samples
run_edger(scaled, "edgerOD")
raw <- catch$counts
colnames(raw) <- samples
run_edger(raw, "edgerNaive")
