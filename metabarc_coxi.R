## cargar librerias y datos
library(dada2)
set.seed(100)
getwd()
fnAs <- sort(list.files(pattern = "_A.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnAs), ".f"), `[`, 1)
sample.names
filtAs <- file.path("filtered", paste0(sample.names, "_A_filt.fastq.gz"))
# plotQualityProfile(fnAs[17:32])
isMulti <- TRUE

# filtrar y trimar
out <- filterAndTrim(fnAs, filtAs,
  maxEE = 2, rm.phix = TRUE,
  multithread = isMulti,
  truncQ = 20,
  compress = TRUE
)
# tasas de error
errA <- learnErrors(filtAs, verbose = TRUE, multithread = isMulti)
# dereplicacion
derepAs <- derepFastq(filtAs, verbose = TRUE)
# inferencia
dadaAs <- dada(derepAs, err = errA, multithread = isMulti)
seqtab <- makeSequenceTable(dadaAs)
# remove quimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = isMulti, verbose = TRUE)
saveRDS(seqtab.nochim, "meta_all.rds")
save.image("meta_all.RData")
# summary table
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaAs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "seq_nonchim")
rownames(track) <- sample.names
write.table(x = track, file = "reads-stats.txt", sep = "\t", quote = FALSE)

# TAXONOMIa

tax <- assignTaxonomy(seqtab.nochim, "database_meta/Metacoxi_tax_new_shey.fasta.gz", multithread = isMulti)
write.table(tax, file = "tax.tsv", sep = "\t")
species <- addSpecies(tax, "database_meta/Metacoxi_spe_new_shey.fasta.gz")
write.table(species, file = "species.tsv", sep = "\t")

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

asv_tax <- tax
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.nov.tsv", sep = "\t", quote = F, col.names = NA)

asv_tax_sp <- species
row.names(asv_tax_sp) <- sub(">", "", asv_headers)
write.table(asv_tax_sp, "ASVs_taxonomy_sp.nov.tsv", sep = "\t", quote = F, col.names = NA)

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.nov.fa")

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.nov.tsv", sep = "\t", quote = F, col.names = NA)
