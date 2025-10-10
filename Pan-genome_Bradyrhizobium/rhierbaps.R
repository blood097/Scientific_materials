library(rhierbaps)
library(phytools)
library(ggtree)

snp.matrix <- load_fasta("Core Gene Alignment.fasta")
hb.results <- hierBAPS(snp.matrix, max.depth = 2,  n.pops = 20, quiet = TRUE)

write.csv(hb.results$partition.df, file = "hierbaps_partition.csv", col.names = TRUE,
          row.names = FALSE)
saveRDS(hb.results, file = "hierbaps.RDS")

tree <- phytools::read.newick("IQTREE.nwk")
hb.results <- read.csv("hierbaps_partition.csv")
colnames(hb.results) <- sub("\\.", " ", colnames(hb.results))

gg <- ggtree(tree,layout="roundrect")
gg <- gg%<+%hb.results#$partition.df
gg <- gg+geom_tippoint(aes(color=factor(gg$data$`level 1`)))
gg + geom_tiplab()

#plot_sub_cluster(hb.results, parsnp, level = 2, sub.cluster = 1)
