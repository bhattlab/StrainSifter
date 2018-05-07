library(ggtree)
library(phangorn)
library(ggplot2)
library(svglite)

# input files
tree.file <- snakemake@input[[1]]
meta.file <- snakemake@input[[2]]
out.file <- snakemake@output[[1]]

# read tree file
tree <- read.tree(tree.file)

# midpoint root tree if there is more than 1 node
if (tree$Nnode > 1){
  tree <- midpoint(tree)
}

# draw tree
p <- ggtree(tree) +
  geom_tiplab(offset = .005) +
  geom_tippoint(size=3) +
  xlim(NA, 1) +
  geom_treescale(width = 0.1)

# save output file
ggsave(out.file, plot = p, height = .5*length(tree$tip.label))
