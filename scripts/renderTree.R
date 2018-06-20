library(ggtree)
library(phangorn)
library(ggplot2)

# input files
tree_file <- snakemake@input[[1]]
out_file <- snakemake@output[[1]]

# read tree file
tree <- read.tree(tree_file)

# midpoint root tree if there is more than 1 node
if (tree$Nnode > 1){
  tree <- midpoint(tree)
}

# draw tree
p <- ggtree(tree, branch.length = 1) +
  geom_tiplab(offset = .05) +
  geom_tippoint(size=3) +
  xlim(0, 4) +
  geom_treescale(width = 0.1, offset=0.25)

# save output file
ggsave(out_file, plot = p, height = 1*length(tree$tip.label))
