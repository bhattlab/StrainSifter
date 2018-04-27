library(ggtree)
library(phangorn)
library(ggplot2)
library(svglite)

# input
tree.file <- snakemake@input[[1]]
meta.file <- snakemake@input[[2]]
out.file <- snakemake@output[[1]]

# input files
# tree.file <- '/Users/Fiona/Desktop/tes.tree'
# meta.file <- '/Users/Fiona/Google Drive/Bhatt lab/Projects/Bacteremia project/bacteremia_metadata_all.csv'

# read tree file
tree <- read.tree(tree.file)

# midpoint root tree if more than 1 node
if (tree$Nnode > 1){
  tree <- midpoint(tree)
}

# read metadata
meta.table <- read.table(meta.file, header = TRUE, sep = ',', comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)

meta.table$`Organism (taxonomy)` <- gsub(' ', '_', meta.table$`Organism (taxonomy)`)

stools <- data.frame(
                Label = meta.table$`Stool label for phylo tree`,
                Sample = paste("S", meta.table$`Stool number`, sep = ''),
                PatientCode = meta.table$`Patient code`
                )
isolates <- unique(data.frame(
                Label = meta.table$`Isolate label for phylo tree`,
                Sample = meta.table$`Organism code`,
                PatientCode = meta.table$`Patient code`
                ))

tree.meta <- rbind(stools, isolates)
tree.meta$PatientCode <- as.character(tree.meta$PatientCode)

# draw tree and scale bar
# ggtree(tree) + geom_treescale(width = 0.1) + geom_tippoint() + geom_tiplab()

# edit tip labels
tree$tip.label <- sub('\\..+', '', tree$tip.label)
tree$tip.label <- as.character(tree.meta$Label[match(tree$tip.label, tree.meta$Sample)])

# draw tree
p <- ggtree(tree) %<+% tree.meta + geom_tiplab(offset = .005) + geom_tippoint(aes(color=PatientCode), size=3) +
  xlim(NA, 1) + geom_treescale(width = 0.1)

# ggsave(out.file, plot = p, height = .5*length(tree$tip.label), useDingbats=FALSE)
ggsave(out.file, plot = p, height = .5*length(tree$tip.label))
