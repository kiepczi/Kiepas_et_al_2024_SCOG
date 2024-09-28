library("ape")
library('ggplot2')
library('ggtree')
library("phytools")


tree <- ape::read.tree('../output/tree/04_tbe.raxml.support')
tree <-midpoint.root(tree)
options(repr.plot.width=8, repr.plot.height=6)

p <- ggtree(tree, layout='c')

Accession <- get_taxa_name(p)

# Create a data frame with a column named "accession"
data <- data.frame(accession = Accession)

# Specify the file path for saving the CSV file
file_path <- "../output/SCOG_tree_node_order.csv"

# Save the data frame to a CSV file without row names
write.csv(data, file = file_path, row.names = FALSE)
