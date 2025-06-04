rm(list=ls())

# GETTING STARTED ---------------------------------------------------------
library(ape)
library(supportScripts)


# LOADING DATA ------------------------------------------------------------
# read species # test without subspecies
species <- read.csv('data/0_traits/birdTraits-merged-subset.csv')

# select columns
species_selected <- species %>%
  select(order_IOC,family_IOC,genus_IOC,species_IOC)

# change the names to underscores  
species_selected$class_IOC <- 'Aves'

# change to factors 
species_selected <- species_selected %>%
  mutate(across(where(is.character), as.factor)) %>% 
  distinct()

# BUILD PDMATRIX ----------------------------------------------------------
# Function to compute pairwise distance
taxonomic_distance <- function(row1, row2, ranks) {
  for (i in seq_along(ranks)) {
    rank <- ranks[i]
    if (row1[[rank]] != row2[[rank]]) {
      return(2 * (length(ranks) - i + 1))  # two paths from each species to MRCA
    }
  }
  return(0)  # same species
}

ranks <- c('class_IOC','order_IOC','family_IOC','genus_IOC','species_IOC')
# Create pairwise distance matrix
species_names <- species_selected$species_IOC
n <- length(species_names)
dist_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(species_names, species_names))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      dist_mat[i, j] <- taxonomic_distance(species_selected[i, ], species_selected[j, ], ranks)
    }
  }
}

# View or save
print(dist_mat)
dist_mat['Motacilla flava','Motacilla alba'] # should be 2? 
dist_mat['Galerida cristata','Lullula arborea'] # should be 4? 
dist_mat['Ciconia ciconia','Phalacrocorax carbo'] # should be 8? 


# Convert matrix to a distance object
dist_obj <- as.dist(dist_mat)

# Use hierarchical clustering (e.g., UPGMA)
hc <- hclust(dist_obj, method = "average")

# Convert to a phylogenetic tree
library(ape)
tree <- as.phylo(hc)

# Plot it
# Create a vector of unique colors
family_colors <- setNames(rainbow(length(unique(species_selected$order_IOC))), unique(species_selected$order_IOC))

# Map colors to each tip based on family
tip_colors <- family_colors[species_selected$order_IOC]

pdf('data/1_preprocessing/Taxonomy/tree_fromPD.pdf',
    width = 8,
    height = 8)
plot(tree,tip.color = tip_colors, cex = 0.6, type = 'fan')
dev.off()

write.tree(tree,'data/1_preprocessing/Taxonomy/tree_fromPD.tre')

# BUILD TREE -----------------------------------------------------
metaGenerator('data/1_preprocessing/Taxonomy/')
