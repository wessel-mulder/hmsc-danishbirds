rm(list=ls())

# GETTING STARTED ---------------------------------------------------------
library(ape)
library(supportScripts)


# LOADING DATA ------------------------------------------------------------
# read species # test without subspecies
species <- read.csv('data/1_preprocessing/Taxonomy/taxonomy.csv')

# select columns
species_selected <- species %>%
  select(Species2_AVONET,genus,Family2_AVONET,Order2_AVONET)

# change the names to underscores  
species_selected$class_IOC <- 'Aves'

# change to factors 
species_selected <- species_selected %>%
  mutate(across(where(is.character), as.factor)) %>% 
  distinct()

# BUILD PDMATRIX ----------------------------------------------------------
# Function to compute pairwise distance
tax2dist <- function(lookup,
                     tax_distance,
                     precompute_dist = TRUE) {
  tax_distance <- sort(tax_distance)
  tax_cols <- names(tax_distance)[-length(tax_distance)]
  if (length(intersect(tax_cols, colnames(lookup))) != length(tax_cols))
    stop("The columns in the taxonomy have to include the columns mentioned in the distance vector")
  lookup <- lookup[, tax_cols]
  
  
  if (is.numeric(precompute_dist)) {
    n <- apply(lookup, 2, function(x) length(unique(x)))
    S <- max(n)
    precompute_dist <- ifelse(S > precompute_dist, FALSE, TRUE)
  }
  
  # Calculate distance matrix
  
  if (precompute_dist) {
    entries <- row.names(lookup)
    n <- length(entries)
    dist <- matrix(NA, nrow = n, ncol = n)
    colnames(dist) <- unlist(lookup[, 1])
    row.names(dist) <- unlist(lookup[, 1])
    other <- tax_distance[length(tax_distance)]
    
    for (i in seq_along(entries)) {
      for (j in seq_along(entries)) {
        row <- as.character(lookup[i, ])
        column <- as.character(lookup[j, ])
        if (any(row == column))
          dist[i, j] <- tax_distance[min(which(row == column))] else
            dist[i, j] <- other
      }
    }
    
    return(dist)
    
    # Don't calculate distance matrix
    
  }
}

hierarchy <- c('species','genus','family','order','class')
colnames(species_selected) <- hierarchy
tax_distance <- c(species=0,genus=1,family=2,order=3,class=4)


dist_mat <-tax2dist(species_selected,tax_distance)



# View or save
print(dist_mat)
dist_mat['Motacilla flava','Motacilla alba'] # should be 1? 
dist_mat['Galerida cristata','Lullula arborea'] # should be 2? 
dist_mat['Nucifraga caryocatactes','Oenanthe oenanthe'] # should be 3? 
dist_mat['Ciconia ciconia','Phalacrocorax carbo'] # should be 4? 


# Convert matrix to a distance object
dist_obj <- as.dist(dist_mat)

# Use hierarchical clustering (e.g., UPGMA)
hc <- hclust(dist_obj, method = "average")

# Convert to a phylogenetic tree
library(ape)
tree <- as.phylo(hc)

# Plot it
# Create a vector of unique colors
family_colors <- setNames(rainbow(length(unique(species_selected$order))), unique(species_selected$order))

# Map colors to each tip based on family
tip_colors <- family_colors[species_selected$order]

pdf('data/1_preprocessing/Taxonomy/tree_fromPD.pdf',
    width = 8,
    height = 8)
plot(tree,tip.color = tip_colors, cex = 0.6, type = 'fan')
dev.off()

write.tree(tree,'data/1_preprocessing/Taxonomy/tree_fromPD.tre')

# BUILD TREE -----------------------------------------------------
metaGenerator('data/1_preprocessing/Taxonomy/')
