
# COMPUTE ASSOCIATIONS  ---------------------------------------------------
OmegaCor <- computeAssociations(fitSepTF)
supportLevel <- 0.95

# effect of site
toPlot <- ((OmegaCor[[1]]$support > supportLevel) +
             (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

# function to convert a name
abbreviate_genus <- function(x) {
  sapply(x, function(nm) {
    parts <- unlist(strsplit(nm, "_"))  # split by underscore
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  })
}

# apply to row and column names
rownames(toPlot) <- abbreviate_genus(rownames(toPlot))
colnames(toPlot) <- abbreviate_genus(colnames(toPlot))


# INIT PLOTTING -----------------------------------------------------------
mar=c(1,1,1,1)

width <- 7
height <- 7
if(all_species ==1){
  width <- 20
  height <- 20
}
pdf(file=file.path(input,'results','posterior-omega-coorplot.pdf'),
    width = width,
    height = height)

toPlot[1:3,1:3]
# reorder using eigenvectors
corrplot(toPlot,
         method = "color",
         type='lower',
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 1,
         tl.col = "black",
         order = "hclust",      # <- this reorders nicely
         #addrect=6,
         title = 'Residual associations',
         mar=mar)
dev.off()


# GUILD-LEVEL?  -----------------------------------------------------------
OmegaCor <- computeAssociations(fitSepTF)
supportLevel <- 0.95

# effect of site
toPlot <- ((OmegaCor[[1]]$support > supportLevel) +
             (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

toPlotTraits <- toPlot
for(i in rownames(toPlotTraits)){
  guild <- fitSepTF$TrData$foraging_guild_consensus[which(rownames(fitSepTF$TrData)==i)]
  rownames(toPlotTraits)[which(rownames(toPlotTraits)==i)] <- guild
}
for(i in colnames(toPlotTraits)){
  guild <- fitSepTF$TrData$foraging_guild_consensus[which(rownames(fitSepTF$TrData)==i)]
  colnames(toPlotTraits)[which(colnames(toPlotTraits)==i)] <- guild
}

toPlotTraits[toPlotTraits==0] <- NA

# Get the unique guild names from the row and column names
guilds_row <- unique(rownames(toPlotTraits))
guilds_col <- unique(colnames(toPlotTraits))

# Initialize an empty matrix to store the collapsed means
collapsed <- matrix(NA, nrow = length(guilds_row), ncol = length(guilds_col),
                    dimnames = list(guilds_row, guilds_col))

# Loop through guild pairs and compute the mean correlation
for (r in guilds_row) {
  for (c in guilds_col) {
    # Extract all values for this guild-guild pair
    vals <- toPlotTraits[rownames(toPlotTraits) == r, colnames(toPlotTraits) == c]
    vals_uptri <- vals[upper.tri(vals)]
    collapsed[r, c] <- mean(vals_uptri, na.rm = TRUE)
  }
}


collapsed[is.na(collapsed)] <- 0


# function to convert a name
abbreviate_genus <- function(x) {
  sapply(x, function(nm) {
    parts <- unlist(strsplit(nm, "_"))  # split by underscore
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  })
}

#if(all_species ==1){
  width <- 10
  height <- 10
#}
pdf(file=file.path(input,'results','posterior-omega-guildlevel-coorplot.pdf'),
    width = width,
    height = height)

corrplot(collapsed,
         method = "color",
         type='lower',
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 1,
         tl.col = "black",
         na.label = ' ',
         order = "hclust",
         #addrect = 15, # <- this reorders nicely
         title = 'Residual associations',
         mar=mar)

dev.off()

collapsed
head(collapsed)
# Perform eigen decomposition (PCA is just eigendecomposition here)
eig <- eigen(collapsed)
pca_scores <- eig$vectors %*% diag(sqrt(eig$values))

### PLOT PCA 
pdf(file=file.path(input,'results','posterior-omega-guildlevel-pca.pdf'),
    width = width,
    height = height)
# first 2 PCs for plotting
plot(pca_scores[,1], pca_scores[,2], 
     xlab = "PC1", ylab = "PC2", main = "Guild cooccurrences PCA",
     pch = 19,col='white')
text(pca_scores[,1], pca_scores[,2], labels = colnames(collapsed), cex = 0.9)

# plot PCA without zeros 
keep <- colSums(collapsed) != 0

# Subset the matrix symmetrically
collapsed_pca <- collapsed[keep, keep]
eig <- eigen(collapsed_pca)
pca_scores <- eig$vectors %*% diag(sqrt(eig$values))
plot(pca_scores[,1], pca_scores[,2], 
     xlab = "PC1", ylab = "PC2", main = "Guild cooccurrences PCA",
     pch = 19,col='white')
text(pca_scores[,1], pca_scores[,2], labels = colnames(collapsed_pca), cex = 0.9)

dev.off()


# 
# # vioplot? 
# guild_list <- list()
# for(i in unique(rownames(collapsed))){
#   print(i)
#   species <- rownames(fitSepTF$TrData)[fitSepTF$TrData$foraging_guild_consensus==i]
#   vals <- c(toPlot[which(rownames(toPlot)%in%species),])
#   vals[vals==0] <- NA
#   guild_list[[i]] <- vals
# }
# 
# guild_list$Owls


# migrate-LEVEL?  -----------------------------------------------------------
OmegaCor <- computeAssociations(fitSepTF)
supportLevel <- 0.95

# effect of site
toPlot <- ((OmegaCor[[1]]$support > supportLevel) +
             (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

toPlotTraits <- toPlot
for(i in rownames(toPlotTraits)){
  guild <- fitSepTF$TrData$Migration_a3_DOF[which(rownames(fitSepTF$TrData)==i)]
  rownames(toPlotTraits)[which(rownames(toPlotTraits)==i)] <- guild
}
for(i in colnames(toPlotTraits)){
  guild <- fitSepTF$TrData$Migration_a3_DOF[which(rownames(fitSepTF$TrData)==i)]
  colnames(toPlotTraits)[which(colnames(toPlotTraits)==i)] <- guild
}

toPlotTraits[toPlotTraits==0] <- NA
toPlotTraits

# Get the unique guild names from the row and column names
guilds_row <- unique(rownames(toPlotTraits))
guilds_col <- unique(colnames(toPlotTraits))

# Initialize an empty matrix to store the collapsed means
collapsed <- matrix(NA, nrow = length(guilds_row), ncol = length(guilds_col),
                    dimnames = list(guilds_row, guilds_col))

# Loop through guild pairs and compute the mean correlation
for (r in guilds_row) {
  for (c in guilds_col) {
    # Extract all values for this guild-guild pair
    vals <- toPlotTraits[rownames(toPlotTraits) == r, colnames(toPlotTraits) == c]
    vals_uptri <- vals[upper.tri(vals)]
    print(paste0(r,c))
    print(vals_uptri)
    collapsed[r, c] <- mean(vals_uptri, na.rm = TRUE)
  }
}

collapsed[is.na(collapsed)] <- 0


# function to convert a name
abbreviate_genus <- function(x) {
  sapply(x, function(nm) {
    parts <- unlist(strsplit(nm, "_"))  # split by underscore
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  })
}

corrplot(collapsed,
         method = "color",
         type='lower',
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 1,
         tl.col = "black",
         order = "hclust",
         #addrect = 15, # <- this reorders nicely
         title = 'Residual associations',
         mar=mar)


collapsed
# Perform eigen decomposition (PCA is just eigendecomposition here)
eig <- eigen(collapsed)
pca_scores <- eig$vectors %*% diag(sqrt(eig$values))

# Optional: first 2 PCs for plotting
plot(pca_scores[,1], pca_scores[,2], 
     xlab = "PC1", ylab = "PC2", main = "Species PCA")
text(pca_scores[,1], pca_scores[,2], labels = colnames(collapsed), cex = 0.7, pos=3)

# look what happens with species
toPlot
# Perform eigen decomposition (PCA is just eigendecomposition here)
eig <- eigen(toPlot)
pca_scores <- eig$vectors %*% diag(sqrt(eig$values))
is.na(pca_scores) <- 0

# Optional: first 2 PCs for plotting
plot(pca_scores[,1], pca_scores[,2], 
     xlab = "PC1", ylab = "PC2", main = "Species PCA")

text(pca_scores[,1], pca_scores[,2], labels = colnames(toPlot), cex = 0.7, pos=3)
