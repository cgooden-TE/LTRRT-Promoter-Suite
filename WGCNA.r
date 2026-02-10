## WGCNA
library(WGCNA)

adj_logCPM <- readRDS("/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/WGCNA/adj_logCPM.rds")
meta <- readRDS("/home/caleb/data/PaperWritingReruns/Dec2025_BugFix/WGCNA/pca_metadata.rds")
#### 1
options(stringsAsFactors = FALSE)
# Enable threading if available (Linux/macOS)
enableWGCNAThreads() # on some systems allowWGCNAThreads() is used instead
# Make sure columns (samples) order matches meta
stopifnot(all(colnames(adj_logCPM) %in% rownames(meta)))
meta <- meta[colnames(adj_logCPM), , drop = FALSE]
stopifnot(identical(colnames(adj_logCPM), rownames(meta)))

# WGCNA expects samples in rows, genes in columns
datExpr0 <- t(adj_logCPM)

# Basic sanity check and removal of bad samples/genes
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]

# Remove genes with zero variance and keep the more variable half (tune as needed)
madVals <- apply(datExpr0, 2, mad, na.rm = TRUE)
keep <- madVals >= quantile(madVals, 0.50, na.rm = TRUE) # top 50% by MAD
datExpr <- datExpr0[, keep]

# Optional: remove extreme outlier samples
sampleTree <- hclust(dist(datExpr), method = "average")
# Plot and decide a static cut height if needed:
# plot(sampleTree); abline(h = someHeight, col=2)
# clust <- cutreeStatic(sampleTree, cutHeight = someHeight, minSize = 10)
# datExpr <- datExpr[clust==1,]

#### 2
traits <- data.frame(
    platform = factor(meta[colnames(adj_logCPM), "platform"]),
    tissue   = factor(meta[colnames(adj_logCPM), "tissue_umap"])
)
rownames(traits) <- rownames(datExpr)

# One-hot encode categorical traits for correlations
traitMat <- model.matrix(~ 0 + platform + tissue, data = traits)
colnames(traitMat) <- sub("^platform", "plat_", sub("^tissue", "tis_", colnames(traitMat)))

#### 3
powers <- c(1:20, seq(22, 40, by = 2))
sft <- pickSoftThreshold(datExpr,
    powerVector = powers,
    corFnc = "bicor",
    corOptions = list(maxPOutliers = 0.1),
    networkType = "signed", verbose = 5
)

fi <- sft$fitIndices

pick_power <- function(fi, r2_min = 0.5, meanK_max = 1000, slope_max = 0) {
    ok <- fi$SFT.R.sq >= r2_min & fi$mean.k. <= meanK_max & fi$slope <= slope_max
    if (any(ok)) {
        return(min(fi$Power[ok]))
    }
    # fallback: smallest power with negative slope and minimal mean connectivity
    ok2 <- fi$slope <= slope_max
    if (any(ok2)) {
        return(fi$Power[ok2][which.min(fi$mean.k.[ok2])])
    }
    fi$Power[which.max(fi$SFT.R.sq)]
}

softPower <- pick_power(fi)
softPower

#### 4
bwnet <- blockwiseModules(
    datExpr,
    power = softPower,
    corType = "bicor", # robust correlation
    maxPOutliers = 0.1,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30, # tune if needed
    mergeCutHeight = 0.25, # merge similar modules
    numericLabels = FALSE, # get color labels
    pamRespectsDendro = TRUE,
    saveTOMs = FALSE,
    verbose = 3
)

moduleColors <- bwnet$colors
mes <- orderMEs(bwnet$MEs)

# Module–trait correlations
mevsTraitCor <- bicor(mes, traitMat, maxPOutliers = 0.1, use = "pairwise.complete.obs")
mevsTraitP <- corPvalueStudent(mevsTraitCor, nrow(datExpr))

# Gene significance (GS) and module membership (kME), for, e.g., a Leaf contrast:
# pick a trait column of interest
leafCol <- grep("^tis_Leaf$", colnames(traitMat), value = TRUE)
gs.leaf <- bicor(datExpr, traitMat[, leafCol], maxPOutliers = 0.1)
kME <- bicor(datExpr, mes, maxPOutliers = 0.1)

# Hub genes per module by eigengene connectivity:
hubs <- chooseTopHubInEachModule(datExpr, moduleColors,
    power = softPower,
    type = "signed", corFnc = "bicor"
)
head(hubs)

saveRDS(list(
    bwnet = bwnet,
    moduleColors = moduleColors,
    mes = mes,
    traits = traits,
    datExpr = datExpr,
    softpwr = softPower
), file = "WGCNA_fullResults.rds")

