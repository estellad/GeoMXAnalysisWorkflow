
# Taken from standR vignette ----------------------------------------------
library(standR)
library(SpatialExperiment)
library(limma)
library(ExperimentHub)
library(ggalluvial)

eh <- ExperimentHub()
query(eh, "standR")

countFile <- eh[["EH7364"]]
sampleAnnoFile <- eh[["EH7365"]]
featureAnnoFile <- eh[["EH7366"]]

spe_standR <- readGeoMx(countFile, sampleAnnoFile, featureAnnoFile = featureAnnoFile, hasNegProbe = TRUE)
# this is different from dkd_spe_subset

# Check the overlap of these
overlap_genes <- intersect(rownames(spe_standR), rownames(dkd_spe_subset)) # 2978

overlap_ROIs <- intersect(colnames(spe_standR), colnames(dkd_spe_subset)) # 70
# so dkd_spe_subset is a subset from spe_standR from ExperimentHub

spe_standR_subset <- spe_standR[rownames(spe_standR) %in% overlap_genes, colnames(spe_standR) %in% overlap_ROIs]

all_equal(counts(spe_standR_subset), counts(dkd_spe_subset))



# Metadata visualization
old_region <- colData(spe_standR)$region
new_region <- paste0(colData(spe_standR)$region,"_",colData(spe_standR)$SegmentLabel) |>
  (\(.) gsub("_Geometric Segment","",.))() |>
  paste0("_",colData(spe_standR)$pathology) |>
  (\(.) gsub("_NA","_ns",.))()
plotSampleInfo(spe_standR, column2plot = c("SlideName","disease_status","regions"))


# Gene level QC
spe_standR <- addPerROIQC(spe_standR, rm_genes = TRUE)
plotGeneQC(spe_standR, ordannots = "regions", col = regions, point_size = 2)

# ROI level QC
plotROIQC(spe_standR, y_threshold = 50000, col = SlideName)
spe_standR <- spe_standR[,rownames(colData(spe_standR))[colData(spe_standR)$lib_size > 50000]]

# RLE
plotRLExpr(spe_standR_standR, ordannots = "SlideName", assay = 2, col = SlideName)




