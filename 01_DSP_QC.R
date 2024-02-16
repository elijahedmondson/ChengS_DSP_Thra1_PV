library(readxl)
library(NanoStringNCTools)
library(GeomxTools)
library(knitr)

output_prefix<-"Cheng_WTA1"
projectname<-"Cheng_WTA1"

datadir<-"F:/GeoMX KPC/Cheng_WTA1/raw_data"
DCCdir<-"DCC-20230817"
PKCfilename<-"Mm_R_NGS_WTA_v1.0.pkc"
WorkSheet<-"final2_full.xlsx"
final <- read_excel("F:/GeoMX KPC/Cheng_WTA1/raw_data/final2_full.xlsx")
output_dir<-"processed_data/"

DCCFiles <- list.files(file.path(datadir , DCCdir), pattern=".dcc$", full.names=TRUE)
PKCFiles <- file.path(datadir, PKCfilename)
SampleAnnotationFile <- file.path(datadir, WorkSheet)

myData<-readNanoStringGeoMxSet(dccFiles = DCCFiles,
                               pkcFiles = PKCFiles,
                               phenoDataFile = SampleAnnotationFile,
                               phenoDataSheet = "Template",
                               phenoDataDccColName = "Sample_ID",
                               protocolDataColNames = c("aoi", "roi"),
                               experimentDataColNames = c("panel"))

#Shift counts to one to mimic how DSPDA handles zero counts
myData <- shiftCountsOne(myData, elt="exprs", useDALogic=TRUE)
pkcs <- annotation(myData)
modules <- gsub(".pkc", "", pkcs)

########
####  Segment QC
#
#
# Before excluding any low-performing ROI/AOI segments, we visualize the distributions of the data for the different QC parameters.  The cutoffs used are:
#
#   * **Raw sequencing reads**: segments with <1000 raw reads are removed.
# * **% Aligned,% Trimmed, or % Stitched sequencing reads**: segments below ~80% for one or more of these QC parameters are removed.
# * **% Sequencing saturation ([1-deduplicated reads/aligned reads]%)**: segments below ~50% require additional sequencing to capture full sample diversity and are not typically analyzed until improved.
# * **Negative Count**: this is the geometric mean of the several unique negative probes in the GeoMx panel that do not target mRNA and establish the background count level per segment; segments with low negative counts (1-10) are not necessarily removed but may be studied closer for low endogenous gene signal and/or insufficient tissue sampling.
# * **No Template Control (NTC) count**: values >1,000 could indicate contamination for the segments associated with this NTC; however, in cases where the NTC count is between 1,000- 10,000, the segments may be used if the NTC data is uniformly low (e.g. 0-2 counts for all probes).
# * **Nuclei**: >100 nuclei per segment is generally recommended; however, this cutoff is highly study/tissue dependent and may need to be reduced; what is most important is consistency in the nuclei distribution for segments within the study.  We use 20 in this case.
# * **Area**: generally correlates with nuclei; a strict cutoff is not generally applied based on area.
#
# _***Note***:  You may need to change these cutoffs depending on your experimental setup._
########

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned to known targets (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of cells observed in a segment (100)
       minArea = 1000)         # Minimum segment area (5000)
myData <-
  setSegmentQCFlags(myData, qcCutoffs = QC_params)

# Collate QC Results
QCResults <- protocolData(myData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))








col_by <- "COMP2"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 200) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 7) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}


## ----plotQCHist, warning = FALSE, message = FALSE-----------------------------
QC_histogram(sData(myData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(myData), "Stitched (%)", col_by, 80)
QC_histogram(sData(myData), "Aligned (%)", col_by, 75)
QC_histogram(sData(myData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(myData), "area", col_by, 10, scale_trans = "log10")
QC_histogram(sData(myData), "Nuclei count", col_by, 10)

#QC_histogram(sData(myData), "Aligned", col_by, 10000)
# calculate the negative geometric means for each module
negativeGeoMeans <-
  esBy(negativeControlSubset(myData),
       GROUP = "Module",
       FUN = function(x) {
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs")
       })
protocolData(myData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(myData)[, negCols] <- sData(myData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(myData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}


# detatch neg_geomean columns ahead of aggregateCounts call
pData(myData) <- pData(myData)[, !colnames(pData(myData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
print("No Template Control (NTC) wells are essential for detecting contamination or non-specific amplification")

kable(table(NTC_Count = sData(myData)$NTC),
      col.names = c("NTC Count", "# of Segments"))



## ----QCSummaryTable, results = "aexprs()## ----QCSummaryTable, results = "asis"-----------------------------------------
kable(QC_Summary, caption = "QC Summary Table for each Segment")

## ----removeQCSampleProbe, eval = TRUE-----------------------------------------
myData <- myData[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(myData)

## ----setbioprobeqcflag,  eval = TRUE------------------------------------------
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to
# FALSE if you do not want to remove local outliers
myData <- setBioProbeQCFlags(myData,
                             qcCutoffs = list(minProbeRatio = 0.1,
                                              percentFailGrubbs = 20),
                             removeLocalOutliers = TRUE)

ProbeQCResults <- fData(myData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

## ----bioprobeQCTable, echo = FALSE, results = "asis"--------------------------
kable(qc_df, caption = "Probes flagged or passed as outliers")


## ----excludeOutlierProbes-----------------------------------------------------
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <-
  subset(myData,
         fData(myData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(myData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
myData <- ProbeQCPassed

## ----aggregateCounts, eval = TRUE---------------------------------------------
# Check how many unique targets the object has
length(unique(featureData(myData)[["TargetName"]]))

# collapse to targets
target_myData <- aggregateCounts(myData)
dim(target_myData)
exprs(target_myData)[100:103, 1:4]


## ----calculateLOQ, eval = TRUE------------------------------------------------
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_myData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_myData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_myData)[, vars[1]] *
             pData(target_myData)[, vars[2]] ^ cutoff)
  }
}
pData(target_myData)$LOQ <- LOQ

head(pData(target_myData)$LOQ)

## ----LOQMat, eval = TRUE------------------------------------------------------
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_myData)$Module == module
  Mat_i <- t(esApply(target_myData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_myData)$TargetName, ]

## ----segDetectionBarplot------------------------------------------------------
# Save detection rate information to pheno data
pData(target_myData)$GenesDetected <-
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_myData)$GeneDetectionRate <-
  pData(target_myData)$GenesDetected / nrow(target_myData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_myData)$DetectionThreshold <-
  cut(pData(target_myData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15,1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_myData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = genotype)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")


library(DT)
tmp_data<-pData(target_myData)
tmp_data<-tmp_data[,c("class","region","GenesDetected","GeneDetectionRate")]
datatable(tmp_data,
          'Table showing the genes detected for each ROI',
          options=list(pageLength = 10)
) %>%  formatRound(c(4), 2)

# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_myData)$DetectionThreshold,
            pData(target_myData)$class))

## ----filterSegments-----------------------------------------------------------
target_myData <-
  target_myData[, pData(target_myData)$GeneDetectionRate >= .05]                    ########EFE excludes samples with low gene detection
pData(target_myData)[,24:27]

dim(target_myData)

target_myData@phenoData@data$genotype









library(dplyr)
library(ggforce)



## ----goi detection------------------------------------------------------------


library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_myData)]
fData(target_myData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_myData)$DetectionRate <-
  fData(target_myData)$DetectedSegments / nrow(pData(target_myData))

# Gene of interest detection table

goi <- c("Kras", "Trp53", "Cd274", "Cd8a", "Cd68", "Epcam","Cre",
         "Tgfb1","Ccnd1","Cdk4","Thra","Thrb","Il33",	"Ctnnb1",
         "Cdk6","Myc","Smad3","Smad4","Cdkn1a","Agr2","Foxa2")


goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_myData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_myData)[goi, "DetectionRate"]))

## ----tableGOI, echo = FALSE, results = "asis"---------------------------------
kable(goi_df, caption = "Detection rate for Genes of Interest", align = "c",
      col.names = c("Gene", "Detection, # Segments", "Detection Rate, % of Segments"))

## ----plotDetectionRate, eval = TRUE-------------------------------------------
#Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 3, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_myData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_myData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

## ----subsetGenes, eval = TRUE-------------------------------------------------
# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_myData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_myData <-
  target_myData[fData(target_myData)$DetectionRate >= 0.035 |                    ########EFE change to include additional genes?
                  fData(target_myData)$TargetName %in% neg_probes, ]
dim(target_myData)

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_myData)]

## ----previewNF, fig.width = 8, fig.height = 8, fig.wide = TRUE, eval = TRUE, warning = FALSE, message = FALSE----
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "class"
Stat_data <-
  data.frame(row.names = colnames(exprs(target_myData)),
             Segment = colnames(exprs(target_myData)),
             Annotation = pData(target_myData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_myData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_myData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) +
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

## ----normalizeObject, eval = TRUE---------------------------------------------
# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_myData <- normalize(target_myData , data_type = "RNA",
                           norm_method = "quant",
                           desiredQuantile = .75,
                           toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_myData <- normalize(target_myData , data_type = "RNA",
                           norm_method = "neg",
                           fromElt = "exprs",
                           toElt = "neg_norm")

## ----normplot, fig.small = TRUE-----------------------------------------------
#visualize the first 10 segments with each normalization method
# boxplot(exprs(target_myData)[,1:79],
#         col = "#9EDAE5", main = "Raw Counts",
#         log = "y", names = 1:79, xlab = "Segment",
#         ylab = "Counts, Raw")
#
# boxplot(assayDataElement(target_myData[,1:79], elt = "q_norm"),
#         col = "#2CA02C", main = "Q3 Norm Counts",
#         log = "y", names = 1:79, xlab = "Segment",
#         ylab = "Counts, Q3 Normalized")
#
# boxplot(assayDataElement(target_myData[,1:79], elt = "neg_norm"),
#         col = "#FF7F0E", main = "Neg Norm Counts",
#         log = "y", names = 1:79, xlab = "Segment",
#         ylab = "Counts, Neg. Normalized")

## ----dimReduction, eval = TRUE------------------------------------------------
library(umap)
library(Rtsne)

shapes = c(15,16,17,18,19,20,21,22,23,24,25)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_myData , elt = "q_norm"))),
       config = custom_umap)
pData(target_myData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

umapplot <-ggplot(pData(target_myData),
                  aes(x = UMAP1, y = UMAP2, label=animal, shape = genotype, color = region_geno)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  theme_bw() #+
  #theme(text = element_text(size = 10)) +
  #theme(legend.position="none")
umapplot
ggsave(umapplot, file="C:/Users/edmondsonef/Desktop/umap.png", width = 18, height = 10, units = "in", bg = "white")


# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_myData , elt = "q_norm"))),
        perplexity = ncol(target_myData)*.15)
pData(target_myData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
tsneplot <- ggplot(pData(target_myData),
       aes(x = tSNE1, y = tSNE2, color = age_young_old, shape = genotype, label=`region`, size = 5)) +
  geom_point(size = 3) +geom_text(hjust=1.1, vjust=0.2)+
  theme_bw()
  #theme(legend.position="none")
tsneplot

ggsave(tsneplot, file="C:/Users/edmondsonef/Desktop/tsne.png", width = 18, height = 10, units = "in", bg = "white")


## run PCA
PCAx<-1
PCAy<-2
PCAxy <- c(as.integer( PCAx ),as.integer( PCAy) ) # selected principal components


pca.object <- prcomp(t(log2(assayDataElement(target_myData , elt = "q_norm"))))
pcaData = as.data.frame(pca.object$x[, PCAxy]);
pData(target_myData)[, c("PC1", "PC2")] <- pcaData[,c(1,2)]
percentVar=round(100*summary(pca.object)$importance[2, PCAxy],0)


pcaplot <- ggplot(pData(target_myData),
       aes(x = PC1, y = PC2, color=region, shape = age_young_old, label=`genotype`)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  xlab(paste0("PC", PCAx ,": ", percentVar[1], "% variance")) +
  ylab(paste0("PC", PCAy ,": ", percentVar[2], "% variance")) +
  theme_bw()
  #theme(legend.position="none")
pcaplot

ggsave(pcaplot, file="C:/Users/edmondsonef/Desktop/pcaKPC2.png", width = 18, height = 10, units = "in", bg = "white")




## ----CVheatmap, eval = TRUE, echo = TRUE, fig.width = 8, fig.height = 6.5, fig.wide = TRUE----
library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_myData,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:50]

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.9949)]


pheatmap(assayDataElement(target_myData[GOI, ], elt = "log_q"),
         scale = "row",
         main="Heatmap of high CV (>0.99) genes",
         show_rownames = T, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average", #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid"
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_myData)[, c("age_young_old", "region", "genotype")],
         cutree_rows = 4,
         cutree_cols = 6)


save(final, target_myData, neg_probes, file = "F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_2_16_2024.RData")

