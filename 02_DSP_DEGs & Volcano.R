library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)
library(enrichplot)
library(data.table)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationHub)
library(GOSemSim)
library(clusterProfiler)
library(GOSemSim)
library(ggwordcloud)
library(ggplot2)
library(cowplot)
library(ReactomePA)
library(DOSE)
library(msigdbr)
library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)
library(topGO)
library(scales) # for percent
library(reshape2)
library(cowplot)
library(umap)
library(Rtsne)

load("F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_10_25_2023.RData")

pData(target_myData)$region_geno

# convert test variables to factors
pData(target_myData)$testRegion <- factor(pData(target_myData)$genotype)#, c("glands"))
pData(target_myData)[["slide"]] <-  factor(pData(target_myData)[["MHL number"]])
assayDataElement(object = target_myData, elt = "log_q") <- assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
results <- c()
for(region in c("glands", "mucosa", "stroma")) {
  ind <- pData(target_myData)$region == region
  mixedOutmc <-
    mixedModelDE(target_myData[, ind], elt = "log_q",
                 modelFormula = ~ testRegion + (1 + testRegion | slide), ##co-exist in a given tissue section
                 #modelFormula = ~ testRegion + (1 | slide),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  r_test$Gene <-  unlist(lapply(colnames(mixedOutmc), rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}


results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color, levels = c("NS or FC < 0.5", "P < 0.05", "FDR < 0.05", "FDR < 0.001"))
dplyr::count(results, FDR < 0.001)
dplyr::count(results, `Pr(>|t|)` < 0.05)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)

write <- dplyr::filter(results, FDR < 0.05 & abs(results$Estimate) > 0.5)
write.csv(write, "F:/GeoMX KPC/Cheng_WTA1/processed_data/DEG.csv")

top_g <- c()
for(genotype in c("PV/N", "N/N")) {
  ind <- results$Subset == genotype
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:50]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:50]])
}
top_g <- unique(top_g)
top_g

#reverse log fold change to fit with label
results$Estimate1 <- results$Estimate*(-1)
# Graph results
volc_plot <- ggplot(results,                                                             ###CHANGE
       aes(x = Estimate1, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Glands <- log2(FC) -> Stroma",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.001 & abs(results$Estimate1) > 0.7),
                  size = 6, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, #lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom")+
  facet_wrap(~Subset, scales = "free_y")

volc_plot




library(VennDiagram)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)
library(tidyverse)
library(plyr)
library(GGally)
library(ggplot2)
library(tidyverse)
library(gapminder)
library(dplyr)
library(RColorBrewer)

myCol <- brewer.pal(4, "Pastel2")


write <- dplyr::filter(results, FDR < 0.001 & abs(results$Estimate) > 0.5)
head(write)
PVN_stroma <- dplyr::filter(write, Subset == "PV/N" & Estimate > 0.5)
PVN_glands <- dplyr::filter(write, Subset == "PV/N" & Estimate < -0.5)
NN_stroma <- dplyr::filter(write, Subset == "N/N" & Estimate > 0.5)
NN_glands <- dplyr::filter(write, Subset == "N/N" & Estimate < -0.5)

gene_list <- list(PVN_glands = PVN_glands$Gene,
                  NN_glands = NN_glands$Gene,
                  PVN_stroma = PVN_stroma$Gene,
                  NN_stroma = NN_stroma$Gene)
VennDiagram <- venn.diagram(x = gene_list,
                            fill = myCol,
                            #cat.col = c("red", "green"),
                            cex = 2,lty = "blank",
                            cat.cex = 2,
                            filename = NULL)
cowplot::plot_grid(VennDiagram)








# ggplot(pData(target_myData), aes(x = p53_liver, y = assayDataElement(target_myData["Anxa1", ], elt = "q_norm"))) +
#   geom_violin() +
#   geom_jitter(width = .2) +
#   labs(y = "- Expression") +
#   scale_y_continuous(trans = "log2") +
#   facet_wrap(~Strain) +
#   theme_bw()
#
# ggplot(pData(target_myData),
#        aes(x = assayDataElement(target_myData["Anxa1", ],
#                                 elt = "q_norm"),
#            y = assayDataElement(target_myData["Msln", ],
#                                 elt = "q_norm"),
#            color = p53_liver, label=class)) +
#   geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
#   geom_point(size = 3) +
#   theme_bw() +
#   scale_x_continuous(trans = "log2") +
#   scale_y_continuous(trans = "log2") +
#   labs(x = "Anxa1 Expression", y = "Msln Expression")

#
# ## ----heatmap, eval = TRUE, fig.width = 8, fig.height = 6.5, fig.wide = TRUE----
# # select top significant genes based on significance, plot with pheatmap
# GOI <- unique(subset(results, `FDR` < 0.01)$Gene)
# pheatmap(log2(assayDataElement(target_myData[GOI, ], elt = "q_norm")),
#          scale = "row",
#          show_rownames = FALSE, show_colnames = FALSE,
#          border_color = NA,
#          #clustering_method = "average",
#          #clustering_distance_rows = "correlation",
#          #clustering_distance_cols = "correlation",
#          cutree_cols = 4,
#          cutree_rows = 4,
#          breaks = seq(-3, 3, 0.05),
#          color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
#          annotation_col = pData(target_myData)[, c("class", "Strain")])

## ----maPlot, fig.width = 8, fig.height = 12, fig.wide = TRUE, warning = FALSE, message = FALSE----
results$MeanExp <-  rowMeans(assayDataElement(target_myData, elt = "q_norm"))

top_g2 <- results$Gene[results$Gene %in% top_g &
                         results$FDR < 0.01 &
                         abs(results$Estimate) > .5 &
                         results$MeanExp > quantile(results$MeanExp, 0.9)]

ggplot(subset(results, !Gene %in% neg_probes),
       aes(x = MeanExp, y = Estimate,
           size = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
  scale_x_continuous(trans = "log2") +
  geom_point(alpha = 0.5) +
  labs(y = "Enriched in Mucosa <- log2(FC) -> Enriched in Glands",
       x = "Mean Expression",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray")) +
  geom_text_repel(data = subset(results, Gene %in% top_g),
  #geom_text_repel(data = subset(results, Gene %in% top_g), size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2, max.overlaps = 75) +
  theme_bw(base_size = 16) +
  facet_wrap(~Subset, nrow = 2, ncol = 1)







###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###


results <- Null

#### FORMAT DATA
head(results)
names(results)[1] <- 'SYMBOL'
eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)

universe <- distinct(results, SYMBOL, .keep_all = T)
head(universe)

resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5 & results$`Pr(>|t|)` < 0.05)
#resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5 & results$FDR < 0.05)
summary(resultsGO)


ego <- enrichGO(gene          = resultsGO$ENTREZID,
                keyType       = "ENTREZID",
                universe      = universe$ENTREZID, ##list of all genes??
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", #"BP", "MF", "CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

goplot(ego2)
dotplot(ego)
upsetplot(ego)
#plotGOgraph(ego, useFullNames = T, useInfo = "names")


wcdf<-read.table(text=ego$GeneRatio, sep = "/")[1]
wcdf$term<-ego[,2]
wcdf$p.adjust<-ego$p.adjust
wcdf <- dplyr::filter(wcdf, V1 > 20)
wcdf <- dplyr::top_n(wcdf, 100, V1)
wcdf

selected_pathways <- c("synapse organization",
                       "axon guidance",
                       "regulation of synapse assembly",
                       "regulation of trans-synaptic signaling",
                       "positive regulation of synapse assembly",
                       "regulation of membrane potential",
                       "transmission of nerve impulse",
                       "axon development",
                       "synapse assembly",
                       "axonogenesis",
                       "chemical synaptic transmission, postsynaptic",
                       "synaptogenesis",
                       "gliogenesis",
                       "axonogenesis",
                       "cell-substrate adhesion",
                       "oligodendrocyte development",
                       "neurogenesis",
                       "cell-substrate adhesion",
                       "regulation of actin cytoskeleton organization",
                       "regulation of neurotransmitter transport")
dotplot(ego, showCategory = selected_pathways, font.size=10)


edox <- setReadable(ego, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)



#Create geneList
head(results)
results <- distinct(results, SYMBOL, .keep_all = T)
geneList = results$Estimate1
names(geneList) = as.character(results$ENTREZID)
geneList = sort(geneList, decreasing = T)

###gseGO PATHWAY###
ego2 <- gseGO(geneList      = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "MF", #"BP", "MF", and "CC"
              minGSSize    = 10,
              maxGSSize    = 1000,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
goplot(ego2)
dotplot(ego2)
upsetplot(ego2, 10)




wcdf<-read.table(text=ego2$ID, sep = "/")[1]
wcdf$term<-ego2[,2]
wcdf$p.adjust<-ego2$p.adjust
wcdf <- dplyr::filter(wcdf, V1 > 20)
wcdf <- dplyr::top_n(wcdf, 100, V1)
wcdf



















