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

load("F:/GeoMX KPC/WTA_11232022/processed_data/KPC_geoMX_exp2.RData")

#test <- "Step1_KPC"
#test <- "Step2_KPC"
#test <- "Step3_KPC_allmet"
#test <- "Step3_KPC_lung"
#test <- "Step3_KPC_liver"

#test <- "Step3_ortho"

#test <- "Step1_R172H"
#test <- "Step2_R172H"
#test <- "Step3_R172H"

#test <- "Step1_R270H"
#test <- "Step2_R270H"
#test <- "Step3_R270H"

#test <- "p53_panin"
test <- "p53_PDAC"
#test <- "p53_liver"

# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$p53_PDAC, c("KPC_PDAC","R172H_PDAC", "R270H_PDAC"))                           
pData(target_myData)[["slide"]] <-                                            ### Control for 
  factor(pData(target_myData)[["MHL"]])
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("Full ROI")) {
  ind <- pData(target_myData)$segment == status
  mixedOutmc <-
    mixedModelDE(target_myData[, ind], elt = "log_q",
                 #modelFormula = ~ testRegion + (1 + testRegion | slide),        
                 modelFormula = ~ testRegion + (1 | slide),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- status
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}



results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color, levels = c("NS or FC < 0.5", "P < 0.05", "FDR < 0.05", "FDR < 0.001"))
dplyr::count(results, FDR < 0.05)
dplyr::count(results, `Pr(>|t|)` < 0.05)

top <- dplyr::filter(results, `Pr(>|t|)` < 0.05)
# top <- dplyr::filter(results, results$FDR < 0.05)
write.csv(top, "F:/GeoMX KPC/WTA_11232022/processed_data/STRAIN_p05_KPC_PDAC.csv")

head(results)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:50]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:50]])
}
top_g <- unique(top_g)
top_g




features <- c("Rock2","Cybrd1","Nr1d1","Bsg","Tmprss4","Tm9sf3","Mmp23","Rhof","Sftpd", "Aqp5","Ccna1",
              "Muc3","Muc5ac","Muc3a","Kif12","Calml4","Dbp", "Mrtfb", "Rplp0","Dnajc10","Rps12",
              "Pdzd8", "Mtch2", "Msln", "Prom1", "Vars2","Porcn","Rpl6","Ybx1","Wfdc2","Tpi1",
              "Golim4","Otop3","F3", "Id2","Adamtsl5","Bag1","Rnf186","Glis2","Slc35f5","Tspan12",
              "Slc9a4", "Ephb2", "Tmem45b","Tmprss2","Pdxdc1","Lgals2", "Esrp1", "Tmem54", "Ptprf", "Ccnd2",
              "Ern2","Sult1c2", "Gltp","Spock3","Sgms2","Rasgrf1","St8sia3",
              "Rap1gap","Rbms3","Ccdc92","Ncald","Ppp1r1b","Gabbr2","Alb",
              "Nt5c2","Cdkn2a","Atrnl1","Camk2n1","Sem1", "Ctnnd1","Adgre5", "Dennd4c",
              "Setbp1","Dennd4c","Hs3st1","Shf","Nfib", "Tuba1b", "Net1", "Ncald","Spock3",
              "Smad4", "Flna", "Cntn1", "Cntn6","Sgms2","Nrxn1","Nrxn2","Nrxn3","Lamb2","Rasgrf1",
              "Lama5", "Rtn4", "Picalm","Efnb2", "Rbms3", "Rock2","Ephb2","Efnb2", "Adam10", "Mmp2", "Mmp9", 
              "Rac1", "St8sia3", "Camk2n1", "Cdc42", "Spock3", "Rasgrf1",
              "Lama5", "Itgb1", "Ezr","S100a6", "Gsto1", "Gkn1","Ezr","Sema3d", "Sema4b","Sema4g","Sema5a","St8sia3",
              "Lypd8l", "Anxa2", "Cdh1", "Myrf", "Flna", "Slc12a2", "Actn1", "Fn1", "Hnf1b",
              "Vasp","Vdac2", "Syncrip", "Rpl5", "Pard3","Dync1i2", "Calm1", "Calm2", "Calm3", "Itgb1","Kras","Trp53","Net1","Nt5c2",
              "Clu","S100a6", "Anxa2", "Myrf", "Sema4b","Sema4g","Efnb2", "Flna", "Slc12a2", "Actn1", "Actb","Tuba1b",
              "Vasp", "Syncrip", "Pard3","Rac1", "Rhoa", "Cdc42", "Dync1i2", "Calm1", "Calm2", "Calm3","Lama5", "Itgb1")


#reverse log fold change to fit with label
results$Estimate1 <- results$Estimate*(-1)
# Graph results
volc_plot <- ggplot(results,                                                             ###CHANGE
       aes(x = Estimate1, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = " <- log2(FC) -> ",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g),# & FDR < 0.01),
  #geom_text_repel(data = subset(results, Gene %in% features & `Pr(>|t|)` < 0.05),
                  size = 6, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") 

volc_plot


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
results$MeanExp <-
  rowMeans(assayDataElement(target_myData,
                            elt = "q_norm"))

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
  labs(y = "Enriched in XXX <- log2(FC) -> Enriched in xxx",
       x = "Mean Expression",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray")) +
  geom_text_repel(data = subset(results, Gene %in% features),
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
names(results)[2] <- 'SYMBOL'
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
              ont          = "BP", #"BP", "MF", and "CC"
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

selected_pathways <- c("synapse organization",
                       "axon guidance",
                       "axon development",
                       "axonogenesis",
                       "synaptogenesis",
                       "gliogenesis",
                       "axonogenesis", 
                       "cell-substrate adhesion",
                       "oligodendrocyte development", 
                       "neurogenesis",
                       "cell-substrate adhesion",
                       "regulation of actin filament organization",
                       "anterograde trans-synaptic signaling",
                       "trans-synaptic signaling",
                       "synaptic signaling",
                       "regulation of dendrite extension",
                       "dendrite extension",
                       "regulation of trans-synaptic signaling")

dotplot(ego2, showCategory = selected_pathways, font.size=10)




















