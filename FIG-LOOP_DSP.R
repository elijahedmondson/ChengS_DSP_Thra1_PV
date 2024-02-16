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
#library(org.Hs.eg.db)
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

datadir <-"C:/Users/edmondsonef/Desktop/R-plots/"
setwd(datadir)

results <- read.csv("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Lists/DEG_1-3-24_withIntercept.csv")
head(results)
names(results)[2] <- 'SYMBOL'

universe <- read.csv("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/DEG_universe.csv")
universe <- dplyr::select(universe,SYMBOL,ENTREZID)

m_t2g <- "C:/Users/edmondsonef/Desktop/DSP GeoMx/data/WTA_04122022/raw_data/msigdb.v7.5.1.entrez.gmt"


eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)
names(results)[6] <- 'Pr(>|t|)'
mt_list = split(results, f = results$Contrast)
names(mt_list)



names(mt_list)


##FOR LOOP
i = 12
for(i in 8:15){
  suffix <- names(mt_list[i])
  outname <-paste0(suffix, "_NEW_comps_MHL_with_int")

  gene <- mt_list[[i]]
  gene <- distinct(gene, SYMBOL, .keep_all = T)
  #gene <- dplyr::filter(gene, gene$FDR < 0.05)

  top <- dplyr::filter(gene, gene$FDR < 0.05)
  #top <- dplyr::filter(gene, gene$FDR < 0.001)
  count = count(top)
  print(paste(suffix, ":",count, "genes with FDR < 0.05."))

  top_g <- c()
  for(cond in c("Full ROI")) {
    ind <- gene$Subset == cond
    top_g <- c(top_g,
               gene[ind, 'SYMBOL'][
                 order(gene[ind, 'invert_P'], decreasing = TRUE)[1:30]],
               gene[ind, 'SYMBOL'][
                 order(gene[ind, 'invert_P'], decreasing = FALSE)[1:30]])
  }
  top_g <- unique(top_g)
  top_g

  #reverse log fold change to fit with label
  gene$Estimate1 <- gene$Estimate*(-1)


  pVP <- ggplot(gene,
                aes(x = Estimate1, y = -log10(`Pr(>|t|)`),
                    color = Color, label = SYMBOL)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x =  suffix,#"<- log2(FC) ->",
         y = "Significance, -log10(P)",
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(gene, SYMBOL %in% top_g & `Pr(>|t|)` < 0.05),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom")
  pVP
  # volcano <- paste0(datadir, "_", outname, "_volcano.png")
  # ggsave(pVP, file=volcano, width = 8, height = 8, units = "in", bg = "white")
  # rm(volcano)

  ego <- enrichGO(gene          = gene$ENTREZID,
                  keyType       = "ENTREZID",
                  universe      = as.character(universe$ENTREZID),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "MF", #"BP", "MF", and "CC"
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

  pDP <- dotplot(ego)
  upsetplot(ego)


  #####gseGO()
  head(gene)
  gene <- distinct(gene, SYMBOL, .keep_all = T)
  geneList = gene$Estimate1
  names(geneList) = as.character(gene$ENTREZID)
  geneList = sort(geneList, decreasing = T)

  ego2 <- gseGO(geneList      = geneList,
                OrgDb        = org.Mm.eg.db,
                ont          = "MF", #"BP", "MF", and "CC"
                minGSSize    = 5,
                maxGSSize    = 1000,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
  p1 <- goplot(ego2)
  p2 <- dotplot(ego2)
  p3 <- upsetplot(ego2, 10)
  gg_all <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

  multiplot <- paste0(datadir, "_", outname, "_gseGO_BP.png")
  ggsave(gg_all, file=multiplot, width = 15, height = 10, units = "in", bg = "white")
  rm(ego2,gg_all, p1, p2, p3)


  ###MSigDb analysis
  ###MSigDb analysis
  ###MSigDb analysis


  # H: hallmark gene sets       --good
  # C1: positional gene sets
  # C2: curated gene sets       --good
  # C3: motif gene sets         --good
  # C4: computational gene sets
  # C5: GO gene sets            --good
  # C6: oncogenic signatures    --good
  # C7: immunologic signatures  --good
  msig_list <- c("C2", "C6", "C7", "H", "C5", "C3")

  j = 5
  for(j in 5:5){
    Msig <- msig_list[j]
    m_t2g <- msigdbr(species = "Mus musculus", category = Msig) %>%
      dplyr::select(gs_name, entrez_gene)

    edo <- enricher(gene$ENTREZID, TERM2GENE=m_t2g)
    edo2 <- GSEA(geneList, TERM2GENE = m_t2g, eps=0)

    p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
    p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
    gg_all <- cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])

    multiplot <- paste0(datadir, "_", outname, "_", Msig, "_GSEA_ORA_GSEA.png")
    ggsave(gg_all, file=multiplot, width = 15, height = 10, units = "in", bg = "white")
    rm(multiplot, gg_all, p1, p2)

}
}
