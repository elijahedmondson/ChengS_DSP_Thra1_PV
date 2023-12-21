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
library(ggplot2)
library(cowplot)
library(ReactomePA)
library(DOSE)
library(msigdbr)
#####
#####GO
#####



#load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_geoMX_new.RData")
results <- read.csv("F:/GeoMX KPC/Cheng_WTA1/processed_data/Seurat/PVglands_Nglands.csv")
#results <- read.csv("F:/GeoMX KPC/Cheng_WTA1/processed_data/Seurat/PVmucosa_Nmucosa.csv")
#results <- read.csv("F:/GeoMX KPC/Cheng_WTA1/processed_data/Seurat/PVstroma_Nstroma.csv")
results <- read.csv("F:/GeoMX KPC/Cheng_WTA1/processed_data/DEG_12-21-23.csv")

head(results)
names(results)[2] <- 'SYMBOL'
names(results)[6] <- 'Pr(>|t|)'
head(results)

eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)
universe <- distinct(results, SYMBOL, .keep_all = T)
head(universe)


resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5 & results$FDR < 0.05)
summary(resultsGO)

gene <- resultsGO
head(gene)

ego <- enrichGO(gene          = gene$ENTREZID,
                keyType       = "ENTREZID",
                universe      = universe$ENTREZID, ##list of all genes??
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", #"BP", "MF", and "CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
p1 <- dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")
p1
gene <- distinct(results, SYMBOL, .keep_all = T)
head(gene)
## assume 1st column is ID
## 2nd column is FC
## feature 1: numeric vector
geneList = gene[,4] #which column?
head(geneList)

names(geneList) = as.character(gene[,9])
head(geneList)
geneList = sort(geneList, decreasing = T)
head(geneList)

#"BP" = biological process
#"MF" = molecular function
#"CC" = cellular component

ego2 <- gseGO(geneList     = geneList, ##??
              OrgDb        = org.Mm.eg.db,
              ont          = "BP", #"BP", "MF", and "CC"
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = T)


p2 <- dotplot(ego2, showCategory=30) + ggtitle("dotplot for GSEA")
fig <- cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])
fig
setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("fig.tiff", units="in", width=18, height=18, res=250)
fig
dev.off()






edox <- setReadable(ego2, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox, categorySize="pvalue",foldChange=geneList)

edox2 <- pairwise_termsim(edox)
enrichplot::cnetplot(edox2)
gseaplot(ego2, geneSetID = 1, title = ego2$Description[1])




#### FORMAT DATA

results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/11-2-22_MHL_class_with_int.csv")

results <- dplyr::filter(results, abs(results$Estimate) > 0.5)

head(results)
names(results)[6] <- 'Pr(>|t|)'
head(results)

results <- dplyr::filter(results, results$'Pr(>|t|)' < 0.05)

names(results)[2] <- 'SYMBOL'
eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)
universe <- distinct(results, SYMBOL, .keep_all = T)



results.up <- dplyr::filter(results, results$Estimate > 0)
results.down <- dplyr::filter(results, results$Estimate < 0)
names(results)
head(results)
results_list = split(results, f = results$Contrast)
names(results_list)

results.down_list = split(results.down, f = results.down$Contrast)
names(results.down_list)

results.up_list = split(results.up, f = results.up$Contrast)
#results <- distinct(results, SYMBOL, .keep_all = T)

write.csv(genelist, "C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/GENELIST_11-2-22_MHL_class_with_int.csv")
#results.up1 <- read.csv("C:/Users/edmondsonef/Desktop/results.up.csv")

genelist <- dplyr::bind_rows(results.up1, results.down1)




#
#
# #####
# ##### CompareCluster()
# #####
#
# #load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_geoMX_new.RData")
# #results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_no.int.csv")
# #results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_WITH.int.csv")
# #results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.08.22_class_MHL_no_int.csv")
# #results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/GENELIST_11-2-22_MHL_class_with_int.csv")
# results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/GENE LIST 07.06.22_comps_MHL_WITH.int.csv")
# #results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/11-3-22_MHL_progression2_with_int.csv")
#
# head(results)
# #names(results)[2] <- 'SYMBOL'
# names(results)[8] <- 'Pr(>|t|)'
# head(results)
#
# eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"), #toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
#            OrgDb="org.Mm.eg.db")
# results <- dplyr::left_join(results, eg, by = "SYMBOL")
# rm(eg)
# universe <- distinct(results, SYMBOL, .keep_all = T)
# head(universe)
#
# resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5)# & results$FDR < 0.1)
# #resultsGO <- dplyr::filter(results, results$Estimate > 0.5 & results$'Pr(>|t|)' < 0.05)
# #resultsGO <- dplyr::filter(results, results$Estimate < -0.5 & results$FDR < 0.05)
# #resultsGO <- dplyr::filter(results, results$Estimate < -0.5 & results$'Pr(>|t|)' < 0.05)
#
# resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5 & results$FDR < 0.1)
# #resultsGO <- dplyr::filter(results, abs(results$Estimate) > 1.0 & results$'Pr(>|t|)' < 0.05)
# #resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5 & results$FDR < 0.05)
# #resultsGO <- dplyr::filter(results, abs(results$Estimate) > 1.0 & results$FDR < 0.01)
# summary(resultsGO)
#
#
# # ###Create files for up and down regulated
# # resultsGO.up <- dplyr::filter(results, results$Estimate > 0.5)
# # resultsGO.down <- dplyr::filter(results, results$Estimate < -0.5)
# # write.csv(resultsGO.up,"C:/Users/edmondsonef/Desktop/resultsGO.up.csv")
# #
# # resultsGO.up <- read.csv("C:/Users/edmondsonef/Desktop/resultsGO.up.csv")
# # unique(resultsGO.up$Contrast)
# # resultsGO.down <- read.csv("C:/Users/edmondsonef/Desktop/resultsGO.down.csv")
# # unique(resultsGO.down$Contrast)
# # finalGL <- dplyr::bind_rows(resultsGO.down, resultsGO.up)
# #
# # head(finalGL)
# # names(finalGL)[7] <- 'Pr(>|t|)'
# # write.csv(finalGL, "C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/GENE LIST 07.06.22_comps_MHL_WITH.int.csv")
# #
#
# #gcSample =  list of different samples
# resultsCC <- dplyr::select(resultsGO, Estimate, Contrast, ENTREZID)
# unique(resultsCC$Contrast)
#
# # resultsCC <- dplyr::filter(resultsCC, Contrast == c("1-Normal acini - 4-PanINlo",
# #                                                   "1-Normal acini - 5-PanINhi",
# #                                                   "5-PanINhi - 6-PDAC",
# #                                                   "4-PanINlo - 6-PDAC"))
#
#
# # resultsCC <- dplyr::filter(resultsCC, Contrast == c("ADM", "Bystander",
# #                                                   "Acinar",
# #                                                   "PanIN",
# #                                                   "Carcinoma",
# #                                                   "Islet",
# #                                                   "IPMN",
# #                                                   "Stroma"))
# resultsCC <- dplyr::filter(resultsCC, Contrast == c("ADM",
#                                                     "Normal acini",
#                                                     "PanIN",
#                                                     "PDAC"))
#
# #Create a list containing gene IDs for each category
# ids <- unique(resultsCC$Contrast)
# mt_list<-list()
# for(i in 1:length(ids)){
#   id <- ids[i]
#   df <- dplyr::filter(resultsCC, resultsCC$Contrast == id)
#   mt_list[[i]]<-  as.character(df$ENTREZID)
# }
# ###CHECK THAT ROW NAMES ARE ACCURATE
# names(mt_list) <- c(ids)
# str(mt_list)
# mmu_kegg = download_KEGG(species = 'mmu', keggType = "KEGG", keyType = "kegg")
# ck <- compareCluster(geneCluster = mt_list,
#                      #fun = "enrichKEGG", organism = "mmu")
#                      #"groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
#                      fun = "enrichGO", ont = "BP", OrgDb = org.Mm.eg.db)
# ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
# ck <- pairwise_termsim(ck)
# ck2 <- simplify(ck, cutoff=0.7, by="p.adjust", select_fun=min)
#
# dotplot(ck2)
# head(ck2)
# goplot(ck2)
# cnetplot(ck2)
#
# #Create geneList
# head(resultsGO)
# gene <- distinct(resultsGO, SYMBOL, .keep_all = T)
# geneList = gene$Estimate
# names(geneList) = as.character(gene$ENTREZID)
# geneList = sort(geneList, decreasing = T)
#
#
#
#
# set.seed(2022-11-2)
# selected_pathways <- c("synapse organization",
#                        "synaptogenesis",
#                        #"alpha-amino acid metabolic process",
#                        "gliogenesis",
#                        "axonogenesis",
#                        "cell-substrate adhesion",
#                        "oligodendrocyte development",
#                        "neurogenesis",
#                        "actin filament organization",
#                        "regulation of translation",
#                        "cell-substrate adhesion",
#                        "regulation of actin cytoskeleton organization",
#                        "negative regulation of neural precursor cell")
# dotplot(ck2, showCategory = selected_pathways, font.size=10)
#
#
# cnetplot(ck2, node_label="gene",showCategory = selected_pathways,
#          cex_label_category = 1.2)
# cnetplot(ck2, node_label="all", showCategory = selected_pathways,
#          cex_label_category = 3.2,cex_label_gene = 1.7)#, foldChange=geneList)
#
#
#
#
# emapplot(ck, legend_n=2)
# emapplot(ck, pie="count", cex_category=1.5, layout="kk")
#
# cnetplot(ck, circular = T, colorEdge = TRUE, showCategory = selected_pathways)
#
#
# p1 <- treeplot(ck2)
# p2 <- treeplot(ck2, hclust_method = "average", showCategory = selected_pathways)
# aplot::plot_list(p1, p2, tag_levels='A')
#
# cnet <- cnetplot(ck2, node_label="all", showCategory = selected_pathways,
#                  cex_label_category = 2.5,cex_label_gene = 1.3, foldChange=geneList)
# cnet
#
# ggsave(cnet, file="C:/Users/edmondsonef/Desktop/cnet.png", width = 15, height = 15, units = "in", bg = "white")
#
#
#
#
# head(resultsCC)
# mt_list = split(resultsCC, f = resultsCC$Contrast)
#
# mydf <- data.frame(Entrez=names(geneList), FC=geneList)
# mydf <- mydf[abs(mydf$FC) > 1,]
# mydf$group <- "upregulated"
# mydf$group[mydf$FC < 0] <- "downregulated"
# mydf$othergroup <- "A"
# mydf$othergroup[abs(mydf$FC) > 2] <- "B"
#
# formula_res <- compareCluster(Entrez~class, data=mydf, fun="enrichKEGG")
#
# head(formula_res)
#



