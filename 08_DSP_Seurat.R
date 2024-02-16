library(GeomxTools)
library(Seurat)
library(SpatialDecon)
library(patchwork)
library(DESeq2)
library(MAST)


# #load("F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_12_21_2023.RData")
# #load("F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_1_3_2024.RData")
# load("F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_2_16_2024.RData")
# assayDataElementNames(target_myData)
#
# mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "q_norm")
# mySeurat
# save(mySeurat, final, target_myData, as.Seurat.NanoStringGeoMxSet, file="F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_seurat3.RData")

load("F:/GeoMX KPC/Cheng_WTA1/processed_data/Cheng_WTA1_seurat3.RData")

# head(mySeurat, 3)
# mySeurat@misc[1:8]
# head(mySeurat@misc$sequencingMetrics)# sequencing metrics
# head(mySeurat@misc$QCMetrics$QCFlags) # QC metrics
# head(mySeurat@assays$GeoMx@meta.features) # gene metadata
# VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)

mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "q_norm", ident = "COMP1")
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)
mySeurat <- as.Seurat.NanoStringGeoMxSet(target_myData, normData = "exprs", ident = "COMP1")
VlnPlot(mySeurat, features = "nCount_GeoMx", pt.size = 5)




mySeurat <- FindVariableFeatures(mySeurat, ident = "COMP5_age")
top50 <- head(VariableFeatures(mySeurat), 100)
LabelPoints(plot = plot1, points = top50, repel = T)

features <- c("Egln3", "Lrp2", "Slc34a2","Arg1","Hdc","Cxcl5","Ccl2",
              "Calb1","Csf3","Rnf186","Sgk1","Lrg1","Postn","Cxcl10",
              "Thrsp","Ier3", "Hr", "Klf9", "Dio3", "Shh",
              "Aldh1a1","Aldh1a3", "Adamtsl4","Htra1", "Epas1", "Fto", "Crmp1",
              "Pfkfb3", "Gbp3", "Gar22", "Desi1", "Trp53inp2","Stat5a","Aldoc",
              "Cxcl1","Il1a","Sprr2b","Spp1","Ppp1r1b","Icam1","Clca1","Iglc1",
              "Cxcl15","Aspg","Muc4","Spink12","Il17rb","Cfb","Cyp2f2",
              "Thra", "Thrab", "Il33","Wnt7a","Uox","Nppc","Cxcl2",
              "Slc26a7","Ncoa7","Trim15","Spink1", "Apobec2","Galm",
              "Agr2","Aldh1a3","Bcl3","Ccn3","Add2","Cyp21a1","Sprr2e","Pkdcc")
LabelPoints(plot = plot1, points = features, repel = T)


mySeurat <- ScaleData(mySeurat)

mySeurat <- RunPCA(mySeurat, assay = "GeoMx", verbose = FALSE)
#mySeurat <- RunPCA(mySeurat, features = VariableFeatures(object = mySeurat))
print(mySeurat[["pca"]], dims = 1:2, nfeatures = 3)
DimPlot(mySeurat, reduction = "pca", pt.size = 5, label = TRUE, group.by = "COMP5_age")


mySeurat <- FindNeighbors(mySeurat, reduction = "pca", dims = seq_len(30))

#mySeurat <- FindClusters(mySeurat, verbose = FALSE)

mySeurat <- RunUMAP(mySeurat, reduction = "pca", dims = seq_len(30))
DimPlot(mySeurat, reduction = "umap", pt.size = 5, label = TRUE, group.by = "COMP2")


levels(mySeurat)
levels(x = mySeurat) <- c("PV/N_glands_<5 months",
                          "PV/N_mucosa_<5 months",
                          "PV/N_stroma_<5 months",
                          "N/N_glands_<5 months",
                          "N/N_mucosa_<5 months",
                          "N/N_stroma_<5 months",
                          "PV/N_glands_>5 months",
                          "PV/N_mucosa_>5 months",
                          "PV/N_stroma_>5 months",
                          "N/N_glands_>5 months",
                          "N/N_mucosa_>5 months",
                          "N/N_stroma_>5 months")



levels(mySeurat)

features <- c("Klf9", "Dio3","Hr", "Thy1", "Sult1d1", "Lrp2")
features <- c("Klf4", "Prap1","Muc4", "Il33", "Agr2", "Wnt7a")

fig <- RidgePlot(mySeurat, sort = F, features = features,
                 idents = c("epithelium_PV/N",
                            "epithelium_N/N"),
                 ncol = 1)
fig


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("fig.tiff", units="in", width=20, height=15, res=300)
fig
dev.off()

###DEG
###DEG
###DEG
###DEG
###DEG
###DEG
###DEG
###DEG
###DEG
###DEG
###DEG


#NormlizeData() before "FindMarkers()
levels(mySeurat)
de_markers <- FindMarkers(mySeurat, ident.1 = "stroma_PV/N", ident.2 = "stroma_N/N",test.use = "negbinom")


results <- de_markers
library(tibble)
results <- tibble::rownames_to_column(results, "Gene")


library(ggrepel)
# Categorize Results based on P-value & FDR for plotting
results$Color <- "NS or FC < 0.5"
results$Color[results$p_val < 0.05] <- "P < 0.05"
results$Color[results$p_val_adj < 0.05] <- "FDR < 0.05"
results$Color[results$p_val_adj < 0.001] <- "FDR < 0.001"
results$Color[abs(results$avg_log2FC) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))


results$invert_P <- (-log10(results$p_val)) * sign(results$avg_log2FC)
top_g <- c()
top_g <- c(top_g,
           results[, 'Gene'][order(results[,'invert_P'], decreasing= T)[1:40]],
           results[, 'Gene'][order(results[,'invert_P'], decreasing= F)[1:40]])
top_g <- unique(top_g)
top_g


features <- c("Slc2a3","Hyal1","Lbp","Nexmif","Il17rb","Sult1d1","Znhit6",
              "Asph","Angptl7","St6galnac5","Pim3","Krt83","Gng12","Faim2",
              "C3","Prap1","Pglyrp1","Ltf","Lcn2","Fcgbp","Cfb","Gjb2", "Aoc1",
              "Trpv6", "Ppp1r1b","Egln3", "Lrp2", "Slc34a2","Arg1","Hdc","Cxcl5","Ccl2",
              "Calb1","Csf3","Rnf186","Sgk1","Lrg1","Postn","Cxcl10",
              "Thrsp","Ier3", "Hr", "Klf9", "Dio3", "Shh",
              "Aldh1a1","Aldh1a3", "Adamtsl4","Htra1", "Epas1", "Fto", "Crmp1",
              "Pfkfb3", "Gbp3", "Gar22", "Desi1", "Trp53inp2","Stat5a","Aldoc",
              "Cxcl1","Il1a","Sprr2b","Spp1","Ppp1r1b","Icam1","Clca1","Iglc1",
              "Cxcl15","Aspg","Muc4","Spink12","Il17rb","Cfb","Cyp2f2",
              "Thra", "Thrab", "Il33","Wnt7a","Uox","Nppc","Cxcl2",
              "Slc26a7","Ncoa7","Trim15","Spink1", "Apobec2","Galm",
              "Agr2","Aldh1a3","Bcl3","Ccn3","Add2","Cyp21a1","Sprr2e","Pkdcc", "Wnt7a", "Thy1",
              "Igha","Igkc","Ptn","Krt15","Chil1","Jchain","Krt5","Gas2l3","Krt14","Serpinb11","Aspg","Ctla2a",
              "Col15a1", "Sprr2d")


# Graph results
vplot <- ggplot(results,                                                             ###CHANGE
       aes(x = avg_log2FC, y = -log10(p_val),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "N/N <- log2(FC) -> PV/N",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% features & p_val < 0.05 & abs(avg_log2FC) >0.5),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
vplot

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("fig.tiff", units="in", width=8, height=8, res=200)
vplot
dev.off()


results <- dplyr::filter(results, p_val_adj < 0.05)
write.csv(results, "F:/GeoMX KPC/Cheng_WTA1/processed_data/Seurat DEG/stroma.csv")







levels(mySeurat)
de_ALLmarkers <- FindAllMarkers(mySeurat, grouping.var = "indent", verbose = T, test.use = "negbinom")


library(dplyr)

DEG_all <- de_ALLmarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
write.csv(DEG_all,"C:/Users/edmondsonef/Desktop/de_ALLmarkers.csv")

top10 <- DEG_all %>% group_by(cluster) %>% dplyr::top_n(20, avg_log2FC)
DefaultAssay(mySeurat) <- 'integrated'
S2 <- ScaleData(mySeurat, features=rownames(mySeurat))
DefaultAssay(S2) <- 'integrated'

idents.1 <- WhichCells(S2, idents = c("epithelium_PV/N",
                                      "stroma_PV/N",
                                      "epithelium_N/N",
                                      "stroma_N/N"))

idents.2 <- WhichCells(S2, idents = c("PV/N_glands_5-6months",
                                      "PV/N_mucosa_5-6months",
                                      "PV/N_stroma_5-6months",
                                      "N/N_glands_5-6months",
                                      "N/N_mucosa_5-6months",
                                      "N/N_stroma_5-6months"))
#"N/N_glands_5-6months",
#"N/N_mucosa_5-6months",
#"N/N_stroma_5-6months"))
hm1 <- DoHeatmap(S2, cells = idents.1)
hm2 <- DoHeatmap(S2, cells = cells.2)

fig <- DoHeatmap(object = S2, cells = idents.1, features = top10, group.by = "COMP1", group.bar = TRUE, group.colors = NULL,
                 disp.min = -2.5, disp.max = NULL, slot = "scale.data", assay = NULL, label = T,
                 size = 5.5, hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE, lines.width = NULL,
                 group.bar.height = 0.02, combine = TRUE)
fig

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("fig.tiff", units="in", width=20, height=22, res=300)
fig
dev.off()




mySeurat <- RunPCA(mySeurat, features = VariableFeatures(object = mySeurat))
print(mySeurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mySeurat, dims = 1:5, reduction = "pca")

#
# test.use =
# "wilcox" "bimod" "roc" "t" "poisson" "LR" "MAST"
# "negbinom" -- appropriate for count
# "DESeq2" -- appropriate for count (need to initialize the mySeurat with count matrix)
#
# Identify genes in the top 3rd of the CV values

mySeurat <- FindVariableFeatures(mySeurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mySeurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mySeurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# I partly figured this out and thought I would update in case someone else comes looking with the same issue.It seems that I was getting crazy LogFC's and weird looking volcanoes because of a problem with the normalization of my data. For these analyses I had used the Seurat "SCT integration" pipeline as outlined here: https://satijalab.org/seurat/v3.2/integration.html  I probably did something wrong, but after triple checking the workflow, I couldn't find it. When I switched back to an integration workflow that includes log normalization, rather than SCT normalization, everything appears to be fixed. So in the end, I am unsure of why the SCT integration pipeline failed on me, but switching to log normalization was the solution to my problem.

















###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###




#### FORMAT DATA
names(results)[1] <- 'SYMBOL'
eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)

universe <- distinct(results, SYMBOL, .keep_all = T)
head(universe)

resultsGO <- dplyr::filter(results, abs(results$avg_log2FC) > 0.5 & results$p_val_adj < 0.05)
summary(resultsGO)


ego <- enrichGO(gene          = resultsGO$ENTREZID,
                keyType       = "ENTREZID",
                universe      = universe$ENTREZID, ##list of all genes??
                OrgDb         = org.Mm.eg.db,
                ont           = "MF", #"BP", "MF", and "CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

dotplot(ego)
upsetplot(ego)
plotGOgraph(ego, useFullNames = T, useInfo = "names")

selected_pathways <- c("thyroid hormone receptor binding",
                       "thyroid hormone receptor activity",
                       "thyroid-stimulating hormone receptor activity",
                       "thyroid hormone receptor coactivator activity",
                       "thyroid hormone metabolic process",
                       "oligodendrocyte development",
                       "neurogenesis",
                       "cell-substrate adhesion",
                       "regulation of actin cytoskeleton organization")
dotplot(ego, showCategory = selected_pathways, font.size=10)

head(ego,10)

ggplot(ck2[1:20], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))

str(ego)

#####




#####
##### CompareCluster()
#####

#load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_geoMX_new.RData")

results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/GENE LIST 07.06.22_comps_MHL_WITH.int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/11-3-22_MHL_progression2_with_int.csv")

head(results)
#names(results)[2] <- 'SYMBOL'
names(results)[8] <- 'Pr(>|t|)'
head(results)

eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"), #toType=c("ENSEMBL", "ENTREZID", "UNIPROT"),
           OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)
universe <- distinct(results, SYMBOL, .keep_all = T)
head(universe)

resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5)# & results$FDR < 0.1)
#resultsGO <- dplyr::filter(results, results$Estimate > 0.5 & results$'Pr(>|t|)' < 0.05)

summary(resultsGO)


# ###Create files for up and down regulated
# resultsGO.up <- dplyr::filter(results, results$Estimate > 0.5)
# resultsGO.down <- dplyr::filter(results, results$Estimate < -0.5)
# write.csv(resultsGO.up,"C:/Users/edmondsonef/Desktop/resultsGO.up.csv")
#
# resultsGO.up <- read.csv("C:/Users/edmondsonef/Desktop/resultsGO.up.csv")
# unique(resultsGO.up$Contrast)
# resultsGO.down <- read.csv("C:/Users/edmondsonef/Desktop/resultsGO.down.csv")
# unique(resultsGO.down$Contrast)
# finalGL <- dplyr::bind_rows(resultsGO.down, resultsGO.up)
#
# head(finalGL)
# names(finalGL)[7] <- 'Pr(>|t|)'
# write.csv(finalGL, "C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/GENE LIST 07.06.22_comps_MHL_WITH.int.csv")
#

#gcSample =  list of different samples
resultsCC <- dplyr::select(resultsGO, Estimate, Contrast, ENTREZID)
unique(resultsCC$Contrast)
resultsCC <- dplyr::filter(resultsCC, Contrast == c("ADM",
                                                    "Normal acini",
                                                    "PanIN",
                                                    "PDAC"))

#Create a list containing gene IDs for each category
ids <- unique(resultsCC$Contrast)
mt_list<-list()
for(i in 1:length(ids)){
  id <- ids[i]
  df <- dplyr::filter(resultsCC, resultsCC$Contrast == id)
  mt_list[[i]]<-  as.character(df$ENTREZID)
}
###CHECK THAT ROW NAMES ARE ACCURATE
names(mt_list) <- c(ids)
str(mt_list)
mmu_kegg = download_KEGG(species = 'mmu', keggType = "KEGG", keyType = "kegg")
ck <- compareCluster(geneCluster = mt_list,
                     #fun = "enrichKEGG", organism = "mmu")
                     #"groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
                     fun = "enrichGO", ont = "BP", OrgDb = org.Mm.eg.db)
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
ck <- pairwise_termsim(ck)
ck2 <- simplify(ck, cutoff=0.7, by="p.adjust", select_fun=min)

dotplot(ck2)
head(ck2)
goplot(ck2)
cnetplot(ck2)

#Create geneList
head(resultsGO)
gene <- distinct(resultsGO, SYMBOL, .keep_all = T)
geneList = gene$Estimate
names(geneList) = as.character(gene$ENTREZID)
geneList = sort(geneList, decreasing = T)




set.seed(2022-11-2)
selected_pathways <- c("synapse organization",
                       "synaptogenesis",
                       #"alpha-amino acid metabolic process",
                       "gliogenesis",
                       "axonogenesis",
                       "cell-substrate adhesion",
                       "oligodendrocyte development",
                       "neurogenesis",
                       "actin filament organization",
                       "regulation of translation",
                       "cell-substrate adhesion",
                       "regulation of actin cytoskeleton organization",
                       "negative regulation of neural precursor cell")
dotplot(ck2, showCategory = selected_pathways, font.size=10)


cnetplot(ck2, node_label="gene",showCategory = selected_pathways,
         cex_label_category = 1.2)
cnetplot(ck2, node_label="all", showCategory = selected_pathways,
         cex_label_category = 3.2,cex_label_gene = 1.7)#, foldChange=geneList)




emapplot(ck, legend_n=2)
emapplot(ck, pie="count", cex_category=1.5, layout="kk")

cnetplot(ck, circular = T, colorEdge = TRUE, showCategory = selected_pathways)


p1 <- treeplot(ck2)
p2 <- treeplot(ck2, hclust_method = "average", showCategory = selected_pathways)
aplot::plot_list(p1, p2, tag_levels='A')

cnet <- cnetplot(ck2, node_label="all", showCategory = selected_pathways,
                 cex_label_category = 2.5,cex_label_gene = 1.3, foldChange=geneList)
cnet

ggsave(cnet, file="C:/Users/edmondsonef/Desktop/cnet.png", width = 15, height = 15, units = "in", bg = "white")




head(resultsCC)
mt_list = split(resultsCC, f = resultsCC$Contrast)

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~class, data=mydf, fun="enrichKEGG")

head(formula_res)



names(mt_list)
gene <- resultsCC
head(gene)


gene <- distinct(gene, SYMBOL, .keep_all = T)

#####gseGO()
head(gene)
## assume 1st column is ID
## 2nd column is FC
## feature 1: numeric vector
geneList = gene[,1] #which column?
head(geneList)

names(geneList) = as.character(gene[,3])
head(geneList)
geneList = sort(geneList, decreasing = T)
head(geneList)

#"BP" = biological process
#"MF" = molecular function
#"CC" = cellular component

ego3 <- gseGO(geneList     = geneList, ##??
              OrgDb        = org.Mm.eg.db,
              ont          = "BP", #"BP", "MF", and "CC"
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

#head(ego3)
goplot(ego3)
dotplot(ego3)
upsetplot(ego3, 10)







###
###
#GENE CONCEPT NETWORK
ego <- pairwise_termsim(ego)
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

selected_pathways <- c("synapse organization",
                       "axonogenesis",
                       "cell-substrate adhesion")

selected_pathways <- c("synapse organization",
                       "synaptogenesis",
                       "gliogenesis",
                       "axonogenesis",
                       "cell-substrate adhesion",
                       "regulation of neurogenesis",
                       "neurogenesis",
                       "cell-substrate adhesion",
                       "regulation of actin cytoskeleton organization")
dotplot(ego2, showCategory = selected_pathways, font.size=10)
cnetplot(ego2, node_label="all", categorySize="pvalue", showCategory = selected_pathways, foldChange=geneList)
cnetplot(edox, foldChange=geneList, circular = T, colorEdge = TRUE, showCategory = "axonogenesis")
heatplot(ego, foldChange=geneList, showCategory=selected_pathways)

library(clusterProfiler)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)
p1 <- emapplot(xx)
p2 <- emapplot(xx, legend_n=2)
p3 <- emapplot(xx, pie="count")
p4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
###
###














#####
#####PLOTTING GENES
#####
#####
## ----targetTable, eval = TRUE, as.is = TRUE-----------------------------------

#load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_DSP.RData")

head(gene)
names(gene)[1] <- 'Gene'
head(gene)

kable(subset(gene, Gene %in% c("Pdzd8", "Mtch2", "Spock3", "Serpina3k", "Cybrd1", "Vars2")),
      row.names = FALSE)
kable(subset(gene, Gene %in% c("Pdzd8")), row.names = FALSE)


## ----targetExprs, eval = TRUE-------------------------------------------------
# show expression for a single target: PDHA1
ggplot(pData(target_myData),
       aes(x = progression2, fill = progression2,
           y = assayDataElement(target_myData["Kras", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "") +
  scale_y_continuous(trans = "log2") +
  #facet_wrap(~Strain) +
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(), text = element_text(size = 24))


## ----targetExprs2, fig.width = 8, fig.wide = TRUE, eval = TRUE----------------
glom <- pData(target_myData)$progression1# == "Metastasis"

# show expression of PDHA1 vs ITGB1
ggplot(pData(target_myData),
       aes(x = assayDataElement(target_myData["Trp53", ],
                                elt = "q_norm"),
           y = assayDataElement(target_myData["Msln", ],
                                elt = "q_norm"),
           color = comps, label=dsxf)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  # geom_vline(xintercept =
  #              max(assayDataElement(target_myData["Net1", ],
  #                                   elt = "q_norm")),
  #            lty = "dashed", col = "darkgray") +
  # geom_hline(yintercept =
  #              max(assayDataElement(target_myData["Rock2", ],
  #                                   elt = "q_norm")),
  #            lty = "dashed", col = "darkgray") +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  labs(x = "Trp53 Expression", y = "Msln Expression")
#+
#facet_wrap(~class)

## ----heatmap, eval = TRUE, fig.width = 8, fig.height = 6.5, fig.wide = TRUE----
# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(gene, `FDR` < 0.001)$Gene)
pheatmap(log2(assayDataElement(target_myData[GOI, ], elt = "q_norm")),
         scale = "row",
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 3, cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_myData)[, c("progression1", "progression1")])

## ----maPlot, fig.width = 8, fig.height = 12, fig.wide = TRUE, warning = FALSE, message = FALSE----
gene$MeanExp <-
  rowMeans(assayDataElement(target_myData,
                            elt = "q_norm"))

top_g2 <- gene$Gene[gene$Gene %in% top_g &
                      gene$FDR < 0.001 &
                      abs(gene$Estimate) > .5 &
                      gene$MeanExp > quantile(gene$MeanExp, 0.9)]

ggplot(subset(gene, !Gene %in% neg_probes),
       aes(x = MeanExp, y = Estimate,
           size = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
  scale_x_continuous(trans = "log2") +
  geom_point(alpha = 0.5) +
  labs(y = "Enriched in XXX <- log2(FC) -> Enriched in XXX",
       x = "Mean Expression",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray")) +
  geom_text_repel(data = subset(gene, Gene %in% top_g2),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2) +
  theme_bw(base_size = 16) +
  facet_wrap(~Subset, nrow = 2, ncol = 1)







###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###
###PATHWAY###




#### FORMAT DATA
names(results)[1] <- 'SYMBOL'
eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),#, "UNIPROT","ENSEMBL", ),
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













































as.Seurat.NanoStringGeoMxSet <- function(x, ident = NULL, normData = NULL,
                                         coordinates = NULL,
                                         forceRaw = FALSE, ...){

  if (!try(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Please install Seurat from CRAN before converting to a Seurat object")
  }else{
    requireNamespace("Seurat", quietly = TRUE)
  }


  if(featureType(x) == "Probe"){
    stop("Data must be on Target level before converting to a Seurat Object")
  }

  if(is.null(normData)){
    stop("normData must not be NULL, please provide name of normalized counts matrix")
  }

  if(!normData %in% assayDataElementNames(x)){
    stop(paste0("The normData name \"", normData, "\" is not a valid assay name. Valid names are: ",
                paste(assayDataElementNames(x), collapse = ", ")))
  }

  normFactor_names <- "normFactors|qFactors|negFactors|hkFactors|hknormFactors"

  if(length(grep(pattern = normFactor_names, names(sData(x)))) == 0 &
     forceRaw == FALSE){
    stop("It is NOT recommended to use Seurat's normalization for GeoMx data.
             Normalize using GeomxTools::normalize() or set forceRaw to TRUE if you want to continue with Raw data")
  }


  sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID",
                         "Well", "SeqSetId", "Raw", "Trimmed", "Stitched",
                         "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads",
                         "NTC_ID", "NTC", "Trimmed (%)", "Stitched (%)",
                         "Aligned (%)", "Saturated (%)")

  QCMetrics <- "QCFlags"

  seuratConvert <- suppressWarnings(Seurat::CreateSeuratObject(counts = assayDataElement(x, normData),
                                                               assay = "GeoMx",
                                                               project = expinfo(experimentData(x))[["title"]]))
  seuratConvert <- suppressWarnings(Seurat::AddMetaData(object = seuratConvert,
                                                        metadata = sData(x)[,!colnames(sData(x)) %in%
                                                                              c(sequencingMetrics,
                                                                                QCMetrics)]))
  seuratConvert@assays$GeoMx <- Seurat::AddMetaData(object = seuratConvert@assays$GeoMx,
                                                    metadata = fData(x))

  if(!is.null(ident)){
    if(!ident %in% colnames(seuratConvert@meta.data)){
      stop(paste0("ident \"", ident, "\" not found in GeoMxSet Object"))
    }

    Seurat::Idents(seuratConvert) <- seuratConvert[[ident]]
  }


  seuratConvert@misc <- otherInfo(experimentData(x))
  seuratConvert@misc[["sequencingMetrics"]] <- sData(x)[colnames(sData(x)) %in%
                                                          sequencingMetrics]
  seuratConvert@misc[["QCMetrics"]] <- sData(x)[colnames(sData(x)) %in%
                                                  QCMetrics]

  if(ncol(seuratConvert@misc[["QCMetrics"]]) == 0){
    seuratConvert@misc[["QCMetrics"]] <- NULL
  }

  if(!is.null(coordinates)){
    xcoord <- coordinates[1]
    ycoord <- coordinates[2]

    if(xcoord %in% colnames(seuratConvert@meta.data) &
       ycoord %in% colnames(seuratConvert@meta.data)){
      coord.df <- data.frame(x=seuratConvert@meta.data[[xcoord]],
                             y=seuratConvert@meta.data[[ycoord]])
      colnames(coord.df) <- coordinates
      seuratConvert@meta.data <- seuratConvert@meta.data[!colnames(seuratConvert@meta.data) %in%
                                                           coordinates]
    }else{
      if(!xcoord %in% colnames(seuratConvert@meta.data) &
         !ycoord %in% colnames(seuratConvert@meta.data)){
        stop(paste0("xcoord \"", xcoord, "\" and ycoord \"",
                    ycoord, "\" not found in GeoMxSet Object"))
      }

      if(!xcoord %in% colnames(seuratConvert@meta.data)){
        stop(paste0("xcoord \"", xcoord,
                    "\" not found in GeoMxSet Object"))
      }

      if(!ycoord %in% colnames(seuratConvert@meta.data)){
        stop(paste0("ycoord \"", ycoord,
                    "\" not found in GeoMxSet Object"))
      }
    }

    rownames(coord.df) <- rownames(seuratConvert@meta.data)

    # need to create DSP specific image class
    seuratConvert@images$image =  new(
      Class = 'SlideSeq',
      assay = "GeoMx",
      key = "image_",
      coordinates = coord.df
    )
  }

  return(seuratConvert)
}

#' Convert Object to SpatialExperiment
#'
#' @param x GeoMxSet object to convert
#' @param ... arguments to be passed to other methods
#'
#' @export


