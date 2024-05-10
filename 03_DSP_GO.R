#library(knitr)
#library(dplyr)
#library(ggforce)
#library(GeoMxWorkflows)
#library(NanoStringNCTools)
#library(GeomxTools)
library(readxl)
library(enrichplot)
#library(data.table)
#library(fgsea)
library(ggplot2)
#library(ggrepel)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#library(AnnotationHub)
library(GOSemSim)
library(clusterProfiler)
#library(cowplot)
#library(ReactomePA)
library(DOSE)
#library(msigdbr)
#####
#####GO
#####



results <- read.csv("F:/GeoMX KPC/Cheng_WTA1/processed_data/Seurat DEG/epithelium.csv")
results <- read_excel("C:/Users/edmondsonef/Desktop/Pathology Reports/ChengS/ChengS_DSP Results/Cheng DEG list.xlsx",
                      sheet = "All Epithelium")
                      #sheet = "All Stroma")
                      #sheet = "Epithelium <5mo")
                      #sheet = "Epithelium >5mo")
                      #sheet = "Stroma <5mo")
                      #sheet = "Stroma >5mo")




head(results)
names(results)[1] <- 'SYMBOL'
#names(results)[6] <- 'Pr(>|t|)'
head(results)
eg <- bitr(results$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
results <- dplyr::left_join(results, eg, by = "SYMBOL")
rm(eg)

universe <- read.csv("C:/Users/edmondsonef/Desktop/Pathology Reports/ChengS/ChengS_DSP Results/DEG/DEG_universe.csv")
universe <- dplyr::select(universe,SYMBOL,ENTREZID)
head(universe)
universe <- distinct(universe, SYMBOL, .keep_all = T)
head(universe)

# results_up <- results %>%  dplyr::filter(avg_log2FC >0)
# results_down <- results %>%  dplyr::filter(avg_log2FC <0)

gene <- results
gene <- dplyr::filter(results, results$avg_log2FC > 0.5)
gene <- dplyr::filter(results, results$avg_log2FC < 0.5)
head(gene)

ego <- enrichGO(gene          = gene$ENTREZID,
                keyType       = "ENTREZID",
                universe      = as.character(universe$ENTREZID), ##list of all genes??
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", #"BP", "MF", and "CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
p1 <- dotplot(ego, showCategory=20) + ggtitle("dotplot for ORA")
p1

write.csv(ego@result,"C:/Users/edmondsonef/Desktop/writeSTR.csv")


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("EGO_epithelium_combined.tiff", units="in", width=8, height=8, res=250)
p1
dev.off()

selected_pathways <- c("epithelial cell proliferation",
                       "negative regulation of epithelial cell proilferation",
                       "response to corticosteroid",
                       "response to glucocorticoid",
                       "response to steroid hormone",
                       "female pregnancy",
                       "multi-multicellular organism process",
                       "negative regulation of keratinocyte proliferation",
                       "keratinocyte proliferation",
                       "regulation of keratinocyte proliferatoin",
                       "thyroid hormone receptor activity",
                       "thyroid hormone receptor binding",
                       "cellular response to lipid",
                       "embryo implantation",
                       "reproductive structure development",
                       "reproductive system development",
                       "epidermis development",
                       "epidermal cell differentiation",
                       "leukocyte migration",
                       "placenta development",
                       "tumor necrosis factor production",
                       "decidualization",
                       "maternal process involved in female pregnancy",
                       "keratinization",
                       "epidermis development",
                       "skin development",
                       "keratinocyte differentiation",
                       "epidermal cell differentiation",
                       "negative regulation of keratinocyte proliferation",
                       "regulation of cell-cell adhesion",
                       "adaptive immune response",
                       "gland development",
                       "epithelial tube morphogenesis",
                       "proximal/distal pattern formation",
                       "connective tissue development",
                       "urogenital system development")

selected_pathways <- c("positive regulation of inflammatory response",
                       "positive regulation of cell activation",
                       "positive regulation of defense response",
                       "regulation of macrophage activation",
                       "positive regulation of response to external stimulus",
                       "macrophage activation",
                       "positive regulation of tumor necrosis factor production",
                       "positive regulation of leukocyte activation",
                       "nitric-oxide synthase biosynthetic process",
                       "tumor necrosis factor production",
                       "positive regulation of cytokine production",
                       "positive regulation of cell activation",
                       "positive regulation of defense response",
                       "regulation of macrophage activation",
                       "regulation of chemokine production",
                       "adaptive immune response",
                       "regulation of tumor necrosis factor production",
                       "immunoglobulin production",
                       "positive regulation of cytokine production",
                       "leukocyte migration",
                       "cytokine-mediated signaling pathway")


p1 <- dotplot(ego, showCategory=selected_pathways) + ggtitle("ORA: Il-33 pathways upregulated in Thra1PV")
p1

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("EGO_epithelium_combined.tiff", units="in", width=8, height=8, res=250)
p1
dev.off()




gene <- distinct(results, SYMBOL, .keep_all = T)
gene <- as.matrix(gene)
geneList = as.numeric(gene[,3])
#geneList = as.numeric(gene[,3]*-1)
names(geneList) = as.character(gene[,7])
geneList = sort(geneList, decreasing = T)
head(geneList)


#"BP" = biological process
#"MF" = molecular function
#"CC" = cellular component

ego2 <- gseGO(geneList     = geneList, ##??
              OrgDb        = org.Mm.eg.db,
              ont          = "BP", #"BP", "MF", and "CC"
              minGSSize    = 5,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = T)
head(ego2)
p2 <- dotplot(ego2, showCategory=20) + ggtitle("dotplot for GSEA")
p2


ego2 <- setReadable(ego2, 'org.Mm.eg.db', 'ENTREZID')
write.csv(ego2@result,"C:/Users/edmondsonef/Desktop/writeSTR.csv")


p2 <- dotplot(ego2, showCategory=20) + ggtitle("dotplot for GSEA")
fig <- cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])
fig
setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("fig.tiff", units="in", width=15, height=10, res=250)
fig
dev.off()




goplot(ego2)

edox <- setReadable(ego2, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox)
enrichplot::cnetplot.enrichResult(edox)

edox2 <- pairwise_termsim(edox)
enrichplot::cnetplot(edox2)
gseaplot(ego2, geneSetID = 1, title = ego2$Description[1])



heatplot(edox2, foldChange=geneList, showCategory=5)
treeplot(edox2, hclust_method = "average")
emapplot(edox2, layout="kk")
ridgeplot(edox2)
gseaplot(edox2, geneSetID = 337, by = "runningScore", title = edox2$Description[337])
gseaplot2(edox2, geneSetID = 49, title = edox2$Description[49])

dotplot(edox2, showCategory=20)









#
#
# #####
# ##### CompareCluster()
# #####
results <- read_excel("C:/Users/edmondsonef/Desktop/Pathology Reports/ChengS/ChengS_DSP Results/Cheng DEG list.xlsx",sheet = "All Epithelium")
#sheet = "All Stroma")
#sheet = "Epithelium <5mo")
#sheet = "Epithelium >5mo")
#sheet = "Stroma <5mo")
#sheet = "Stroma >5mo")

#load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_geoMX_new.RData")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_no.int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.06.22_comps_MHL_WITH.int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/07.08.22_class_MHL_no_int.csv")
#results <- read.csv("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/GENELIST_11-2-22_MHL_class_with_int.csv")
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
#resultsGO <- dplyr::filter(results, results$Estimate < -0.5 & results$FDR < 0.05)
#resultsGO <- dplyr::filter(results, results$Estimate < -0.5 & results$'Pr(>|t|)' < 0.05)

resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5 & results$FDR < 0.1)
#resultsGO <- dplyr::filter(results, abs(results$Estimate) > 1.0 & results$'Pr(>|t|)' < 0.05)
#resultsGO <- dplyr::filter(results, abs(results$Estimate) > 0.5 & results$FDR < 0.05)
#resultsGO <- dplyr::filter(results, abs(results$Estimate) > 1.0 & results$FDR < 0.01)
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

# resultsCC <- dplyr::filter(resultsCC, Contrast == c("1-Normal acini - 4-PanINlo",
#                                                   "1-Normal acini - 5-PanINhi",
#                                                   "5-PanINhi - 6-PDAC",
#                                                   "4-PanINlo - 6-PDAC"))


# resultsCC <- dplyr::filter(resultsCC, Contrast == c("ADM", "Bystander",
#                                                   "Acinar",
#                                                   "PanIN",
#                                                   "Carcinoma",
#                                                   "Islet",
#                                                   "IPMN",
#                                                   "Stroma"))
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




