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

# universe <- read.csv("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/DEG_universe.csv")
# universe <- dplyr::select(universe,SYMBOL,ENTREZID)
# head(universe)
ALLepi <- read_excel("C:/Users/edmondsonef/Desktop/KLF9 manuscript/Supp. Table 1, DEG list.xlsx", sheet = "All Epithelium")
Allstr <- read_excel("C:/Users/edmondsonef/Desktop/KLF9 manuscript/Supp. Table 1, DEG list.xlsx", sheet = "All Stroma")
YOUNGepi <- read_excel("C:/Users/edmondsonef/Desktop/KLF9 manuscript/Supp. Table 1, DEG list.xlsx", sheet = "Epithelium <5mo")
OLDepi <- read_excel("C:/Users/edmondsonef/Desktop/KLF9 manuscript/Supp. Table 1, DEG list.xlsx", sheet = "Epithelium >5mo")
YOUNGstroma <- read_excel("C:/Users/edmondsonef/Desktop/KLF9 manuscript/Supp. Table 1, DEG list.xlsx", sheet = "Stroma <5mo")
OLDstroma <- read_excel("C:/Users/edmondsonef/Desktop/KLF9 manuscript/Supp. Table 1, DEG list.xlsx", sheet = "Stroma >5mo")

gene_list <- list(Epithelium = ALLepi$Gene,
                  Stroma = Allstr$Gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red"),
                            cat.col = c("blue", "red"),
                            cex = 4,lty = "blank",
                            cat.cex = 2)
cowplot::plot_grid(VennDiagram)


t=get.venn.partitions(gene_list, keep.elements = T, force.unique = T)


ALLepi_WT <- ALLepi %>%  dplyr::filter(avg_log2FC >0)
ALLepi_PV <- ALLepi %>%  dplyr::filter(avg_log2FC <0)
Allstr_WT <- Allstr %>%  dplyr::filter(avg_log2FC >0)
Allstr_PV <- Allstr %>%  dplyr::filter(avg_log2FC <0)

gene_list <- list(WT_Epithelium = ALLepi_WT$Gene,
                  WT_Stroma = Allstr_WT$Gene,
                  PV_Epithelium = ALLepi_PV$Gene,
                  PV_Stroma = Allstr_PV$Gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red","green","yellow"),
                            cat.col = c("blue", "red","green","yellow"),
                            cex = 2,lty = "blank",
                            cat.cex = 1.5)
cowplot::plot_grid(VennDiagram)

gene_list <- list(#WT_Epithelium = ALLepi_WT$Gene,
                  WT_Stroma = Allstr_WT$Gene,
                  #PV_Epithelium = ALLepi_PV$Gene)#,
                  PV_Stroma = Allstr_PV$Gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red"),
                            cat.col = c("blue", "red"),
                            cex = 4,lty = "blank",
                            cat.cex = 2)
cowplot::plot_grid(VennDiagram)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=10, height=10, res=200)
cowplot::plot_grid(VennDiagram)
dev.off()





gene_list <- list(All = Allstr$Gene,
                  Young = YOUNGstroma$Gene,
                  Old = OLDstroma$Gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red", "green"),
                            cat.col = c("blue", "red", "green"),
                            cex = 2,lty = "blank",
                            cat.cex = 2)
cowplot::plot_grid(VennDiagram)



library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

gene_list <- list(All_Stroma = Allstr$Gene,
                  Young_Stroma = YOUNGstroma$Gene,
                  Old_Stroma = OLDstroma$Gene)
                  #All_Epitheium = ALLepi$Gene,
                  #Young_Epitheium = YOUNGepi$Gene,
                  #Old_Epitheium = OLDepi$Gene)

VennDiagram <- venn.diagram(x = gene_list, filename = NULL,output=F,
                            fill = myCol,
                            cat.col = myCol,
                            cex = 2,lty = "blank",
                            cat.cex = 2)
cowplot::plot_grid(VennDiagram)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=10, height=10, res=200)
cowplot::plot_grid(VennDiagram)
dev.off()



ALLepi_up <- ALLepi %>%  dplyr::filter(avg_log2FC >0)
ALLepi_down <- ALLepi %>%  dplyr::filter(avg_log2FC <0)
Allstr_up <- Allstr %>%  dplyr::filter(avg_log2FC >0)
Allstr_down <- Allstr %>%  dplyr::filter(avg_log2FC <0)

gene_list <- list(WT_Epi = ALLepi_up$Gene,
                  PV_Epi = ALLepi_down$Gene,
                  WT_Str = Allstr_up$Gene,
                  PV_Str = Allstr_down$Gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red", "green", "yellow"),
                            cat.col = c("blue", "red", "green", "yellow"),
                            cex = 2,lty = "blank",
                            cat.cex = 2)
cowplot::plot_grid(VennDiagram)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=10, height=10, res=200)
cowplot::plot_grid(VennDiagram)
dev.off()



write <- dplyr::full_join(MMstr, SEURstr, by = "Gene")
write.csv(write,"C:/Users/edmondsonef/Desktop/writeSTR.csv")



N_epithelium <- de_ALLmarkers %>%  dplyr::filter(cluster == "epithelium_N/N" & p_val_adj < 0.05 & avg_log2FC >1)
PV_epithelium <- de_ALLmarkers %>%  dplyr::filter(cluster == "epithelium_PV/N" & p_val_adj < 0.05 & avg_log2FC >1)
N_stroma <- de_ALLmarkers %>%  dplyr::filter(cluster == "stroma_N/N" & p_val_adj < 0.01 & avg_log2FC >1)
PV_stroma <- de_ALLmarkers %>%  dplyr::filter(cluster == "stroma_PV/N" & p_val_adj < 0.01 & avg_log2FC >1)


gene_list <- list(PV = PV_epithelium$gene,
                  N = N_epithelium$gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red"),
                            cat.col = c("blue", "red"),
                            cex = 2,lty = "blank",
                            cat.cex = 2)

cowplot::plot_grid(VennDiagram)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=10, height=10, res=200)
cowplot::plot_grid(VennDiagram)
dev.off()

gene_list <- list(PV_epithelium = PV_epithelium$gene,
                  N_epithelium = N_epithelium$gene,
                  PV_stroma = PV_stroma$gene,
                  N_stroma = N_stroma$gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red", "green", "yellow"),
                            cat.col = c("blue", "red", "green", "yellow"),
                            cex = 2,lty = "blank",
                            cat.cex = 2)

cowplot::plot_grid(VennDiagram)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=10, height=10, res=200)
cowplot::plot_grid(VennDiagram)
dev.off()



N_epithelium_young <- de_ALLmarkers %>%  dplyr::filter(cluster == "epithelium_N/N_<5 months" & abs(avg_log2FC) >1)
N_epithelium_old <- de_ALLmarkers %>%  dplyr::filter(cluster == "epithelium_N/N_5+ months" & abs(avg_log2FC) >1)
PV_epithelium_young <- de_ALLmarkers %>%  dplyr::filter(cluster == "epithelium_PV/N_<5 months" & abs(avg_log2FC) >1)
PV_epithelium_old <- de_ALLmarkers %>%  dplyr::filter(cluster == "epithelium_PV/N_5+ months" & abs(avg_log2FC) >1)

N_stroma_young <- de_ALLmarkers %>%  dplyr::filter(cluster == "stroma_N/N_<5 months" & abs(avg_log2FC) >1)
N_stroma_old <- de_ALLmarkers %>%  dplyr::filter(cluster == "stroma_N/N_5+ months" & abs(avg_log2FC) >1)
PV_stroma_young <- de_ALLmarkers %>%  dplyr::filter(cluster == "stroma_PV/N_<5 months" & abs(avg_log2FC) >1)
PV_stroma_old <- de_ALLmarkers %>%  dplyr::filter(cluster == "stroma_PV/N_5+ months" & abs(avg_log2FC) >1)


gene_list <- list(PV_young = PV_stroma_young$gene,
                  PV_old = PV_stroma_old$gene,
                  N_young = N_stroma_young$gene,
                  N_old = N_stroma_old$gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red", "green", "yellow"),
                            cat.col = c("blue", "red", "green", "yellow"),
                            cex = 2,lty = "blank",
                            cat.cex = 2)

cowplot::plot_grid(VennDiagram)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=10, height=10, res=200)
cowplot::plot_grid(VennDiagram)
dev.off()




seurat <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list seurat.xlsx", sheet = "PV_glands_seurat")
mmDE <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list mixedmodelDE.xlsx", sheet = "PV_glands")

PV_glands <- dplyr::full_join(seurat, mmDE, by = "Gene")
PV_glands <- distinct(PV_glands, Gene, .keep_all = T)
write.csv(PV_glands, "C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/FULL_PV_glands.csv")

rm(seurat, mmDE)
seurat <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list seurat.xlsx", sheet = "N_glands_seurat")
mmDE <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list mixedmodelDE.xlsx", sheet = "N_glands")

N_glands <- dplyr::full_join(seurat, mmDE, by = "Gene")
N_glands <- distinct(N_glands, Gene, .keep_all = T)
write.csv(N_glands, "C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/FULL_N_glands.csv")

rm(seurat, mmDE)
seurat <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list seurat.xlsx", sheet = "PV_mucosa_seurat")
mmDE <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list mixedmodelDE.xlsx", sheet = "PV_mucosa")

PV_mucosa <- dplyr::full_join(seurat, mmDE, by = "Gene")
PV_mucosa <- distinct(PV_mucosa, Gene, .keep_all = T)
write.csv(PV_mucosa, "C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/FULL_PV_mucosa.csv")

rm(seurat, mmDE)
seurat <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list seurat.xlsx", sheet = "N_mucosa_seurat")
mmDE <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list mixedmodelDE.xlsx", sheet = "N_mucosa")

N_mucosa <- dplyr::full_join(seurat, mmDE, by = "Gene")
N_mucosa <- distinct(N_mucosa, Gene, .keep_all = T)
write.csv(N_mucosa, "C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/FULL_N_mucosa.csv")

rm(seurat, mmDE)
seurat <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list seurat.xlsx", sheet = "PV_stroma_seurat")
mmDE <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list mixedmodelDE.xlsx", sheet = "PV_stroma")

PV_stroma <- dplyr::full_join(seurat, mmDE, by = "Gene")
PV_stroma <- distinct(PV_stroma, Gene, .keep_all = T)
write.csv(PV_stroma, "C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/FULL_PV_stroma.csv")


rm(seurat, mmDE)
seurat <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list seurat.xlsx", sheet = "N_stroma_seurat")
mmDE <- read_excel("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/Cheng DEG list mixedmodelDE.xlsx", sheet = "N_stroma")

N_stroma <- dplyr::full_join(seurat, mmDE, by = "Gene")
N_stroma <- distinct(N_stroma, Gene, .keep_all = T)
write.csv(N_stroma, "C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/FULL_N_stroma.csv")



gene_list <- list(PV_glands = PV_glands$Gene,
                  N_glands = seurat$Gene,
                  PV_mucosa = PV_mucosa$Gene,
                  N_mucosa = N_mucosa$Gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red", "green", "yellow"),
                            cat.col = c("blue", "red", "green", "yellow"),
                            cex = 2,lty = "blank",
                            cat.cex = 2)

cowplot::plot_grid(VennDiagram)


gene_list <- list(PV_stroma = PV_stroma$Gene,
                  N_stroma = N_stroma$Gene)
VennDiagram <- venn.diagram(x = gene_list, filename = NULL,
                            fill = c("blue", "red"),
                            cat.col = c("blue", "red"),
                            cex = 2,lty = "blank",
                            cat.cex = 2)

cowplot::plot_grid(VennDiagram)



