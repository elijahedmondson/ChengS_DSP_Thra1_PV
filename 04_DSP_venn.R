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

universe <- read.csv("C:/Users/edmondsonef/Desktop/ChengS_DSP Results/DEG/DEG_universe.csv")
universe <- dplyr::select(universe,SYMBOL,ENTREZID)
head(universe)


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



