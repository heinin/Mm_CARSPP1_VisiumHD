#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 11/30/2023
# Description: Plotting mouse Visium HD data
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(googlesheets4)
library(gprofiler2)
library(scCustomize)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/CART_plot_functions.R")
source("/home/hnatri/13384_CART/13384_Tumors/SPP1_ms/HGG_SPP1/13384_tumor_ms_themes.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Import data
#==============================================================================#

# The final object
seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/Mm_VisiumHD/tumor_obj_w_all_tacco_labels_09_09_2026.rds")

head(seurat_data@meta.data)
table(seurat_data$Treatment)

# Fixing TMA3 treatment group
seurat_data$Treatment <- ifelse(seurat_data$TMA == 3 & seurat_data$Sample %in% c(166, 159), "Anti-Spp1", seurat_data$Treatment)
seurat_data$Treatment <- gsub(" ", "", seurat_data$Treatment)
table(seurat_data$Treatment, seurat_data$TMA)

seurat_data$Treatment <- gsub("Anti-Spp1/CAR", "Anti-Spp1+CAR", seurat_data$Treatment)

seurat_data$Treatment_Day <- paste0(seurat_data$Treatment, "_", seurat_data$Day)

treatment_cols <- c("Anti-Spp1+CAR" = "aquamarine3",
                    "CAR" = "cyan4",
                    "Anti-Spp1" = "mediumpurple",
                    "Tumor" = "deeppink3")

day_cols <- c("16" = "gray40",
              "20" = "azure3")

DimPlot(seurat_data,
        reduction = "sp",
        group.by = "projected_snn_res.0.8")

p1 <- DimPlot(seurat_data,
        reduction = "sp",
        group.by = "Treatment",
        cols = treatment_cols) +
  theme_classic() +
  coord_fixed() +
  NoAxes()

p2 <- DimPlot(seurat_data,
        reduction = "sp",
        group.by = "Day",
        cols = day_cols) +
  theme_classic() +
  coord_fixed() +
  NoAxes()

p1 / p2

filename <- "/home/hnatri/VisiumHD_dimplot_treatment.pdf"
pdf(file = filename,
    width = 20,
    height = 8)

p1

dev.off()

filename <- "/home/hnatri/VisiumHD_dimplot_day.pdf"
pdf(file = filename,
    width = 20,
    height = 8)

p2

dev.off()


FeaturePlot(seurat_data,
            reduction = "sp",
            features = c("nCount_Spatial", "nFeature_Spatial"),
            split.by = "TMA") &
  theme_classic() &
  coord_fixed()

filename <- "/home/hnatri/VisiumHD_counts_features.pdf"
pdf(file = filename,
    width = 20,
    height = 12)

FeaturePlot(seurat_data,
            reduction = "sp",
            features = c("nCount_Spatial", "nFeature_Spatial"),
            ncol = 1,
            raster = T) &
  theme_classic() &
  coord_fixed()

dev.off()

FeaturePlot_scCustom(seurat_object = seurat_data,
                     features =  c("nCount_Spatial", "nFeature_Spatial"),
                     reduction = "sp",
                     num_columns = 1,
                     #raster = F,
                     #raster.dpi = c(1024, 1024),
                     #slot = "data",
                     #pt.size = 0.1,
                     #na_cutoff = 1,
                     order = T)

# GBM genes
gbm_genes <- c("MGMT", "IDH1", "IDH2", "EGFR", "PTEN", "TP53", "TERT")
gbm_genes <- gorth(gbm_genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

FeaturePlot(seurat_data,
            reduction = "sp",
            features = gbm_genes,
            ncol = 1,
            raster = T,
            raster.dpi = c(1024, 1024)) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
  theme_classic() &
  coord_fixed()

# "Mgmt"  "Idh1"  "Idh2"  "Egfr"  "Trp53" "Tert" 
FeaturePlot_scCustom(seurat_object = seurat_data,
                     features = "C1qa",
                     reduction = "sp",
                     num_columns = 1,
                     raster = F,
                     #raster.dpi = c(1024, 1024),
                     #slot = "data",
                     #pt.size = 0.1,
                     #na_cutoff = 20,
                     order = T)

filename <- "/home/hnatri/VisiumHD_GBM_genes.png"
ggsave(file = filename,
    width = 12,
    height = 24)
    #dpi = 300)
#filename <- "/home/hnatri/VisiumHD_GBM_genes.pdf"
#pdf(file = filename,
#    width = 12,
#    height = 24)


FeaturePlot_scCustom(seurat_object = seurat_data,
                     features = c(gbm_genes, c("C1qa")),
                     reduction = "sp",
                     num_columns = 1,
                     raster = F,
                     #raster.dpi = c(1024, 1024),
                     #slot = "data",
                     #pt.size = 1,
                     #na_cutoff = 1,
                     order = T)

dev.off()


filename <- "/home/hnatri/VisiumHD_TME_genes.png"
ggsave(file = filename,
       width = 10,
       height = 26)
#filename <- "/home/hnatri/VisiumHD_TME_genes.pdf"
#pdf(file = filename,
#    width = 8,
#    height = 12)
FeaturePlot_scCustom(seurat_object = seurat_data,
                     features = c("Tigit", "Cd8a", "Cd3e", "Spp1", "Tnfrsf12a", "Tnfsf12", "C1qa", "Cd44", "Col1a1", "Pten"),
                     reduction = "sp",
                     num_columns = 1,
                     #raster = F,
                     #raster.dpi = c(1024, 1024),
                     slot = "data",
                     pt.size = 1,
                     na_cutoff = 1,
                     order = T) &
  coord_fixed() &
  NoAxes()

dev.off()

FeaturePlot_scCustom(seurat_object = seurat_data,
                     features = c("Spp1", "C1qa", "Tnfrsf12a", "Tnfsf12"),
                     reduction = "sp",
                     num_columns = 1,
                     #raster = F,
                     #raster.dpi = c(1024, 1024),
                     slot = "data",
                     pt.size = 1,
                     na_cutoff = 1,
                     order = T) &
  coord_fixed() &
  NoAxes()

grep("Ifn", rownames(seurat_data), value = T)

FeaturePlot_scCustom(seurat_object = seurat_data,
                     features = c("Spp1", "Jak1", "Jak2", "Ifna", "Ifnb1"),
                     reduction = "sp",
                     num_columns = 1,
                     #raster = F,
                     #raster.dpi = c(1024, 1024),
                     slot = "data",
                     pt.size = 1,
                     na_cutoff = 1,
                     order = T) &
  coord_fixed() &
  NoAxes()

VlnPlot(seurat_data,
        features = c("Spp1", "C1qa", "Tnfrsf12a", "Tnfsf12"),
        split.by = "Treatment",
        group.by = "Day",
        pt.size = 0,
        ncol = 1,
        cols = treatment_cols) &
  theme_classic()

DotPlot(seurat_data,
        features = c(c("Spp1", "Jak1", "Jak2"), grep("Ifn", rownames(seurat_data), value = T)),
        group.by = "predicted.id_class_immune") &
  theme_classic() &
  coord_flip() &
  RotatedAxis() &
  xlab("") &
  ylab("")

DotPlot(seurat_data,
        features = c("Spp1", "Jak1", "Jak2"),
        #group.by = "predicted.id_class_immune",
        group.by = "Treatment_Day") &
  theme_classic() &
  coord_flip() &
  RotatedAxis() &
  xlab("") &
  ylab("")

VlnPlot(seurat_data,
        features = c(c("Spp1", "Jak1", "Jak2"), grep("Ifn", rownames(seurat_data), value = T))[1:14],
        group.by = "Treatment",
        split.by = "Day",
        pt.size = 0,
        ncol = 2,
        cols = treatment_cols) &
  theme_classic() &
  xlab("")

VlnPlot(seurat_data,
        features = c("Spp1", "C1qa", "Tnfrsf12a", "Tnfsf12"),
        group.by = "predicted.id_class_immune",
        pt.size = 0,
        ncol = 1) &
  theme_classic() &
  RotatedAxis() &
  NoLegend() &
  xlab("")

seurat_data$Treatment_Day <- paste0(seurat_data$Treatment, "_", seurat_data$Day)
tams <- subset(seurat_data, subset = predicted.id_class_immune %in% c("TAM 1", "TAM 2", "prol. TAM"))

DotPlot(tams,
        features = c("Spp1", "C1qa", "Tnfrsf12a", "Tnfsf12"),
        group.by = "Treatment_Day") &
  theme_classic()

VlnPlot(tams,
        features = c("Spp1", "C1qa", "Tnfrsf12a", "Tnfsf12"),
        group.by = "Treatment_Day",
        split.by = "predicted.id_class_immune",
        pt.size = 0,
        ncol = 1,
        cols = treatment_cols) &
  theme_classic() &
  RotatedAxis() &
  xlab("")

FeaturePlot(subset(seurat_data, subset = Sample %in% c(172)),
            features = c("Tnfrsf12a", "Tnfsf12"),
            #split.by = "Sample",
            blend = T) &
  theme_classic() &
  RotatedAxis() &
  xlab("")

FeaturePlot(tams,
        features = c("Tnfrsf12a", "Tnfsf12"),
        split.by = "Sample",
        blend = T,
        reduction = "umap.sketch") &
  theme_classic() &
  RotatedAxis() &
  xlab("")

VlnPlot(seurat_data,
        features = c("Tigit", "Cd8a", "Cd3e", "Spp1", "Tnfrsf12a", "Tnfsf12", "C1qa", "Cd44", "Col1a1", "Pten"),
        split.by = "Treatment",
        group.by = "Day",
        pt.size = 0,
        ncol = 1,
        cols = treatment_cols) &
  theme_classic()

seurat_data$Treatment_Day <- paste0(seurat_data$Treatment, "_", seurat_data$Day)
seurat_data$Treatment_Day <- factor(seurat_data$Treatment_Day,
                                    levels = c("Tumor_16", "Tumor_20",
                                               "Anti-Spp1_16", "Anti-Spp1_20",
                                               "CAR_16", "CAR_20",
                                               "Anti-Spp1+CAR_16", "Anti-Spp1+CAR_20"))

dotp <- DotPlot(seurat_data,
        features = c("Tigit", "Cd8a", "Cd3e", "Spp1", "C1qa", "Cd44", "Col1a1", "Pten"),
        group.by = "Treatment_Day",
        cols = c("gray89", "tomato3")) +
  coord_flip() +
  theme_classic() +
  RotatedAxis() +
  xlab("") +
  ylab("")

dotp

filename <- "/home/hnatri/VisiumHD_dotplot.pdf"
pdf(file = filename,
    width = 5,
    height = 3)

dotp

dev.off()

#==============================================================================#
# CAR T phenotypes
#==============================================================================#

#counts_data <- LayerData(seurat_data, layer = "counts")
#
#counts_data[1:3,1:3]
#tail(rownames(counts_data))
#
#
#car_bins <- counts_data["mIL13",]
#
#identical(colnames(seurat_data), names(car_bins))
#seurat_data$CAR_counts <- car_bins
#
#seurat_data$CARpos <- ifelse(seurat_data$CAR_counts > 0, "CAR", "nonCAR")
#
#table(seurat_data$CARpos, seurat_data$Sample)
#
#carpos <- names(car_bins[car_bins > 0])
#carpos_list <- list("CARpos" = carpos)
#
#pps <- DimPlot_scCustom(seurat_object = seurat_data,
#                 cells.highlight = carpos_list,
#                 cols.highlight = "tomato3",
#                 colors_use = "gray90",
#                 sizes.highlight = 2,
#                 pt.size = 1,
#                 reduction = "sp",
#                 order = T,
#                 raster = F) &
#  coord_fixed() &
#  NoAxes() &
#  NoLegend()
#
#filename <- "/home/hnatri/VisiumHD_CARpos.png"
#ggsave(file = filename,
#       width = 10,
#       height = 4)
#
#pps
#
#dev.off()

# Mouse Immune Cell Marker Summary (Mus musculus)
#
# - Total Leukocytes: Ptprc (CD45)
# - T Cells (Pan): Cd3e, Cd3g, Cd247 (CD3)
# - Helper T Cells: Cd4
# - Cytotoxic T Cells: Cd8a, Cd8b1
# - Regulatory T Cells (Tregs): Foxp3, Cd25 (Il2ra), Ctla4
# - B Cells: Cd19, Ms4a1 (CD20), Cd79a
# - NK Cells: Ncr1 (NKp46), Klrb1c (NK1.1), Cd226
# - Macrophages (Pan): Adgre1 (F4/80), Cd68, Cd11b (Itgam)
# - M1 Macrophages (Pro-inflammatory): Nos2 (iNOS), Cd80, Cd86
# - M2 Macrophages (Anti-inflammatory): Cd163, Mrc1 (CD206), Arg1
# - Conventional Dendritic Cells (cDCs): Itgax (CD11c), H2-Ab1 (MHCII)
# - Plasmacytoid DCs (pDCs): Siglech, Bst2 (CD317)
# - Neutrophils: Ly6g, Ly6c, S100a8
# - Monocytes: Ly6c, Cd11b (Itgam)
# - Hematopoietic Stem/Progenitor Cells: Cd34, Kit (c-Kit)
# 
# Key Activation/Memory Markers
# 
# - Naïve T Cells: Sell (CD62L), Ccr7
# - Memory T Cells: Cd44, Il7r (CD127)
# - Exhaustion Markers: Pdcd1 (PD-1), Cd274 (PD-L1), Havcr2 (TIM-3)

mouse_immune <- c("Ptprc", "Cd3e", "Cd3g", "Cd247", "Cd4", "Cd8a", "Cd8b1", 
                  "Foxp3", "Cd25", "Il2ra", "Ctla4", "Cd19", "Ms4a1", "Cd79a",
                  "Ncr1", "NKp46", "Klrb1c", "Cd226",
                  "Adgre1", "Cd68", "Cd11b", "Itgam",
                  "Nos2", "Cd80", "Cd86",
                  "Cd163", "Mrc1", "Arg1", "Itgax", "Cd11c",
                  "Siglech", "Bst2", "Ly6g", "Ly6c", "S100a8",
                  "Ly6c", "Cd34", "Kit")

lymph_markers <- c("Ptprc", "Cd3e", "Cd3g", "Cd247", "Cd4", "Cd8a", "Cd8b1", 
                   "Foxp3", "Cd25", "Il2ra", "Ctla4", "Cd19", "Ms4a1", "Cd79a",
                   "Ncr1", "NKp46", "Klrb1c", "Cd226")

#==============================================================================#
# Cell type annotations
#==============================================================================#

seurat_data$TACCO_myeloid_label <- ifelse(is.na(seurat_data$TACCO_myeloid_label), "NA", seurat_data$TACCO_myeloid_label)
seurat_data$projected_snn_res.1.2 <- ifelse(is.na(seurat_data$projected_snn_res.1.2), "NA", seurat_data$projected_snn_res.1.2)
seurat_data$predicted.id_tam <- ifelse(is.na(seurat_data$predicted.id_tam), "NA", seurat_data$predicted.id_tam)

DimPlot(seurat_data,
        group.by = "TACCO_myeloid_label",
        reduction = "sp") &
  coord_fixed()

projected_snn_res.1.2_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(unique(seurat_data$projected_snn_res.1.2)))
names(projected_snn_res.1.2_col) <- sort(unique(seurat_data$projected_snn_res.1.2))

TACCO_myeloid_label_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(unique(seurat_data$TACCO_myeloid_label)))
names(TACCO_myeloid_label_col) <- sort(unique(seurat_data$TACCO_myeloid_label))

predicted.id_tam_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(unique(seurat_data$predicted.id_tam)))
names(predicted.id_tam_col) <- sort(unique(seurat_data$predicted.id_tam))

create_barplot(seurat_data,
               group_var = "Treatment_Day",
               plot_var = "TACCO_myeloid_label",
               plot_levels = sort(unique(seurat_data$TACCO_myeloid_label)),
               group_levels = sort(unique(seurat_data$Treatment_Day)),
               plot_colors = TACCO_myeloid_label_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Cell type")

prop_df <- seurat_data@meta.data %>%
  dplyr::count(Sample, Treatment_Day, TACCO_myeloid_label, name = "n") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  # fill in zeros for celltypes absent in a sample
  tidyr::complete(tidyr::nesting(Sample, Treatment_Day), TACCO_myeloid_label,
                  fill = list(n = 0, prop = 0))

ggplot(prop_df, aes(x = Treatment_Day, y = prop, fill = Treatment_Day)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 0.9, alpha = 0.8) +
  facet_wrap(~ TACCO_myeloid_label, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

day16 <- subset(seurat_data, subset = Day == 16)

prop_df <- day16@meta.data %>%
  dplyr::count(Treatment_Day, TACCO_myeloid_label, name = "n") %>%
  dplyr::group_by(Treatment_Day) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

ggplot(prop_df, aes(x = TACCO_myeloid_label, y = prop, fill = Treatment_Day)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

day20 <- subset(seurat_data, subset = Day == 20)

prop_df <- day20@meta.data %>%
  dplyr::count(Treatment_Day, TACCO_myeloid_label, name = "n") %>%
  dplyr::group_by(Treatment_Day) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

ggplot(prop_df, aes(x = TACCO_myeloid_label, y = prop, fill = Treatment_Day)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Cluster
create_barplot(seurat_data,
               group_var = "Treatment_Day",
               plot_var = "projected_snn_res.1.2",
               plot_levels = sort(unique(seurat_data$projected_snn_res.1.2)),
               group_levels = sort(unique(seurat_data$Treatment_Day)),
               plot_colors = projected_snn_res.1.2_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Cell type")

prop_df <- seurat_data@meta.data %>%
  dplyr::count(Sample, Treatment_Day, projected_snn_res.1.2, name = "n") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  # fill in zeros for celltypes absent in a sample
  tidyr::complete(tidyr::nesting(Sample, Treatment_Day), projected_snn_res.1.2,
                  fill = list(n = 0, prop = 0))

ggplot(prop_df, aes(x = Treatment_Day, y = prop, fill = Treatment_Day)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 0.9, alpha = 0.8) +
  facet_wrap(~ projected_snn_res.1.2, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

day16 <- subset(seurat_data, subset = Day == 16)

prop_df <- day16@meta.data %>%
  dplyr::count(Treatment_Day, projected_snn_res.1.2, name = "n") %>%
  dplyr::group_by(Treatment_Day) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

ggplot(prop_df, aes(x = as.character(projected_snn_res.1.2), y = prop, fill = Treatment_Day)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

day20 <- subset(seurat_data, subset = Day == 20)

prop_df <- day20@meta.data %>%
  dplyr::count(Treatment_Day, projected_snn_res.1.2, name = "n") %>%
  dplyr::group_by(Treatment_Day) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

ggplot(prop_df, aes(x = as.character(projected_snn_res.1.2), y = prop, fill = Treatment_Day)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


create_barplot(seurat_data,
               group_var = "Treatment_Day",
               plot_var = "projected_snn_res.1.2",
               plot_levels = sort(unique(seurat_data$projected_snn_res.1.2)),
               group_levels = sort(unique(seurat_data$Treatment_Day)),
               plot_colors = projected_snn_res.1.2_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Cell type")


#==============================================================================#
# Annotation concordance
#==============================================================================#

DimPlot(seurat_data,
        reduction = "sp",
        group.by = "TACCO_myeloid_label")

DimPlot(seurat_data,
        reduction = "sp",
        group.by = "projected_snn_res.0.8")

table(seurat_data$projected_snn_res.0.8,
      seurat_data$TACCO_myeloid_label)

table(seurat_data$projected_snn_res.0.8,
      seurat_data$predicted.id_tam)

table(seurat_data$predicted.id_tam,
      seurat_data$TACCO_myeloid_label)

table(seurat_data$predicted.id_tam) %>%
  as.data.frame()

p1 <- create_barplot(seurat_data,
               group_var = "TACCO_myeloid_label",
               plot_var = "projected_snn_res.1.2",
               plot_levels = sort((unique(seurat_data$projected_snn_res.1.2))),
               group_levels = sort((unique(seurat_data$TACCO_myeloid_label))),
               plot_colors = projected_snn_res.1.2_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Seurat cluster")

p2 <- create_barplot(seurat_data,
               group_var = "projected_snn_res.1.2",
               plot_var = "TACCO_myeloid_label",
               plot_levels = sort((unique(seurat_data$TACCO_myeloid_label))),
               group_levels = sort((unique(seurat_data$projected_snn_res.1.2))),
               plot_colors = TACCO_myeloid_label_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "TACCO")

p1 + p2

p3 <- create_barplot(seurat_data,
               group_var = "predicted.id_tam",
               plot_var = "projected_snn_res.1.2",
               plot_levels = sort((unique(seurat_data$projected_snn_res.1.2))),
               group_levels = sort((unique(seurat_data$predicted.id_tam))),
               plot_colors = projected_snn_res.1.2_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Seurat cluster")

p4 <- create_barplot(seurat_data,
               group_var = "projected_snn_res.1.2",
               plot_var = "predicted.id_tam",
               plot_levels = sort((unique(seurat_data$predicted.id_tam))),
               group_levels = sort((unique(seurat_data$projected_snn_res.1.2))),
               plot_colors = predicted.id_tam_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Seurat label transfer")

p3 + p4

p5 <- create_barplot(seurat_data,
               group_var = "TACCO_myeloid_label",
               plot_var = "predicted.id_tam",
               plot_levels = sort((unique(seurat_data$predicted.id_tam))),
               group_levels = sort((unique(seurat_data$TACCO_myeloid_label))),
               plot_colors = predicted.id_tam_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "Seurat label transfer")

p6 <- create_barplot(seurat_data,
               group_var = "predicted.id_tam",
               plot_var = "TACCO_myeloid_label",
               plot_levels = sort((unique(seurat_data$TACCO_myeloid_label))),
               group_levels = sort((unique(seurat_data$predicted.id_tam))),
               plot_colors = TACCO_myeloid_label_col,
               var_names =  c("Frequency (%)", ""),
               legend_title = "TACCO")

p5 + p6

#==============================================================================#
# Top marker expression
#==============================================================================#

Idents(seurat_data) <- seurat_data$TACCO_myeloid_label
tacco_markers <- FindAllMarkers(seurat_data,
                                only.pos = T)
tacco_markers_sig <- tacco_markers %>%
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC)

DotPlot(seurat_data,
        features = unique(tacco_markers_sig$gene),
        group.by = "TACCO_myeloid_label") &
  theme_classic() &
  #coord_flip() &
  RotatedAxis() &
  xlab("") &
  ylab("")

Idents(seurat_data) <- seurat_data$predicted.id_tam
id_markers <- FindAllMarkers(seurat_data,
                                only.pos = T)
id_markers_sig <- id_markers %>%
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC)

DotPlot(seurat_data,
        features = unique(id_markers_sig$gene),
        group.by = "predicted.id_tam") &
  theme_classic() &
  #coord_flip() &
  RotatedAxis() &
  xlab("") &
  ylab("")

Idents(seurat_data) <- seurat_data$projected_snn_res.1.2
cluster_markers <- FindAllMarkers(seurat_data,
                                only.pos = T)
cluster_markers_sig <- cluster_markers %>%
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

DotPlot(seurat_data,
        features = unique(cluster_markers_sig$gene),
        group.by = "projected_snn_res.1.2") &
  theme_classic() &
  #coord_flip() &
  RotatedAxis() &
  xlab("") &
  ylab("")

DotPlot(seurat_data,
        features = c("Olig1", "Olig2", "Idh1", "Egfr", "Pten", "Trp53", "Ptprc", "Col1a1", "Col1a2", "Vim"),
        group.by = "projected_snn_res.1.2") &
  theme_classic() &
  #coord_flip() &
  RotatedAxis() &
  xlab("") &
  ylab("")

DimPlot(seurat_data,
        reduction = "sp",
        group.by = "projected_snn_res.1.2") +
  theme_classic() +
  coord_fixed() +
  NoAxes()

# Cluster 0
mk0 <- FindMarkers(seurat_data, ident.1 = "0", only.pos = TRUE,
                   min.pct = 0.25, logfc.threshold = 0.25)
mk0 <- mk0 %>% arrange(plyr::desc(avg_log2FC))
head(mk0, 30)          # <- the honest answer to "what is cluster 0"

mg_markers <- c("P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Siglech",
                "Sall1", "Selplg", "Csf1r", "Aif1", "Fcrls", "C1qa", "C1qb")
mg_markers <- mg_markers[mg_markers %in% rownames(seurat_data)]  # drop absent genes

seurat_data$projected_snn_res.1.2 <- factor(seurat_data$projected_snn_res.1.2, levels = sort(unique(seurat_data$projected_snn_res.1.2)))

DotPlot(seurat_data,
        features = mg_markers,
        group.by = "projected_snn_res.1.2") + coord_flip() +
  ggplot2::ggtitle("Microglia")

DotPlot(seurat_data, features = mg_markers)$data %>%
  filter(id == "0") %>%
  select(gene = features.plot, avg.exp.scaled, pct.exp) %>%
  arrange(desc(pct.exp))

qc_cols <- intersect(c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
                     colnames(seurat_data@meta.data))

seurat_data@meta.data %>%
  group_by(cluster = Idents(seurat_data)) %>%
  dplyr::summarise(across(all_of(qc_cols), median), n = dplyr::n()) %>%
  arrange(cluster)

VlnPlot(seurat_data, features = qc_cols, pt.size = 0) &
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7))

#==============================================================================#
# Myeloid expression
#==============================================================================#

DotPlot(seurat_data,
        features = unique(c(lymph_markers, mouse_immune, c("Spp1", "Jak1", "Jak2"), grep("Ifn", rownames(seurat_data), value = T))),
        group.by = "projected_snn_res.1.2") &
  theme_classic() &
  coord_flip() &
  RotatedAxis() &
  xlab("") &
  ylab("")

immune_only <- subset(seurat_data, subset = projected_snn_res.1.2 %in% c(0, 1, 3, 5, 12, 13))

DotPlot(seurat_data,
        features = unique(c(lymph_markers, mouse_immune, c("Spp1", "Jak1", "Jak2"))),
        group.by = "Treatment_Day") &
  theme_classic() &
  RotatedAxis() &
  xlab("") &
  ylab("")

VlnPlot(immune_only,
        features = unique(c(lymph_markers, mouse_immune, c("Spp1", "Jak1", "Jak2")))[1:20],
        group.by = "Treatment_Day",
        ncol = 8,
        pt.size = 0) &
  theme_classic() &
  NoLegend() &
  RotatedAxis() &
  xlab("")

VlnPlot(immune_only,
        features = unique(c(lymph_markers, mouse_immune, c("Spp1", "Jak1", "Jak2")))[21:40],
        group.by = "Treatment_Day",
        ncol = 8,
        pt.size = 0) &
  theme_classic() &
  NoLegend() &
  RotatedAxis() &
  xlab("")

# Top markers for each group
Idents(immune_only) <- immune_only$Treatment_Day
treatment_day_immune_markers <- FindAllMarkers(immune_only,
                                               only.pos = T)

treatment_day_immune_markers_top <- treatment_day_immune_markers %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

DotPlot(immune_only,
        features = unique(treatment_day_immune_markers_top$gene),
        group.by = "Treatment_Day") &
  theme_classic() &
  RotatedAxis() &
  xlab("") &
  ylab("")


#==============================================================================#
# CAF expression
#==============================================================================#

#==============================================================================#
# Lymphoid expression
#==============================================================================#




