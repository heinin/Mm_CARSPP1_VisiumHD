library(Seurat)

dataa <- readRDS("/scratch/rmudunuri/cart_project/13384_tumor/MouseVisiumHD/objects/full_obj_no_CAR/merged_obj_unfiltered_16um_bins.rds")

dataa <- JoinLayers(dataa)

head(dataa@meta.data)
unique(dataa$Treatment)


dataa@meta.data %>%
  select(c("Sample", "Treatment", "Day", "TMA")) %>%
  distinct()

dataas <- FetchData(dataa,
                    vars = "mIL13",
                    layer = "counts")

dataa$CAR_counts <- dataas
dataa$isCAR <- ifelse(dataa$CAR_counts>=1, "CARpos", "nonCAR")
dataa <- subset(dataa, subset = Treatment %in% c("Anti-Spp1/CAR", "CAR"))
dataa <- NormalizeData(dataa)

mouse_immune <- unique(c("Ptprc", "Cd3d", "Cd3e", "Cd3g", "Cd4a", "Cd247", "Cd4", "Cd8a", "Cd8b1", "Tigit",
                  "Foxp3", "Cd25", "Cd27", "Cd28", "Il2ra", "Ctla4", "Cd19", "Ms4a1", "Cd79a",
                  "Ncr1", "NKp46", "Klrb1c", "Cd226", "Cd127",
                  "Adgre1", "Cd68", "Cd11b", "Itgam",
                  "Nos2", "Cd80", "Cd86",
                  "Cd163", "Mrc1", "Arg1", "Itgax", "Cd11c",
                  "Siglech", "Bst2", "Ly6g", "Ly6c", "S100a8",
                  "Ly6c", "Cd34", "Kit", "Sell", "Ccr7", "Cd44",
                  "Pdcd1", "Gitr", "Cd45ra"))


VlnPlot(subset(dataa, subset = Day == 16),
        features = mouse_immune,
        group.by = "Treatment",
        split.by = "isCAR",
        pt.size = 0,
        ncol = 6) &
  xlab("") &
  theme_classic(base_size = 8)

DotPlot(subset(dataa, subset = Day == 16),
        features = c("Cd3d", "Cd8a", "Cd4a", "Tigit", "Foxp3", "Cd25", "Cd27", "Cd28", "Cd127", "Ctla4", "Ccr7", "Gitr", "Cd45ra", "Cd44"),
        split.by = "isCAR",
        group.by = "Treatment") +
  coord_flip() +
  RotatedAxis() +
  ggtitle("Day 16")

VlnPlot(subset(dataa, subset = Day == 16),
        features = c("Cd3d", "Cd8a", "Cd4a", "Tigit", "Foxp3", "Cd25", "Cd27", "Cd28", "Cd127", "Ctla4", "Ccr7", "Gitr", "Cd45ra", "Cd44"),
        split.by = "isCAR",
        group.by = "Treatment") +
  coord_flip() +
  RotatedAxis() +
  ggtitle("Day 16")

VlnPlot(subset(dataa, subset = Day == 16 & isCAR == "CARpos"),
        features = mouse_immune[1:25],
        ncol = 7,
        #split.by = "isCAR",
        group.by = "Treatment",
        pt.size = 0) &
  theme_classic() &
  theme(base.size = 5) &
  NoLegend() &
  RotatedAxis() &
  xlab("")

VlnPlot(subset(dataa, subset = Day == 16 & isCAR == "CARpos"),
        features = mouse_immune[26:49],
        ncol = 6,
        #split.by = "isCAR",
        group.by = "Treatment",
        pt.size = 0) &
  theme_classic() &
  theme(base.size = 5) &
  NoLegend() &
  RotatedAxis() &
  xlab("")
  #ggtitle("Day 16")

VlnPlot(subset(dataa, subset = Day == 16),
        features = c("Cd3d", "Cd8a"),
        split.by = "isCAR",
        group.by = "Treatment",
        pt.size = 0)


dataa2 <- subset(dataa, subset = Day == 16)
table(dataa2$Treatment, dataa2$isCAR)

dataa3 <- subset(dataa2, subset = isCAR == "CARpos")

Idents(dataa3) <- dataa3$Treatment
des <- FindMarkers(dataa3,
                  ident.1 = "CAR",
                  ident.2 = "Anti-Spp1/CAR")

min(des$p_val)

tops <- des %>% filter(p_val < 0.01) %>% arrange(-(abs(avg_log2FC))) %>% head(n = 20)

VlnPlot(subset(dataa, subset = Day == 16 & isCAR == "CARpos"),
        features = rownames(tops),
        ncol = 6,
        #split.by = "isCAR",
        group.by = "Treatment",
        pt.size = 0) &
  theme_classic() &
  theme(base.size = 5) &
  NoLegend() &
  RotatedAxis() &
  xlab("")

des %>% filter(p_val < 0.01) %>% arrange(-(abs(avg_log2FC))) %>% rownames()
