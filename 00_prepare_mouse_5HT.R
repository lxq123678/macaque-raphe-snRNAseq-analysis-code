library(Seurat)
library(harmony)
mouse_Allen <- readRDS(file = "./data/Allen_sero.rds")
mouse_sero <- readRDS(file = "./data/Pet1_sero.rds")
#harmony去批次

mouse_sero$group <- "mouse1"
mouse_Allen$group <- "mouse2"

counts1 <- GetAssayData(mouse_sero, assay = "RNA", slot = "counts")
mouse1 <- CreateSeuratObject(counts = counts1,meta.data = mouse_sero@meta.data) 


counts2 <- GetAssayData(mouse_Allen, assay = "RNA", slot = "counts")
mouse2 <- CreateSeuratObject(counts = counts2,meta.data = mouse_Allen@meta.data)

merged <- merge(mouse1,mouse2)
table(merged$group)
options(future.globals.maxSize= 891289600 * 160)


merged = SCTransform(
  merged, assay="RNA",
  ncells=ncol(merged[["RNA"]]),
  variable.features.n=5000,
  #return.only.var.genes = T,
  #vars.to.regress="percent.mt",
  method="glmGamPoi"
)

merged <- RunPCA(merged,features = VariableFeatures(object = merged), npcs = 50)
merged = RunHarmony(merged,
                    group.by.vars = "group",
                    assay.use = "SCT",
                    dims.use = 1:10,
                    plot_convergence = TRUE
)
merged <- FindNeighbors(merged,reduction = "harmony",dims = 1:10)
library(clustree)

merged <- FindClusters(merged, resolution = seq(0.1, 1, by = 0.1))

clustree(merged, prefix = "SCT_snn_res.")

merged <- FindClusters(merged,resolution = 0.5)
merged <- RunUMAP(merged,reduction="harmony", dims = 1:10)
DimPlot(merged,,label = T,group.by = "seurat_clusters")
DimPlot(merged,,label = T,group.by = "group")
sce <- PrepSCTFindMarkers(merged, assay = "SCT", verbose = TRUE)

all_markers <- FindAllMarkers(
  sce, test.use = "wilcox", only.pos = F,logfc.threshold = 0.25,min.pct = 0.25
)
sig_genes <- subset(all_markers, p_val_adj < 0.05 & avg_log2FC > 1)
top5 <- sig_genes %>%
  subset(avg_log2FC > 1) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
DotPlot(sce, features =unique(top5$gene),
        group.by = "seurat_clusters",
        cols = c("lightgrey", "#E64B35FF"),  # 颜色梯度：灰→红
        dot.scale = 6,                       # 点大小缩放因子
        scale = TRUE) +RotatedAxis()+
  theme(
    axis.text.x = element_text(size = 10),  # X轴标签垂直
  )
merged@meta.data$renew_type <- paste0("5-HT",'_',as.numeric(as.character(merged@meta.data$seurat_clusters)) + 1)

mouse_sero.combined <- JoinLayers(merged,assay = 'RNA')
saveRDS(mouse_sero.combined,file = "./data/mouse_sero.combined_res0.5.rds")


