# ---- 1. Load required packages and read input Seurat objects
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
# Read two mouse serotonin-related datasets from RDS files
mouse_Allen <- readRDS(file = "./data/Allen_sero.rds")
mouse_sero <- readRDS(file = "./data/Pet1_sero.rds")

# ---- 2. Add sample/group labels for batch correction 
# Assign a group label to each dataset; this will be used as the batch variable in Harmony 
mouse_sero$group <- "mouse1"
mouse_Allen$group <- "mouse2"

# ---- 3. Rebuild Seurat objects from count matrices while keeping metadata
# Extract raw count matrices from the RNA assay
counts1 <- GetAssayData(mouse_sero, assay = "RNA", slot = "counts")
mouse1 <- CreateSeuratObject(counts = counts1,meta.data = mouse_sero@meta.data) 
# Create new Seurat objects using the extracted counts and the original metadata
counts2 <- GetAssayData(mouse_Allen, assay = "RNA", slot = "counts")
mouse2 <- CreateSeuratObject(counts = counts2,meta.data = mouse_Allen@meta.data)

# ---- 4. Merge the two datasets and check group composition
merged <- merge(mouse1,mouse2)
# Check how many cells are present in each group after merging
table(merged$group)

# ---- 5. Increase the future global memory limit for parallel/large jobs
options(future.globals.maxSize= 891289600 * 160)

# ---- 6. Normalize and stabilize variance using SCTransform
# Use all cells for SCTransform, select 5,000 variable features, and use glmGamPoi as the modeling method 
merged = SCTransform(
  merged, assay="RNA",
  ncells=ncol(merged[["RNA"]]),
  variable.features.n=5000,
  method="glmGamPoi"
)

# ---- 7. Run PCA on the SCTransform variable features
merged <- RunPCA(merged,features = VariableFeatures(object = merged), npcs = 50)
# ---- 8. Run Harmony to remove batch effects between the two groups
merged = RunHarmony(merged,
                    group.by.vars = "group",
                    assay.use = "SCT",
                    dims.use = 1:10,
                    plot_convergence = TRUE
)
# ---- 9. Build a neighbor graph based on the Harmony embedding
merged <- FindNeighbors(merged,reduction = "harmony",dims = 1:10)
library(clustree)

# ---- 10. Perform clustering across multiple resolutions and inspect stability
# Find clusters at resolutions from 0.1 to 1.0 with a step of 0.1
merged <- FindClusters(merged, resolution = seq(0.1, 1, by = 0.1))
# Visualize how clusters split/merge across resolutions
clustree(merged, prefix = "SCT_snn_res.")

# ---- 11. Finalize clustering at resolution 0.5 and visualize with UMAP
merged <- FindClusters(merged,resolution = 0.5)
merged <- RunUMAP(merged,reduction="harmony", dims = 1:10)
DimPlot(merged,label = T,group.by = "seurat_clusters")
DimPlot(merged,label = T,group.by = "group")
# ---- 12. Prepare SCT assay for differential expression analysis
sce <- PrepSCTFindMarkers(merged, assay = "SCT", verbose = TRUE)

# ---- 13. Find marker genes for every cluster
all_markers <- FindAllMarkers(
  sce, test.use = "wilcox", only.pos = F,logfc.threshold = 0.25,min.pct = 0.25
)
# ---- 14. Keep significantly upregulated genes
sig_genes <- subset(all_markers, p_val_adj < 0.05 & avg_log2FC > 1)
# ---- 15. Select the top 5 upregulated genes per cluster 
top5 <- sig_genes %>%
  subset(avg_log2FC > 1) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# ---- 16. Visualize top marker genes with a dot plot 
DotPlot(sce, features =unique(top5$gene),
        group.by = "seurat_clusters",
        cols = c("lightgrey", "#E64B35FF"),  
        dot.scale = 6,                      
        scale = TRUE) +RotatedAxis()+
  theme(
    axis.text.x = element_text(size = 10),  
  )
  # ---- 17. Create a new metadata column for cluster labels
merged@meta.data$renew_type <- paste0("5-HT",'_',as.numeric(as.character(merged@meta.data$seurat_clusters)) + 1)

# ---- 18. Join RNA layers and save the final object
mouse_sero.combined <- JoinLayers(merged,assay = 'RNA')
# Save the combined Seurat object for downstream analysis
saveRDS(mouse_sero.combined,file = "./data/mouse_sero.combined_res0.5.rds")


