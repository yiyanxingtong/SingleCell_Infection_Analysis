library(Seurat)
library(data.table)
library(magrittr)
setwd("E:/xiangmu/date")
aml1012_counts <- fread("GSM3587923_AML1012-D0.dem.txt.gz")
aml1012_counts <- as.data.frame(aml1012_counts)
rownames(aml1012_counts) <- aml1012_counts[[1]]
aml1012_counts <- aml1012_counts[, -1]
aml1012_anno <- fread("GSM3587924_AML1012-D0.anno.txt.gz")
aml1012_anno <- as.data.frame(aml1012_anno)
rownames(aml1012_anno) <- aml1012_anno[[1]]
obj1012 <- CreateSeuratObject(counts = aml1012_counts, 
                              project = "AML1012", 
                              meta.data = aml1012_anno)
aml210_counts <- fread("GSM3587925_AML210A-D0.dem.txt.gz")
aml210_counts <- as.data.frame(aml210_counts)
rownames(aml210_counts) <- aml210_counts[[1]]
aml210_counts <- aml210_counts[, -1]
aml210_anno <- fread("GSM3587926_AML210A-D0.anno.txt.gz")
aml210_anno <- as.data.frame(aml210_anno)
rownames(aml210_anno) <- aml210_anno[[1]]
obj210 <- CreateSeuratObject(counts = aml210_counts, 
                             project = "AML210A", 
                             meta.data = aml210_anno)
combined_obj <- merge(obj1012, y = obj210, add.cell.ids = c("AML1012", "AML210A"))
table(combined_obj$orig.ident)
combined_obj[["percent.mt"]] <- PercentageFeatureSet(combined_obj, pattern = "^MT-")
VlnPlot(combined_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

combined_obj <- subset(combined_obj, 
                       subset = nFeature_RNA > 500 & 
                         nFeature_RNA < 4000 & 
                         percent.mt < 20)
#过滤
VlnPlot(combined_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

combined_obj <- NormalizeData(combined_obj)
combined_obj <- FindVariableFeatures(combined_obj, selection.method = "vst", nfeatures = 2000)
combined_obj <- ScaleData(combined_obj)
combined_obj <- RunPCA(combined_obj, features = VariableFeatures(object = combined_obj))
ElbowPlot(combined_obj)
ncol(combined_obj)

# 1. 运行 PCA
combined_obj <- RunPCA(combined_obj, features = VariableFeatures(object = combined_obj))
#计算邻近距离 (FindNeighbors)
combined_obj <- FindNeighbors(combined_obj, dims = 1:20)
# 2. 聚类 (FindClusters) 
# resolution 参数决定了分的群数，0.5 是常用默认值。数字越大，群分得越细。
combined_obj <- FindClusters(combined_obj, resolution = 0.5)
combined_obj <- RunUMAP(combined_obj, dims = 1:20)
DimPlot(combined_obj, reduction = "umap", label = TRUE, pt.size = 1)

# 寻找所有群的 Marker 基因（建议只找上调且显著的）
# logfc.threshold: 表达差异倍数（0.25 倍以上）
# min.pct: 至少在 25% 的细胞中表达
combined_obj <- JoinLayers(combined_obj)
all_markers <- FindAllMarkers(combined_obj, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# 查看每个群前 2 名的基因
library(dplyr)
top2 <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
print(top2)

# 1. 定义新名字（顺序对应 Cluster 0-8）
new_labels <- c("Stroma", "Neutrophils", "Macrophages", "Monocytes", "Blasts", 
                "T Cells", "NK Cells", "Endothelial", "B Cells")

# 2. 赋予名字
names(new_labels) <- levels(combined_obj)
combined_obj <- RenameIdents(combined_obj, new_labels)

# 3. 画出带名字的 UMAP 图
library(ggplot2)
DimPlot(combined_obj, reduction = "umap", label = TRUE, pt.size = 1) + 
  ggtitle("Cell Type Annotation for AML Dataset")
DotPlot(combined_obj, features = unique(top2$gene)) + RotatedAxis
DotPlot(combined_obj, features = unique(top2$gene)) + RotatedAxis()


# 1. 提取每个样本中各细胞类型的数量
cell_counts <- table(Idents(combined_obj), combined_obj$orig.ident)
cell_props <- prop.table(cell_counts, margin = 2) # 计算比例

# 2. 转换为数据框方便绘图
df <- as.data.frame(cell_props)
colnames(df) <- c("CellType", "Sample", "Proportion")
library(scales)
library(ggsci)

library(ggplot2)
library(dplyr)
library(ggsci)
library(scales)

# 1. 确保 CellType 是因子，并保持你想要的堆叠顺序（这能保证文字和色块匹配）
df$CellType <- factor(df$CellType, levels = rev(c("Stroma", "Neutrophils", "Macrophages", 
                                                  "Monocytes", "Blasts", "T Cells", 
                                                  "NK Cells", "Endothelial", "B Cells")))

# 2. 绘图：直接使用 position_stack 自动计算位置
ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
  # 柱状图
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white", linewidth = 0.2) +
  # 核心：使用 position_stack(vjust = 0.5) 自动将标签放在每个色块的正中心
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            color = "white", 
            size = 3.5,
            fontface = "bold") + 
  # 配色与主题优化
  scale_fill_npg() + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  ) +
  labs(
    title = "Quantitative Composition of Bone Marrow Niche",
    subtitle = "Analysis of Cell Proportion across Patient Samples",
    y = "Cell Proportion (%)", 
    x = "Patient Sample",
    fill = "Cell Lineage"
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0, 0.05, 0))
