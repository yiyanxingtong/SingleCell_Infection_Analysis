基于单细胞转录组测序（scRNA-seq）的急性髓系白血病（AML）患者微环境异质性分析项目。
AML-scRNAseq-Reproduction (GSE116256)

语言与环境: R 4.4.x / Ubuntu 22.04 / RStudio
核心工具包: Seurat v5, data.table, harmony, ggplot2, ggsci
分析方法: 质量控制(QC)、多样本整合、PCA/UMAP降维、无监督聚类、细胞类型手动注释、差异表达分析(DEA)
项目背景：
本项目旨在复现经典单细胞白血病研究数据集 GSE116256。通过对多名 AML 患者骨髓样本的单细胞转录谱进行深度挖掘，识别恶性原始细胞（Blasts）与正常免疫细胞群，并定量分析患者间的细胞组成差异。

分析流程与关键成果
1. 数据预处理与质控 (Quality Control)
标准: 过滤线粒体基因占比 > 20% 及基因数 (nFeature) < 500 或 > 4000 的低质量细胞。
结果: 最终保留了约 1881 个高质量单细胞数据用于下游分析。

2. 降维与聚类 (Clustering & Visualization)
使用 PCA 进行线性降维，并通过 Elbow Plot 确定前 20 个主成分。
应用 UMAP 非线性降维算法将数据映射至二维空间。
通过 Louvain 算法（resolution=0.5）成功识别出 9 个具有显著生物学差异的细胞簇。

3. 细胞类型注释 (Cell Type Annotation)
通过 FindAllMarkers 鉴定各群特征基因，并结合经典标志物完成注释：
T/NK/B Cells: IL7R, KLRD1, MS4A1
Myeloid Lineage: CD14 (Monocytes), MAFB (Macrophages), ELANE (Neutrophils)
Malignant Blasts: GATA2, CD34 (高表达于 Cluster 4，识别为白血病原始细胞)

4. 组间定量对比 (Proportion Analysis)
结论: 通过堆叠柱状图分析发现，患者 AML210A 的原始细胞（Blasts）占比显著高于 AML1012，揭示了不同患者在疾病进展阶段的显著异质性。

项目展示
(建议在此处贴入你生成的 3 张核心图片)
UMAP Plot: 细胞分群与标注展示。
<img width="394" height="411" alt="UMP_name" src="https://github.com/user-attachments/assets/7e70eee3-0324-4f9f-8769-01ad30b79ec2" />
Dot Plot: 关键 Marker 基因在各亚群的表达强度。
<img width="607" height="411" alt="基因气泡图" src="https://github.com/user-attachments/assets/7e731a1c-ab3e-4d0c-afd1-3a21a50d51c5" />
Bar Plot: 不同患者间细胞比例的定量对比图（含百分比标注）。
<img width="667" height="495" alt="恶性原始细胞比例图" src="https://github.com/user-attachments/assets/29ba7cab-026a-4f5d-ab2a-c2e71f96ad72" />
数据来源：下载 GSM3587923 至 GSM3587926 四个原始数据文件。
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256
运行项目中的 analysis_pipeline.R 脚本。
