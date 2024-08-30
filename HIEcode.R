library(Seurat)
library(harmony)
library(Rcpp)
library(clustree)
library(ggraph)
library(RColorBrewer)
library(patchwork)
library(R.methodsS3)
library(R.utils)
library(tidyverse) 


setwd("/Users/wangjingyu/desktop/HIE")
list.files('./')
getwd()

rm(list = ls());gc()

set.seed(123) 

matrix<-'/Users/wangjingyu/desktop/HIE/scRNA'
samples<-list.files(matrix)
samples
out_dir<-'/Users/wangjingyu/desktop/HIE'
project_name<-'HIE'

scelist=lapply(samples,function(pro){
  sce=CreateSeuratObject(counts=Read10X(paste0(matrix,'/',pro)),
                         project=pro,
                         min.cells=3,
                         min.features=200)
  return(sce)
})
scelist

#
sce.all=merge(x=scelist[[1]],
              y=scelist[ -1 ], #scelist[-1]即除了scelist[1]以外的所有，即c（scelist[2],scelist[3],scelist[4],scelist[5]）
              add.cell.ids=samples,
              project=project_name)
head(sce.all@meta.data)
table(sce.all$orig.ident)

#
sce.all=PercentageFeatureSet(sce.all,'^mt-',col.name = "percent.mt")
sce.all<-subset(sce.all,subset=nFeature_RNA>400 & percent.mt<5)
#
test.seu<-sce.all
test.seu<-test.seu%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()
test.seu <- RunPCA(test.seu, npcs = 50, verbose = FALSE)

###
test.seu=test.seu %>% RunHarmony("orig.ident", plot_convergence = TRUE)
test.seu<-RunUMAP(test.seu, dim=1:30,
                  reduction ="harmony")
DimPlot(test.seu,reduction = "umap")



#
sce_cluster=test.seu

sce_cluster<-FindNeighbors(sce_cluster,reduction="harmony",
                           dims=1:30)
scRNAseq<-sce_cluster
scRNAseq <- FindClusters(scRNAseq, resolution = 0.5)
DimPlot(scRNAseq, reduction = "umap")

#cell cycle analysis
mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")
FeaturePlot(scRNAseq,features = c("Pcna","Mki67"),reduction = "umap")

#去除细胞周期的影响
scRNAseq = ScaleData(scRNAseq, vars.to.regress = c("S.Score", "G2M.Score"),
                      features = rownames(scRNAseq),split.by = "orig.ident")
saveRDS(scRNAseq,"scRNAseq.rds")
scRNAseq <- readRDS("scRNAseq.rds")

scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                           dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")

resolution<-"RNA_snn_res.0.1"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.1")

FeaturePlot(scRNAseq,features = c("Pcna","Mki67"),reduction = "umap")

#可以删除周期相关基因（可选）
dim(scRNAseq)
scRNAseq <- scRNAseq[!grepl("^mt-", rownames(scRNAseq)), ]
dim(scRNAseq) 




markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"res0.1markers_genes.csv")


#
new.cluster.ids <- c("Astrocyte","Microglia","Neuroblast","Endotheliocyte","OC","OPC","MP","Mural cell", "Fibroblast","Microglia","Ependymal cell","CPC","RBC","NC", "Neuron")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.01) 

features= c('Aqp4', 'Gja1', "Gfap",'Plpp3','Hexb','P2ry12','Tmem119', 'Cx3cr1','Tubb3', "Dcx",'Stmn2','Dlx1',"Cldn5",
            'Flt1', 'Itm2a', 'Vwf', 'Plp1','Mbp', 'Mag','Mog', 'Pdgfra', "Olig1",'Olig2', 'C1ql1','Pf4', "Mrc1",'Cd68','Lyz2',"Abcc9",
            'Pdgfrb', 'Rgs5', 'Kcnj8', 'Mylk','Myl9', 'Col1a1','Dcn', 'Pdgfra', "Col3a1", 'Dynlrb2', 'Ccdc153', 'Ak9','Ak7', 'Folr1','Kcnj13', 'Clic6', "Hbb",'Hba-a2', 'Hba-a3','Chga', "Cngt1",'Gnb3','Gnat2',"Slc17a6",
            'Nhlh2', 'Map2')
DotPlot(scRNAseq, features = unique(features)) + RotatedAxis()
#
unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.01)

#
#
Idents(scRNAseq)
#
colnames(scRNAseq@meta.data)

#
scRNAseq$celltype <- Idents(scRNAseq)
saveRDS(scRNAseq,"scRNAseq.rds")
#
scRNAseq$celltype.orig.ident <- paste(scRNAseq$celltype, scRNAseq$orig.ident, sep = "_")
scRNAseq$celltype <- Idents(scRNAseq)
Idents(scRNAseq) <- "celltype.orig.ident"

table(scRNAseq$celltype.orig.ident)


mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")




#astrocyte
scRNAseq <- readRDS("scRNAseq.rds")
scRNAseq <- subset(scRNAseq, subset = celltype == "Astrocyte" )
scRNAseq <- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)
scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")

resolution<-"RNA_snn_res.0.4"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.4")
DimPlot(scRNAseq,split.by = "orig.ident")





markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"Astromarkers_genes.csv")


new.cluster.ids <- c("A1", "A2","A3","A4","A5", "A6","A7","A8","A9","A10","A11","A12")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.05) + NoLegend()

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top10$gene),assay="RNA")+RotatedAxis()+ggplot2:::coord_flip()
DoHeatmap(subset(scRNAseq, downsample = 100), features = unique(top10$gene), size = 3,angle = 0)

mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")

unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.05)

saveRDS(scRNAseq,"AS.rds")


scRNAseq$celltype <- Idents(scRNAseq)

scRNAseq <- readRDS("AS.rds")

library(CytoTRACE2)
library(tidyverse)
scRNAseq@meta.data$CB <- rownames(scRNAseq@meta.data)
sample_CB <- scRNAseq@meta.data %>% 
  group_by(celltype) %>% 
  sample_frac(0.65)
sce3 <- subset(scRNAseq,CB %in% sample_CB$CB) 
cytotrace2_result_sce <- cytotrace2(sce3, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'mouse',
                                    seed = 1234)

annotation <- data.frame(phenotype = sce3@meta.data$celltype) %>% 
  set_rownames(., colnames(sce3))
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)

# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno

(p1+p2+p3+p4) + plot_layout(ncol = 2)

library(monocle3)

library(tidyverse)
library(dplyr)
library(patchwork)


test.seu<-scRNAseq


data <- GetAssayData(test.seu, assay ='RNA', layer = 'counts')

cell_metadata <- test.seu@meta.data

gene_annotation <-data.frame(gene_short_name = rownames(data))

rownames(gene_annotation) <-rownames(data)

cds <- new_cell_data_set(data,cell_metadata =cell_metadata, gene_metadata =gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,preprocess_method = "PCA") 

plot_cells(cds)
cds <- cluster_cells(cds) 

colnames(colData(cds))

plot_cells(cds, color_cells_by="celltype") + ggtitle('cds.umap')



cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(test.seu, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=3)

cds <- order_cells(cds)
plot_cells(cds,
           group_label_size = 6,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5
)

dea_res <- graph_test(cds, neighbor_graph = "principal_graph") %>%
  dplyr::filter(q_value < 0.05)
fwrite(dea_res, "Trajectory_genes.csv")

degs <- dea_res$gene_short_name

top_genes <- dea_res %>%
  dplyr::top_n(n = 100, morans_I) %>%
  dplyr::pull(gene_short_name) %>%
  as.character()
plot_genes_in_pseudotime(cds[top_genes, ],
                         color_cells_by = "celltype",
                         min_expr = 0.5, ncol = 6
)

p <- plot_cells(cds,
                genes = top_genes, show_trajectory_graph = FALSE,
                label_cell_groups = FALSE, label_leaves = FALSE
)
p$facet$params$ncol <- 5
p

gene_module_df <- find_gene_modules(cds[degs, ], resolution = 1e-2)
fwrite(gene_module_df, "Genes_Module.csv")

cell_group_df <- tibble::tibble(
  cell = row.names(colData(cds)),
  cell_group = colData(cds)$celltype
)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pheatmap::pheatmap(agg_mat,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   scale = "column", clustering_method = "ward.D2",
                   fontsize = 6
)

plot_cells(cds,
           genes = gene_module_df %>% dplyr::filter(module %in% c(55, 56, 34, 53,67,1,37,4,58,54,69,59,35,68,12,29,42,39,10,44,9,26,65,13,32,21)),
           group_cells_by = "partition",
           color_cells_by = "partition",
           show_trajectory_graph = FALSE
)

saveRDS(cds,file="AScds.rds")

#Microglia
scRNAseq <- readRDS("scRNAseq.rds")
scRNAseq <- subset(scRNAseq, subset = celltype == "Microglia" )
scRNAseq <- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)
scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")

resolution<-"RNA_snn_res.0.4"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.4")
DimPlot(scRNAseq,split.by = "orig.ident")





markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"MGmarkers_genes.csv")


new.cluster.ids <- c("MG1", "MG2","MG3","MG4","MG5", "MG6","MG7","MG8","MG9","MG10")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.05) + NoLegend()

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top10$gene),assay="RNA")+RotatedAxis()+ggplot2:::coord_flip()
DoHeatmap(subset(scRNAseq, downsample = 100), features = unique(top10$gene), size = 3,angle = 0)

mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")

unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.05)

saveRDS(scRNAseq,"MG.rds")


scRNAseq$celltype <- Idents(scRNAseq)

scRNAseq <- readRDS("MG.rds")

library(CytoTRACE2)
library(tidyverse)
scRNAseq@meta.data$CB <- rownames(scRNAseq@meta.data)
sample_CB <- scRNAseq@meta.data %>% 
  group_by(celltype) %>% 
  sample_frac(0.65)
sce3 <- subset(scRNAseq,CB %in% sample_CB$CB) 
cytotrace2_result_sce <- cytotrace2(sce3, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'mouse',
                                    seed = 1234)

annotation <- data.frame(phenotype = sce3@meta.data$celltype) %>% 
  set_rownames(., colnames(sce3))
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)

# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno

(p1+p2+p3+p4) + plot_layout(ncol = 2)

library(monocle3)

library(tidyverse)
library(dplyr)
library(patchwork)


test.seu<-scRNAseq


data <- GetAssayData(test.seu, assay ='RNA', slot = 'counts')

cell_metadata <- test.seu@meta.data

gene_annotation <-data.frame(gene_short_name = rownames(data))

rownames(gene_annotation) <-rownames(data)

cds <- new_cell_data_set(data,cell_metadata =cell_metadata, gene_metadata =gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,preprocess_method = "PCA") 

plot_cells(cds)
cds <- cluster_cells(cds) 

colnames(colData(cds))

plot_cells(cds, color_cells_by="celltype") + ggtitle('cds.umap')



cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(test.seu, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=3)

cds <- order_cells(cds)
plot_cells(cds,
           group_label_size = 6,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5
)

dea_res <- graph_test(cds, neighbor_graph = "principal_graph") %>%
  dplyr::filter(q_value < 0.05)
fwrite(dea_res, "Trajectory_genes.csv")

degs <- dea_res$gene_short_name

top_genes <- dea_res %>%
  dplyr::top_n(n = 100, morans_I) %>%
  dplyr::pull(gene_short_name) %>%
  as.character()
plot_genes_in_pseudotime(cds[top_genes, ],
                         color_cells_by = "celltype",
                         min_expr = 0.5, ncol = 6
)

p <- plot_cells(cds,
                genes = top_genes, show_trajectory_graph = FALSE,
                label_cell_groups = FALSE, label_leaves = FALSE
)
p$facet$params$ncol <- 5
p

gene_module_df <- find_gene_modules(cds[degs, ], resolution = 1e-2)
fwrite(gene_module_df, "Genes_Module.csv")

cell_group_df <- tibble::tibble(
  cell = row.names(colData(cds)),
  cell_group = colData(cds)$celltype
)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

pheatmap::pheatmap(agg_mat,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   scale = "column", clustering_method = "ward.D2",
                   fontsize = 6
)

plot_cells(cds,
           genes = gene_module_df %>% dplyr::filter(module %in% c(55, 56, 34, 53,67,1,37,4,58,54,69,59,35,68,12,29,42,39,10,44,9,26,65,13,32,21)),
           group_cells_by = "partition",
           color_cells_by = "partition",
           show_trajectory_graph = FALSE
)

saveRDS(cds,file="MGcds.rds")



#Neuroblast
scRNAseq <- readRDS("scRNAseq.rds")
scRNAseq <- subset(scRNAseq, subset = celltype == "Neuroblast" )
scRNAseq <- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)
scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")

resolution<-"RNA_snn_res.0.4"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.4")
DimPlot(scRNAseq,split.by = "orig.ident")

markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"NBtromarkers_genes.csv")


new.cluster.ids <- c("NB1", "NB2","NB3","NB4","NB5", "NB6","NB7","NB8","NB9","NB10","NB11","NB12")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.05) + NoLegend()

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top10$gene),assay="RNA")+RotatedAxis()+ggplot2:::coord_flip()
DoHeatmap(subset(scRNAseq, downsample = 100), features = unique(top10$gene), size = 3,angle = 0)

scRNAseq$celltype.orig.ident <- paste(scRNAseq$celltype, scRNAseq$orig.ident, sep = "_")
scRNAseq$celltype <- Idents(scRNAseq)
Idents(scRNAseq) <- "celltype.orig.ident"
table(scRNAseq$celltype.orig.ident)

mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")

unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.05)

saveRDS(scRNAseq,"NB.rds")


scRNAseq$celltype <- Idents(scRNAseq)

scRNAseq <- readRDS("NB.rds")


#OPC olgi
scRNAseq <- readRDS("scRNAseq.rds")
scRNAseq <- subset(scRNAseq, subset = celltype == "OPC"|celltype == "OC" )
scRNAseq <- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)
scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")

resolution<-"RNA_snn_res.0.4"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.4")
DimPlot(scRNAseq,group.by = "celltype")
DimPlot(scRNAseq,split.by = "orig.ident")

markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"OPCOCmarkers_genes.csv")


new.cluster.ids <- c("OPC1", "OC1","OC2","OC3","OPC2", "OC4","OC5","OPC3","OPC4","OPC/OC","OPC5","OPC6")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.05) 

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top5$gene),assay="RNA")+RotatedAxis()+ggplot2:::coord_flip()
DoHeatmap(subset(scRNAseq, downsample = 100), features = unique(top10$gene), size = 3,angle = 0)

scRNAseq$celltype.orig.ident <- paste(scRNAseq$celltype, scRNAseq$orig.ident, sep = "_")
scRNAseq$celltype <- Idents(scRNAseq)
Idents(scRNAseq) <- "celltype.orig.ident"
table(scRNAseq$celltype.orig.ident)

mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")

unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.05)

saveRDS(scRNAseq,"OPCOC.rds")


scRNAseq$celltype <- Idents(scRNAseq)

scRNAseq <- readRDS("OPCOC.rds")

DotPlot(scRNAseq, features = unique(features)) 


#endo
scRNAseq <- readRDS("scRNAseq.rds")
scRNAseq <- subset(scRNAseq, subset = celltype == "Endotheliocyte")
scRNAseq <- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)
scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")



mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")
FeaturePlot(scRNAseq,features = c("Pcna","Mki67"),reduction = "umap")

#去除细胞周期的影响
scRNAseq = ScaleData(scRNAseq, vars.to.regress = c("S.Score", "G2M.Score"),
                     features = rownames(scRNAseq),split.by = "orig.ident")

scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")

resolution<-"RNA_snn_res.0.2"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.2")
DimPlot(scRNAseq,split.by = "orig.ident",label = TRUE)

FeaturePlot(scRNAseq,features = c("Pcna","Mki67"),reduction = "umap")


markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"ECmarkers_genes.csv")


new.cluster.ids <- c("VEC-capillary1", "VEC-arterial","VEC-venous","VEC-inter","LEC-valve1", "VEC-capillary2","LEC-major1","Proliferative LEC","LEC-valve2","LEC-valve3","Mixed cells","Novel VEC","LEC-major2")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.05) 

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top5$gene),assay="RNA")+ggplot2:::coord_flip()
DoHeatmap(subset(scRNAseq, downsample = 100), features = unique(top10$gene), size = 3,angle = 0)

features= c('Ca4', 'Rgcc',"Sox17",
            'Gja4', 'Lyve1', 'Fabp5','Prss23', 'Cldn11', "Scg3", 'Cd24',
            'Ackr4','Ptgs2','Cxcl12','Igfbp5')
DotPlot(scRNAseq, features = unique(features)) +ggplot2:::coord_flip()

scRNAseq$celltype <- Idents(scRNAseq)
scRNAseq$celltype.orig.ident <- paste(scRNAseq$celltype, scRNAseq$orig.ident, sep = "_")
scRNAseq$celltype <- Idents(scRNAseq)
Idents(scRNAseq) <- "celltype.orig.ident"
table(scRNAseq$celltype.orig.ident)


unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.05)

saveRDS(scRNAseq,"EC.rds")


#FIB
scRNAseq <- readRDS("scRNAseq.rds")
scRNAseq <- subset(scRNAseq, subset = celltype == "Fibroblast")
scRNAseq <- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)
scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 


scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")

resolution<-"RNA_snn_res.0.2"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.2")
DimPlot(scRNAseq,split.by = "orig.ident",label = TRUE)



markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"FIBmarkers_genes.csv")

new.cluster.ids <- c("FIB1", "FIB2","FIB3","FIB4","FIB5", "FIB6","FIB7", "FIB8","FIB9")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.25) 

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top5$gene),assay="RNA")+ggplot2:::coord_flip()
DoHeatmap(subset(scRNAseq, downsample = 100), features = unique(top10$gene), size = 3,angle = 0)

features= c('Sparcl1',"Wnt6",
            'Shisa3', 'Depdc1', 'Tnnt2','Sox10', 'Tmem229a', "Gpr34", 'Igfbpl1')
DotPlot(scRNAseq, features = unique(features)) +ggplot2:::coord_flip()
VlnPlot(scRNAseq, features = unique(features)) 

scRNAseq$celltype <- Idents(scRNAseq)
scRNAseq$celltype.orig.ident <- paste(scRNAseq$celltype, scRNAseq$orig.ident, sep = "_")
scRNAseq$celltype <- Idents(scRNAseq)
Idents(scRNAseq) <- "celltype.orig.ident"
table(scRNAseq$celltype.orig.ident)


unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.05)

saveRDS(scRNAseq,"FIB.rds")


#MP
scRNAseq <- readRDS("scRNAseq.rds")
scRNAseq <- subset(scRNAseq, subset = celltype == "MP")
scRNAseq <- NormalizeData(scRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq <- FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = 2000)
scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 


#cell cycle analysis
mouse_cell_cycle_genes <- readRDS("/Users/wangjingyu/desktop/code/mouse_cell_cycle_genes.rds")
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNAseq <- CellCycleScoring(scRNAseq, s.features = s.genes, g2m.features = g2m.genes)
scRNAseq@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
DimPlot(scRNAseq,reduction = "umap" ,group.by = "Phase")
FeaturePlot(scRNAseq,features = c("Pcna","Mki67"),reduction = "umap")

#去除细胞周期的影响
scRNAseq = ScaleData(scRNAseq, vars.to.regress = c("S.Score", "G2M.Score"),
                     features = rownames(scRNAseq),split.by = "orig.ident")

scRNAseq <- RunPCA(scRNAseq, npcs = 50, verbose = FALSE)
scRNAseq=scRNAseq %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNAseq <- RunUMAP(scRNAseq,dims = 1:30,reduction = "harmony") 

scRNAseq<-FindNeighbors(scRNAseq,reduction="harmony",
                        dims=1:30)
for (i in c(0.1,0.2,0.3,0.4,0.5,0.8,1.0)) {
  scRNAseq<-FindClusters(scRNAseq,resolution=i,algorithm=1)
}
head(scRNAseq@meta.data)

clustree(scRNAseq@meta.data,prefix="RNA_snn_res.")


resolution<-"RNA_snn_res.0.3"
Idents(scRNAseq)=resolution
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.3")
DimPlot(scRNAseq,split.by = "orig.ident",label = TRUE)
DimPlot(scRNAseq,group.by = "RNA_snn_res.0.3",label = TRUE)


markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"MPmarkers_genes.csv")

new.cluster.ids <- c("1", "2","3","4","5", "6","7", "8","9", "10","11", "12","3")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", label = TRUE, pt.size = 0.25) 

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top5$gene),assay="RNA")+ggplot2:::coord_flip()
DoHeatmap(subset(scRNAseq, downsample = 100), features = unique(top10$gene), size = 3,angle = 0)

features= c('Mrc1',"Pf4","Aif1",'Spp1',"Treml4",'Clec7a',
            'Cd3e', 'Cd2', 'S100a9','Mmp9', 'Fcer1a', "Cd19", 'Cd74',
            'Cd69', 'Icos', "Ccr9", 'Sox9', 'Olig1','Tubb3')
DotPlot(scRNAseq, features = unique(features)) +ggplot2:::coord_flip()
VlnPlot(scRNAseq, features = unique(features)) 

scRNAseq$celltype <- Idents(scRNAseq)
scRNAseq$celltype.orig.ident <- paste(scRNAseq$celltype, scRNAseq$orig.ident, sep = "_")
scRNAseq$celltype <- Idents(scRNAseq)
Idents(scRNAseq) <- "celltype.orig.ident"
table(scRNAseq$celltype.orig.ident)


unique(scRNAseq$orig.ident)
DimPlot(scRNAseq,split.by = "orig.ident", label = TRUE, pt.size = 0.05)

saveRDS(scRNAseq,"MP.rds")
scRNAseq <- readRDS("MP.rds")
