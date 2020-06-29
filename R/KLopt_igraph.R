#test performance on zachary data
library(igraph)
gzachary <- read_graph('EdgeLists zachary.tsv', format = c("edgelist"), directed = FALSE)
coords_zachary = layout_with_fr(gzachary)
plot(gzachary)

# color by KLoptimization
label_gzachary <- read.table("FoundComms1 zachary.tsv", sep = '\t', header = FALSE)
plot(gzachary, vertex.color=label_gzachary$V2+1, layout=coords_zachary)
# graphjs(gzachary, vertex.color=c("green", "blue")[label_gzachary$V2+1], vertex.shape="sphere")


# color by Leading eigenvector partition
lec_zachary <- cluster_leading_eigen(gzachary)
plot(lec_zachary, gzachary, layout=coords_zachary)


# color by cluster_fast_greedy
c1 = cluster_louvain(gzachary)
# modularity measure
modularity(c1)
# modularity matrix
B = modularity_matrix(gzachary, membership(c1))
plot(c1, gzachary, layout=coords_zachary)
plot(gzachary, vertex.color=membership(c1), layout=coords_zachary)

# #leiden
# adjacency_matrix <- igraph::as_adjacency_matrix(gzachary)
# partition <- leiden(adjacency_matrix)
# table(partition)
# library("RColorBrewer")
# node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
# plot(graph_object, vertex.color = node.cols)
#test performance on zachary data
gblog <- read_graph('EdgeLists.tsv', format = c("edgelist"), directed = FALSE)
coords_blog = layout_with_fr(gblog)
plot(gblog)

# color by KLoptimization
label_blog <- read.table("FoundComms1.tsv", sep = '\t', header = FALSE)
plot(gblog, vertex.color=label_blog$V2+1, layout=coords_blog)

# color by Leading eigenvector partition
lec_blog <- cluster_leading_eigen(gblog)
plot(lec_blog, gblog), layout=coords_blog)
rglplot(gblog, layout=coords_blog)



#test performance on Seurat PBMC data
library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5, method = "igraph")
DimPlot(pbmc)

input <- as.matrix(pbmc[["RNA_snn"]])
rownames(input) = colnames(input) <- 1:2700
graph_pbmc <- graph.adjacency(input, weighted = T, mode="undirected")
pbmc_edgelist <- as.data.frame(as_edgelist(graph_pbmc))
df2 <- mutate_all(pbmc_edgelist, function(x) as.numeric(as.character(x)))
write.table(df2, "pbmc_edgelist.txt", row.names = FALSE, col.names =  FALSE)

graph_pbmc <- read_graph('pbmc_edgelist.txt', format = c("edgelist"), directed = FALSE)
coords_pbmc = layout_with_fr(graph_pbmc)
lec_pbmc <- cluster_leading_eigen(graph_pbmc)
table(lec_pbmc$membership)
#   1   2   3   4   5   6   7   8   9  10  11  12  13 
# 377 173  60 145 381  81 347 265 276  12 322 164  97
plot(lec_pbmc, graph_pbmc, layout=coords_pbmc)

coords_pbmc = layout_with_fr(graph_pbmc)
plot(graph_pbmc)

# color by KLoptimization
label_pbmc <- read.table("FoundComms1_pbmc.tsv", sep = '\t', header = FALSE)
table(label_pbmc$V2)
  0   1   2   3   4   5   6   7   8 
349 147 192 303 430 251 348 258 422
plot(graph_pbmc, vertex.color=label_pbmc$V2+1, layout=coords_pbmc, vertex.size=1, vertex.label=NA)
# graphjs(gzachary, vertex.color=c("green", "blue")[label_gzachary$V2+1], vertex.shape="sphere")
table(pbmc[['seurat_clusters']]$ seurat_clusters, label_pbmc$V2)
     0   1   2   3   4   5   6   7   8
  0   0   0   0   0   1   0 312   0 378
  1   0   0   1   0   0 236   0 256   0
  2   0   0   0   1 429   0   0   0  44
  3 349   0   0   0   0   0   0   0   0
  4   0   7   0 294   0   0  36   0   0
  5   0   0 157   0   0   0   0   0   0
  6   0 140   0   8   0   0   0   0   0
  7   0   0  34   0   0   0   0   2   0
  8   0   0   0   0   0  15   0   0   0
#leiden
# partition <- leiden(graph_pbmc)
# table(partition)
# library("RColorBrewer")
# node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
# plot(graph_object, vertex.color = node.cols)


# splatter simulation


