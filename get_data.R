
# Based on https://satijalab.org/seurat/articles/spatial_vignette

#install.packages("Seurat")
#remotes::install_github('satijalab/seurat-data')

library(Seurat)
library(SeuratData)

library(patchwork)
library(tidyverse)

InstallData("stxBrain")

brain <- LoadData("stxBrain", type = "anterior1")

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

centroids <- brain$anterior1$centroids |> as.data.frame()

brain@images$anterior1@image |> dim()
brain@images$anterior1@image |> range()

brain@images$anterior1@image |> png::writePNG("anterior1.png")

sf <- brain@images$anterior1@scale.factors$lowres
plot(centroids$y * sf, 600-centroids$x *sf, xlim=c(0,600), ylim=c(0,600), asp=1)
range(centroids$x * brain@images$anterior1@scale.factors$lowres)
range(centroids$y * brain@images$anterior1@scale.factors$lowres)

range(brain$umap@cell.embeddings[,1]*25+300)
range(brain$umap@cell.embeddings[,2]*25+300)


sf <- brain@images$anterior1@scale.factors$lowres
centroids <- brain$anterior1$centroids |> as.data.frame()
#centroids <- GetTissueCoordinates(brain)
dataset <- tibble(
    imageX = centroids$y * sf, 
    imageY = centroids$x *sf,
    umapX = brain$umap@cell.embeddings[,1]*25+300,
    umapY = brain$umap@cell.embeddings[,2]*25+300,
    cluster = brain@meta.data$seurat_clusters,
    nCount = brain@meta.data$nCount_Spatial)

write_csv(dataset, "spots.csv")

json_script <- paste0("let spots = ", jsonlite::toJSON(dataset, pretty=TRUE), ";")
writeLines(json_script, "spots.js")