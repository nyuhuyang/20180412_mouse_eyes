########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
mouse_3eyes_raw <- list()
samples <- c("129_B6","129_B6_aged","B6")
projects <- c("EC-IB-4698","EC-IB-4867","EC-IB-4698")
groups <- c("young.129_B6", "aged.129_B6","young.B6")
conditions <- c("young", "aged","young")

for(i in 1:length(samples)){
    mouse_3eyes_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                     samples[i],"/outs/filtered_gene_bc_matrices/mm10/"))
    colnames(mouse_3eyes_raw[[i]]) <- paste0(groups[i],
                                            ".",colnames(mouse_3eyes_raw[[i]]))
}
mouse_3eyes_list <- lapply(mouse_3eyes_raw, CreateSeuratObject,
                      min.cells = 3,
                      min.genes = 200,
                      project = "DropSeq")
mouse_3eyes_list <- lapply(mouse_3eyes_list, FilterCells, 
                      subset.names = "nGene", 
                      low.thresholds = 200, 
                      high.thresholds = Inf)
mouse_3eyes_list <- lapply(mouse_3eyes_list, NormalizeData)
mouse_3eyes_list <- lapply(mouse_3eyes_list, FindVariableGenes, do.plot = FALSE)
mouse_3eyes_list <- lapply(mouse_3eyes_list, ScaleData)
for(i in 1:length(groups)) mouse_3eyes_list[[i]]@meta.data$conditions <- groups[i]

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
genes.use <- lapply(mouse_3eyes_list, function(x) head(rownames(x@hvg.info), 1000))
genes.use <- unique(unlist(genes.use))
for(i in 1:length(samples)){
        genes.use <- intersect(genes.use, rownames(mouse_3eyes_list[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size

#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
mouse_3eyes <- RunMultiCCA(object.list = mouse_3eyes_list, 
                      genes.use = genes.use,
                      niter = 25, num.ccs = 30,
                      standardize =TRUE)
save(mouse_3eyes, file = "./data/mouse_3eyes_alignment.Rda")
remove(mouse_3eyes_list)
remove(mouse_3eyes_raw)
# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = mouse_3eyes, reduction.use = "cca", 
              group.by = "conditions", pt.size =1, 
              do.return = TRUE)
p2 <- VlnPlot(object = mouse_3eyes, features.plot = "CC1", 
              group.by = "conditions", do.return = TRUE)
#png(paste0(pwd,'/output/Dim_violin_plots.png'))
plot_grid(p1, p2)
#dev.off()
#png(paste0(pwd,'/output/MetageneBicorPlot.png'))
p3 <- MetageneBicorPlot(mouse_3eyes, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE)
p3
#dev.off()

#======1.3 QC ==================================
# Run rare non-overlapping filtering
mouse_3eyes <- CalcVarExpRatio(object = mouse_3eyes, reduction.type = "pca",
                          grouping.var = "conditions", dims.use = 1:15)
mouse_3eyes <- SubsetData(mouse_3eyes, subset.name = "var.ratio.pca",accept.low = 0.5)

mito.genes <- grep(pattern = "^mt-", x = rownames(x = mouse_3eyes@data), value = TRUE)
percent.mito <- Matrix::colSums(mouse_3eyes@raw.data[mito.genes, ])/Matrix::colSums(mouse_3eyes@raw.data)
mouse_3eyes <- AddMetaData(object = mouse_3eyes, metadata = percent.mito, col.name = "percent.mito")
mouse_3eyes <- ScaleData(object = mouse_3eyes, genes.use = genes.use, display.progress = FALSE, 
                                vars.to.regress = "percent.mito", do.par = TRUE, num.cores = 4)

#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

mouse_3eyes <- AlignSubspace(object = mouse_3eyes, grouping.var = "conditions", 
                        dims.align = 1:15)
#Now we can run a single integrated analysis on all cells!
mouse_3eyes <- FindClusters(object = mouse_3eyes, reduction.type = "cca.aligned", dims.use = 1:15, 
                       resolution = 0.8, force.recalc = TRUE, save.SNN = TRUE)
mouse_3eyes <- RunTSNE(object = mouse_3eyes, reduction.use = "cca.aligned", dims.use = 1:15, 
                  dim.embed = 2, do.fast = TRUE)
mouse_3eyes <- RunPCA(object = mouse_3eyes, pc.genes = mouse_3eyes@var.genes,pcs.compute = 30, do.fast = TRUE)
VizPCA(object = mouse_3eyes, pcs.use = 1:2)
p1 <- TSNEPlot(mouse_3eyes, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(mouse_3eyes, do.label = F, do.return = T, pt.size = 1)
#png(paste0(pwd,'/output/TSNES_plots.png'))
plot_grid(p1, p2)
#dev.off()
p4 <- PCElbowPlot(object = mouse_3eyes, num.pc = 30)
p5 <- p3 + theme(legend.position="none") #change legend position
plot_grid(p4, p5)
# Now, we annotate the clusters as before based on canonical markers.
#png(paste0(pwd,'/output/TheClusters.png'))
TSNEPlot(object = mouse_3eyes,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
        ggtitle("TSEN Plot of all clusters")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
#dev.off()
save(mouse_3eyes, file = "./data/mouse_3eyes_alignment.Rda")

#======1.5 Compare clusters for each dataset =======================
cell.all <- FetchData(mouse_3eyes,"conditions")
cell.subsets <- lapply(groups, function(x) 
        rownames(cell.all)[cell.all$conditions == x])

mouse_3eyes.subsets <- list()
for(i in 1:length(groups)){
        mouse_3eyes.subsets[[i]] <- SubsetData(mouse_3eyes, cells.use =cell.subsets[[i]])
}

table(mouse_3eyes.subsets[[1]]@ident)

p <- list()
for(i in 1:length(groups)){
        p[[i]] <- TSNEPlot(object = mouse_3eyes.subsets[[i]],do.label = F, group.by = "ident", 
                           do.return = TRUE, no.legend = TRUE,
                           pt.size = 1,label.size = 4 )+
                ggtitle(groups[i])+
                theme(text = element_text(size=20),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
}
do.call(plot_grid, p)

#======1.6 Further subdivisions with FIt-SNE =======================
lnames = load(file = "./data/mouse_3eyes_alignment.Rda")
lnames

mouse_3eyes_FItSNE <- FindVariableGenes(object = mouse_3eyes, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
hv.genes <- head(rownames(mouse_3eyes_FItSNE@hvg.info), 1000)
mouse_3eyes_FItSNE <- RunPCA(object = mouse_3eyes_FItSNE, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = mouse_3eyes_FItSNE, num.pc = 100)
PCHeatmap(mouse_3eyes_FItSNE, pc.use = c(1:3,73:75), cells.use = 500, do.balanced = TRUE)
mouse_3eyes_FItSNE <- FindClusters(object = mouse_3eyes_FItSNE, reduction.type = "pca", dims.use = 1:75, resolution = 3, 
                    save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)
mouse_3eyes_FItSNE <- RunTSNE(object = mouse_3eyes_FItSNE, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE", 
               nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", 
               fast_tsne_path = "/Users/yah2014/src/FIt-SNE/bin/fast_tsne", 
               max_iter = 2000)
library(cowplot)
p1 <- DimPlot(object = mouse_3eyes_FItSNE, reduction.use = "FItSNE", 
              no.legend = TRUE, do.return = TRUE, do.label = TRUE,
              vector.friendly = FALSE, pt.size = 1,label.size = 8 ) + 
        ggtitle("FItSNE plot of all clusters") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p2 <- DimPlot(object = mouse_3eyes_FItSNE, reduction.use = "FItSNE", no.legend = FALSE, group.by = "conditions", 
              do.return = TRUE, vector.friendly = FALSE, pt.size = 1) + 
        ggtitle("FItSNE plot of all conditions") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
FeaturePlot(mouse_3eyes_FItSNE, c("Ihh","Gli1","Ptch1","Hhip"),reduction.use = "FItSNE", dark.theme = TRUE, 
            pt.size = 1, vector.friendly = FALSE)
plot_grid(p1, p2)
table(mouse_3eyes_FItSNE@ident)

save(mouse_3eyes_FItSNE, file = "./data/mouse_3eyes_FIt-SNE.Rda")
