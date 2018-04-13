########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#==3.1 Compare TSNE vs TItSNE =====
# Rename ident
lnames = load(file = "./data/mouse_3eyes_alignment.Rda")
lnames = load(file = "./data/mouse_3eyes_cluster_3.Rda")
lnames
# rename mouse_3eyes ======
table(mouse_3eyes@ident)
idents <- as.data.frame(table(mouse_3eyes@ident))
old.ident.ids1 <- idents$Var1
new.cluster.ids1 <- c("Mesenchymal cells",
                     "Retinal Pigment Epithelium",
                     "Pericytes",
                     "Smooth muscle cells",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Mesenchymal cells",
                     "Endothelial cells",
                     "Myelinating Schwann cells",
                     "Monocytes",
                     "Monocytes",
                     "T cells",
                     "Melanocytes",
                     "Myelinating Schwann cells")
mouse_3eyes@ident <- plyr::mapvalues(x = mouse_3eyes@ident,
                                     from = old.ident.ids1,
                                     to = new.cluster.ids1)
table(mouse_3eyes@ident)

# rename mouse_3eyes_FItSNE ======
table(mouse_3eyes_cluster_3@ident)
idents <- as.data.frame(table(mouse_3eyes_cluster_3@ident))
old.ident.ids2 <- idents$Var1
new.cluster.ids2 <- c("Pericytes",
                     "Mesenchymal cells",
                     "Retinal Pigment Epithelium",
                     "Mesenchymal cells",
                     "Mesenchymal cells",
                     "Mesenchymal cells",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Smooth muscle cells",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Smooth muscle cells",
                     "Smooth muscle cells",
                     "Mesenchymal cells",
                     "Endothelial cells",
                     "Schwann cells",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Monocytes",
                     "Retinal Pigment Epithelium",
                     "CD 14 Monocytes",
                     "Endothelial cells",
                     "Hematopoietic SCs",
                     "T & NK cells",
                     "Retinal Pigment Epithelium",
                     "Smooth muscle cells",
                     "Endothelial cells",
                     "Monocytes & NK cells",
                     "Melanocytes",
                     "Monocytes & NK cells",
                     "Schwann cells")
mouse_3eyes_cluster_3@ident <- plyr::mapvalues(x = mouse_3eyes_cluster_3@ident,
                                            from = old.ident.ids2,
                                            to = new.cluster.ids2)
table(mouse_3eyes_cluster_3@ident)
#====CD 34 =======
p1 <- SingleFeaturePlot.1(object = mouse_3eyes, feature = "Cd34", 
                          reduction.use = "tsne",
                          do.return = TRUE, 
                          cols.use = c("lightgrey","blue"), pt.size = 0.5) + 
        ggtitle("Cd34 in TSNE plot resolution 0.8") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

p2 <- SingleFeaturePlot.1(object = mouse_3eyes_cluster_3, feature = "Cd34", 
                          reduction.use = "tsne",
                          do.return = TRUE, 
                          cols.use = c("lightgrey","blue"), pt.size = 0.5) + 
        ggtitle("Cd34 in TSNE plot resolution 3") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2)
#====3.2 Further subdivisions Mesenchymal cells =======
Mesenchymal <- SubsetData(object = mouse_3eyes,
                        ident.use = old.ident.ids1[(new.cluster.ids1 %in% "Mesenchymal cells")])
table(Mesenchymal@ident)
p3 <- DimPlot(object = Mesenchymal, reduction.use = "tsne", 
        no.legend = TRUE, do.return = TRUE, do.label = TRUE,
        vector.friendly = FALSE, pt.size = 1,label.size = 5 ) + 
        ggtitle("Original Resolution = 0.8") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
Mesenchymal <- FindVariableGenes(object = Mesenchymal, mean.function = ExpMean, dispersion.function = LogVMR, 
                                        do.plot = FALSE)
hv.genes <- head(rownames(Mesenchymal@hvg.info), 1000)
Mesenchymal <- RunPCA(object = Mesenchymal, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, 
                             pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = Mesenchymal, num.pc = 100)
PCHeatmap(Mesenchymal, pc.use = c(1:3,53:55), cells.use = 500, do.balanced = TRUE)
Mesenchymal <- RunTSNE(object = Mesenchymal, reduction.use = "pca", dims.use = 1:55, 
                       dim.embed = 2, do.fast = TRUE)

# find clusters
Mesenchymal_pca <- FindClusters(object = Mesenchymal, reduction.type = "pca", dims.use = 1:55,
                                resolution = 1,save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)

Mesenchymal_tsne <- FindClusters(object = Mesenchymal, reduction.type = "tsne",dims.use = 1:2,
                                 resolution = 1,save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)
Mesenchymal <- RunTSNE(object = Mesenchymal, reduction.use = "pca", dims.use = 1:55,
                       tsne.method = "FIt-SNE",nthreads = 4, reduction.name = "FItSNE",
                       reduction.key = "FItSNE_",fast_tsne_path = "/Users/yah2014/src/FIt-SNE/bin/fast_tsne", 
                       max_iter = 2000)
Mesenchymal_fitsne <- FindClusters(object = Mesenchymal, reduction.type = "FItSNE",dims.use = 1:2,
                                 resolution = 1,save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)


p4 <- DimPlot(object = Mesenchymal_pca, reduction.use = "pca", 
              no.legend = TRUE, do.return = TRUE, do.label = TRUE,
              vector.friendly = FALSE, pt.size = 1,label.size = 5 ) + 
        ggtitle("PCA Resolution = 1") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p5 <- DimPlot(object = Mesenchymal_tsne, reduction.use = "tsne", 
              no.legend = TRUE, do.return = TRUE, do.label = TRUE,
              vector.friendly = FALSE, pt.size = 1,label.size = 5 ) + 
        ggtitle("tSNE Resolution = 1") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p6 <- DimPlot(object = Mesenchymal_fitsne, reduction.use = "FItSNE", 
              no.legend = TRUE, do.return = TRUE, do.label = TRUE,
              vector.friendly = FALSE, pt.size = 1,label.size = 5 ) + 
        ggtitle("FItSNE Resolution = 1") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p3, p4,p5,p6)

#======3.3 Compare clusters for each dataset =======================
cell.all <- FetchData(mouse_3eyes_FItSNE,"conditions")
groups <- c("young.129_B6", "aged.129_B6","young.B6")
cell.subsets <- lapply(groups, function(x) 
        rownames(cell.all)[cell.all$conditions == x])

mouse_3eyes_FItSNE.subsets <- list()
for(i in 1:length(groups)){
        mouse_3eyes_FItSNE.subsets[[i]] <- SubsetData(mouse_3eyes_FItSNE, cells.use =cell.subsets[[i]])
}

table(mouse_3eyes_FItSNE.subsets[[1]]@ident)

p <- list()
for(i in 1:length(groups)){
        p[[i]] <- DimPlot(object = mouse_3eyes_FItSNE.subsets[[i]], reduction.use = "FItSNE", 
                          do.label = F, group.by = "ident", 
                           do.return = TRUE, no.legend = TRUE,
                           pt.size = 1,label.size = 4 )+
                ggtitle(groups[i])+
                theme(text = element_text(size=20),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
}
do.call(plot_grid, p)

DimPlot(object = mouse_3eyes_FItSNE.subsets[[i]], reduction.use = "FItSNE", 
        do.label = F, group.by = "ident", 
        do.return = TRUE, no.legend = FALSE,
        pt.size = 1,label.size = 4 )+
        theme(text = element_text(size=20))     #larger text including legend title
#===== A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.

freq_table <- prop.table(x = table(mouse_3eyes_FItSNE@ident, mouse_3eyes_FItSNE@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table


#==3.2 SplitDotPlotGG======
# The SplitDotPlotGG function can be useful for viewing conserved cell type markers
# across conditions, showing both the expression level and the percentage of cells
# in a cluster expressing any given gene. 
# Here we plot 2-3 strong marker genes for each of our 13 clusters.
markers.to.plot <- c("Pmel","Dcn", "Laptm5","Mbp", "Sfrp1","Cd14", "Flt1", "Kdr", "Vwf",
                     "Rgs5","Rpe65")
sdp <- SplitDotPlotGG(mouse_eyes, genes.plot = rev(markers.to.plot),
                      cols.use = c("grey","blue"), x.lab.rot = T, plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")

#===3.3 Identify differential expressed genes by ages =========
Pericytes <- SubsetData(mouse_eyes, ident.use = "Pericytes", subset.raw = T)
Pericytes <- SetAllIdent(Pericytes, id = "conditions")
avg.Pericytes <- log1p(AverageExpression(Pericytes, show.progress = FALSE))
colnames(avg.Pericytes) <- c("aged_mouse","young_mouse")
avg.Pericytes$gene <- rownames(avg.Pericytes)
avg.Pericytes$age <- avg.Pericytes$aged_mouse - avg.Pericytes$young_mouse
head(avg.Pericytes[order(avg.Pericytes$age,decreasing = T),"gene"])

RPE <- SubsetData(mouse_eyes, ident.use = "Retinal Pigment Epithelium", subset.raw = T)
RPE <- SetAllIdent(RPE, id = "conditions")
avg.RPE <- log1p(AverageExpression(RPE, show.progress = FALSE))
colnames(avg.RPE) <- c("aged_mouse","young_mouse")
avg.RPE$gene <- rownames(avg.RPE)
avg.RPE$age <- avg.RPE$aged_mouse - avg.RPE$young_mouse
head(avg.RPE[order(avg.RPE$age,decreasing = T),"gene"])


genes.to.label1 = head(avg.Pericytes[order(avg.Pericytes$age,decreasing = T),"gene"],5)
genes.to.label2 = head(avg.Pericytes[order(avg.Pericytes$age,decreasing = F),"gene"],5)
genes.to.label3 = head(avg.RPE[order(avg.RPE$age,decreasing = T),"gene"],5)
genes.to.label4 = head(avg.RPE[order(avg.RPE$age,decreasing = F),"gene"],5)

p1 <- ggplot(avg.Pericytes, aes(aged_mouse, young_mouse)) + geom_point() + ggtitle("Pericytes")
p1 <- LabelUR(p1, genes = genes.to.label1, avg.Pericytes, 
              adj.u.t = 0.3, adj.u.s = 0.23,text.size = 4)
p1 <- LabelUL(p1, genes = genes.to.label2, avg.Pericytes, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25,text.size = 4)
p2 <- ggplot(avg.RPE, aes(aged_mouse, young_mouse)) + geom_point() + ggtitle("Retinal Pigment Epithelium")
p2 <- LabelUR(p2, genes = genes.to.label3, avg.RPE, 
              adj.u.t = 0.15, adj.u.s = 0.1,text.size = 4)
p2 <- LabelUL(p2, genes = genes.to.label4, avg.RPE, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25,text.size = 4)
plot_grid(p1, p2)

#3.4  Compare DE across all major cell types
#We would need the data for all clusters, as well the subclusters (RPE and hematopoietic cells).
#detect changes in gene expression between 129_B6 and 129_B6_aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to 129_B6_aged or viceversa. 
print("3.4 Compare DE across all major cell types")
mouse_eyes.markers <- FindAllMarkersInSameAge(object = mouse_eyes)

mouse_eyes.gde <- FindAllMarkersbyAge(object = mouse_eyes)
write.csv(x= mouse_eyes.gde, file="./output/mouse_eyes_young_vs_aged.csv")

# 3.5 Compare differential expression between subcluster within all major cell types
# plus visualize all major cell types.
#http://satijalab.org/seurat/de_vignette.html#perform-de-analysis-using-alternative-tests
# Compare subclusters within aged and young
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
# keep the original ident name intact
print("3.5 Compare DE between subcluster within all major cell types, and visualize all major cell types")
# Myeloid Cells===============
Myeloid.cells <- SubsetData(object = mouse_eyes,
                            ident.use = old.ident.ids[(new.cluster.ids %in% "Myeloid Cells")])
Hematopoietic.cells <- SubsetData(object = mouse_eyes,
                                  ident.use = c(old.ident.ids[(new.cluster.ids %in% "Myeloid Cells")],
                                                old.ident.ids[(new.cluster.ids %in% "Lymphoid Cells")]))

table(Myeloid.cells@ident)
TSNEPlotbyAges(Myeloid.cells)

table(Hematopoietic.cells@ident)
TSNEPlotbyAges(Hematopoietic.cells)

Myeloid.cells.markers <- FindAllMarkersInSameAge(Myeloid.cells, write.csv = TRUE)

# Perictyes==========
Pericytes <- SubsetData(object = mouse_eyes,
                        ident.use = old.ident.ids[(new.cluster.ids %in% "Pericytes")])
table(Pericytes@ident)
TSNEPlotbyAges(Pericytes)
Pericytes.markers <- FindAllMarkersInSameAge(Pericytes, write.csv = TRUE)

# Endothelial Cells==========
Endothelial.Cells <- SubsetData(object = mouse_eyes,
                                ident.use = old.ident.ids[(new.cluster.ids %in% "Endothelial Cells")])
table(Endothelial.Cells@ident)
TSNEPlotbyAges(Endothelial.Cells)
Endothelial.Cells.markers <- FindAllMarkersInSameAge(Endothelial.Cells, write.csv = TRUE)



# RPE cells=========
RPE.cells <- SubsetData(object = mouse_eyes,
                        ident.use = old.ident.ids[(new.cluster.ids %in% "Retinal Pigment Epithelium")])
RPE.cells <- FindClusters(object = RPE.cells, 
                          reduction.type = "tsne", 
                          dims.use = 1:2, 
                          resolution = 0.05,
                          save.SNN = TRUE)
table(RPE.cells@ident)
TSNEPlotbyAges(RPE.cells)

RPE.cells.markers <- FindAllMarkersInSameAge(RPE.cells, write.csv = TRUE)

#3.6 Compare subcluster between aged vs young===============
print("3.6 Compare subcluster between aged vs young")
# keep the original ident name intact
# Myeloid Cells===============

Myeloid.gde <- FindAllMarkersbyAge(object = Myeloid.cells)
write.csv(x= Myeloid.gde, file="./output/Myeloid.cells_young_vs_aged.csv")

# Perictyes==========
Pericytes.gde <- FindAllMarkersbyAge(object = Pericytes)
write.csv(x= Pericytes.gde, file="./output/Pericytes_young_vs_aged.csv")

# Endothelial Cells==========
table(Endothelial.Cells@ident)
Endothelial.gde <- FindAllMarkersbyAge(Endothelial.Cells)
write.csv(x= Endothelial.gde, file="./output/Endothelium_young_vs_aged.csv")

# RPE cells=========
RPE.gde <- FindAllMarkersbyAge(RPE.cells)
write.csv(x= RPE.gde, file="./output/RPE_young_vs_aged.csv")
