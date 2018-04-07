library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/mouse_3eyes_alignment.Rda")
lnames

# all marker genes
Featureplot <- function(x){
        p <- FeaturePlot(object = mouse_3eyes, 
#                         reduction.use = "FItSNE",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), pt.size = 0.5)
        return(p)
}
Adipocytes <- MouseGenes(mouse_3eyes,c("SLC36A2","P2RX5","MYF5","UCP1","TRIP4","ASCC1"))
Endothelium <- MouseGenes(mouse_3eyes,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                        "Vwf","EMCN","Car4"))
Epithelium <- MouseGenes(mouse_3eyes,c("KRT19","Epcam","KRT5",
                         "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                         "Rlbp1","Msln","Upk3b","Lrrn4"))
RPE <- MouseGenes(mouse_3eyes,c("Rpe65","Rlbp1"))
Fibroblast <- MouseGenes(mouse_3eyes,c("FGF1","FGF9","SFRP1"))
Hematopoietic <- MouseGenes(mouse_3eyes,c("PTPRC","LAPTM5","SRGN"))
Myeloid <-  MouseGenes(mouse_3eyes,c("PPBP","GNG11","HBA2","HBB","Cma1","Mcpt4","Tpsb2",
                                     "Cpa3","LYZ","S100A9","CD14","CCL2","FCGR3A","MS4A7","VMO1"))
Lymphoid <- MouseGenes(mouse_3eyes,c("CD3G","CD3D","CD2","Cd19","CD79A","MS4A1",
                                     "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
Melanocytes <- MouseGenes(mouse_3eyes,c("Pmel","Mlana"))
Mesenchymal <- MouseGenes(mouse_3eyes,c("Pdgfrb","Vim","Has2","Dcn"))
Myelinating_Schwann_cells <- MouseGenes(mouse_3eyes,c("MBP","MPZ"))
Pericytes <- MouseGenes(mouse_3eyes,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                      "Myh11","Mylk","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- MouseGenes(mouse_3eyes,c("Acta2","Myh11"))
Stem_cell <- MouseGenes(mouse_3eyes,c("POU5F1","FUT4","CD34","PROM1","ABCG2","Runx1","ATXN1",
                                      "Nes","NCAM","NGFR"))
Stromal_fibroblasts <- MouseGenes(mouse_3eyes,c("DCN","COL6A1","TIMP3","PDGFRA"))
# Featureplot
Featureplot(Adipocytes) # Adipocytes
Featureplot(Endothelium) # Endothelial Cells
Featureplot(Epithelium) # Epithelium
Featureplot(RPE) # RPE
Featureplot(Fibroblast) # Fibroblasts
Featureplot(Hematopoietic) # Hematopoietic cells
Featureplot(Myeloid) # Myeloid cells
Featureplot(Lymphoid) # Lymphoid cells
Featureplot(Melanocytes) # Melanocytes
Featureplot(Mesenchymal) # Mesenchymal cells
Featureplot(Myelinating_Schwann_cells) # Myelinating Schwann cells
Featureplot(Pericytes) # Pericytes
Featureplot(Smooth_muscle_cells)
Featureplot(Stem_cell)
Featureplot(Stromal_fibroblasts)
# The SplitDotPlotGG function can be useful for viewing conserved cell type markers
# across conditions, showing both the expression level and the percentage of cells
# in a cluster expressing any given gene. 
# Here we plot 1-3 strong marker genes for each of our 13 clusters.
markers.to.plot <- c(Myelinating_Schwann_cells,Melanocytes,Hematopoietic[1:2], Lymphoid[1:2],
                     Myeloid[c(8,10)],Endothelium[c(1:3,5,7)], 
                     Stem_cell[c(5,3)],RPE,Mesenchymal[4],Pericytes[c(1,4,6:9)],
                     Smooth_muscle_cells)
markers.to.plot <- unique(markers.to.plot)
DotPlot(mouse_3eyes, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
                      dot.scale = 8, do.return = T)

# Rename ident for DotPlot
table(mouse_3eyes@ident)
idents <- as.data.frame(table(mouse_3eyes@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Mesenchymal cells 0",
                     "RPE 1",
                     "Pericytes 2",
                     "Smooth muscle cells 3",
                     "Endothelial cells 4",
                     "Endothelial cells 5",
                     "RPE 6",
                     "Endothelial cells 7",
                     "Mesenchymal cells 8",
                     "Endothelial cells 9",
                     "Schwann cells 10",
                     "Monocytes 11",
                     "Monocytes 12",
                     "T cells 13",
                     "Melanocytes 14",
                     "Schwann cells 15")

mouse_3eyes@ident <- plyr::mapvalues(x = mouse_3eyes@ident,
                               from = old.ident.ids,
                               to = new.cluster.ids)
DotPlot(mouse_3eyes, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
#===================================================
# mouse_3eyes <- RenameIdentBack(mouse_3eyes)
# How many cells are in each cluster
lnames = load(file = "./data/mouse_3eyes_alignment.Rda")
lnames
table(mouse_3eyes@ident)
idents <- as.data.frame(table(mouse_3eyes@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Mesenchymal cells",
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
                                     from = old.ident.ids,
                                     to = new.cluster.ids)
table(mouse_3eyes@ident)
TSNEPlot(object = mouse_3eyes, no.legend = F, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("Majore cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              legend.position="none",
              plot.title = element_text(hjust = 0.5)) #title in middle

#=====2.2 - A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
freq_table <- prop.table(x = table(mouse_3eyes@ident, mouse_3eyes@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table


#====== 2.3 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
# Compare clusters for each dataset
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
TSNEPlot(object = mouse_3eyes.subsets[[i]],do.label = F, group.by = "ident", 
         do.return = TRUE, no.legend = F,
         pt.size = 1,label.size = 4 )+
        ggtitle(groups[i])+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

#=== 2.4 Further subdivisions with FIt-SNE  ===========
# Rename ident for DotPlot
lnames = load(file = "./data/mouse_3eyes_FIt-SNE.Rda")
lnames
Featureplot <- function(x){
        p <- FeaturePlot(object = mouse_3eyes_FItSNE, 
                         reduction.use = "FItSNE",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), pt.size = 0.5)
        return(p)
}
markers.to.plot <- c(Myelinating_Schwann_cells,Melanocytes,Hematopoietic[1:2], Lymphoid[c(1:2,8)],
                     Myeloid[c(8,10)],Endothelium[c(1:3,5,7)], 
                     Stem_cell[c(5,3)],RPE,Mesenchymal[4],Pericytes[c(1,4,6:9)],
                     Smooth_muscle_cells)
markers.to.plot <- unique(markers.to.plot)
DotPlot(mouse_3eyes_FItSNE, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

table(mouse_3eyes_FItSNE@ident)
idents <- as.data.frame(table(mouse_3eyes_FItSNE@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Pericytes 0",
                     "Mesenchymal cells 1",
                     "RPE 2",
                     "Mesenchymal cells 3",
                     "Mesenchymal cells 4",
                     "Mesenchymal cells 5",
                     "Endothelial cells 6",
                     "RPE 7",
                     "Smooth muscle cells 8",
                     "Endothelial cells 9",
                     "RPE 10",
                     "Endothelial cells 11",
                     "Endothelial cells 12",
                     "Smooth muscle cells 13",
                     "Smooth muscle cells 14",
                     "Mesenchymal cells 15",
                     "Endothelial cells 16",
                     "Schwann cells 17",
                     "Endothelial cells 18",
                     "RPE 19",
                     "Endothelial cells 20",
                     "RPE 21",
                     "Monocytes 22",
                     "RPE 23",
                     "CD 14 Monocytes 24",
                     "Endothelial cells 25",
                     "Hematopoietic SCs 26",
                     "T & NK cells 27",
                     "RPE 28",
                     "Smooth muscle cells 29",
                     "Endothelial cells 30",
                     "Monocytes & NK 31",
                     "Melanocytes 32",
                     "Monocytes & NK 33",
                     "Schwann cells 34")

mouse_3eyes_FItSNE@ident <- plyr::mapvalues(x = mouse_3eyes_FItSNE@ident,
                                     from = old.ident.ids,
                                     to = new.cluster.ids)
DotPlot(mouse_3eyes_FItSNE, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
#===================================================
# mouse_3eyes <- RenameIdentBack(mouse_3eyes)
# How many cells are in each cluster
lnames = load(file = "./data/mouse_3eyes_FIt-SNE.Rda")
lnames
table(mouse_3eyes_FItSNE@ident)
idents <- as.data.frame(table(mouse_3eyes_FItSNE@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Pericytes",
                     "Mesenchymal cells",
                     "Retinal Pigment Epithelium",
                     "Mesenchymal cells",
                     "Mesenchymal cells",
                     "Mesenchymal cells",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Smooth\n muscle cells",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Smooth\n muscle cells",
                     "Smooth\n muscle cells",
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
                     "Smooth\n muscle cells",
                     "Endothelial cells",
                     "Monocytes & NK cells",
                     "Melanocytes",
                     "Monocytes & NK cells",
                     "Schwann cells")
mouse_3eyes_FItSNE@ident <- plyr::mapvalues(x = mouse_3eyes_FItSNE@ident,
                                     from = old.ident.ids,
                                     to = new.cluster.ids)
table(mouse_3eyes_FItSNE@ident)
DimPlot(object = mouse_3eyes_FItSNE, reduction.use = "FItSNE", 
        no.legend = TRUE, do.return = TRUE, do.label = TRUE,
        vector.friendly = FALSE, pt.size = 1,label.size = 5 ) + 
        ggtitle("FItSNE plot of all clusters") + 
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
