library(Seurat)
library(dplyr)
library(plyr)
library(janitor)
library(pheatmap)
source("./R/Seurat_functions.R")

#====== 4.1 load data  ==========================================
lnames = load(file = "./data/mouse_3eyes_cluster_3.Rda")
lnames
table(mouse_3eyes_cluster_3@ident)
idents <- as.data.frame(table(mouse_3eyes_cluster_3@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Pericytes",
                     "Mesenchymal cells",
                     "Mesenchymal cells",
                     "Mesenchymal cells",
                     "Retinal Pigment Epithelium",
                     "Mesenchymal cells",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Smooth muscle cells",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Smooth muscle cells",
                     "Smooth muscle cells",
                     "Endothelial cells",
                     "Mesenchymal cells",
                     "Endothelial cells",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Myelinating Schwann cells",
                     "Retinal Pigment Epithelium",
                     "Retinal Pigment Epithelium",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Monocytes",
                     "Smooth muscle cells",
                     "CD 14 Monocytes",
                     "Hematopoietic SCs",
                     "Mesenchymal & EC",
                     "T & NK cells",
                     "NK & Monocytes",
                     "Endothelial cells",
                     "Melanocytes",
                     "Myelinating Schwann cells",
                     "NK & Monocytes",
                     "Pericytes & RPE",
                     "EC & RPE",
                     "Endothelial cells")
mouse_3eyes_cluster_3@ident <- plyr::mapvalues(x = mouse_3eyes_cluster_3@ident,
                                               from = old.ident.ids,
                                               to = new.cluster.ids)
TSNEPlot(object = mouse_3eyes_cluster_3, no.legend = F, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("Majore cell types")+
        theme(text = element_text(size=20),     #larger text including legend title							
              legend.position="none",
              plot.title = element_text(hjust = 0.5)) #title in middle

#====== 4.2 split data by conditions  ==========================================
cell.all <- FetchData(mouse_3eyes_cluster_3,"conditions")
groups <- c("young.129_B6", "aged.129_B6","young.B6")
cell.subsets <- lapply(groups, function(x)
        rownames(cell.all)[cell.all$conditions == x])

mouse_3eyes_cluster_3.subsets <- list()
for(i in 1:length(groups)){
        mouse_3eyes_cluster_3.subsets[[i]] <- SubsetData(mouse_3eyes_cluster_3, cells.use =cell.subsets[[i]])
}

table(mouse_3eyes_cluster_3.subsets[[1]]@ident)

p <- list()
for(i in 1:2){
        p[[i]] <- TSNEPlot(object = mouse_3eyes_cluster_3.subsets[[i]],do.label = F, group.by = "ident", 
                           do.return = TRUE, no.legend = TRUE,
                           pt.size = 1,label.size = 4 )+
                ggtitle(groups[i])+
                theme(text = element_text(size=20),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
}
do.call(plot_grid, p)
TSNEPlot(object = mouse_3eyes_cluster_3.subsets[[i]],do.label = F, group.by = "ident", 
         do.return = TRUE, no.legend = F,
         pt.size = 1,label.size = 4 )+
        ggtitle(groups[i])+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle

#====== 4.3 Organize FPKM data ===========================================
FPKM.ori <- readr::read_table2("./data/TDI_FPKM_Quantile_Norm_11262017.txt")
FPKM <- as.data.frame(FPKM.ori)

Major_CellTypes <- sub('\\..*', '', colnames(FPKM)) # remove all num after.
Major_CellTypes <- Major_CellTypes[-1] # remove "gene"
length(Major_CellTypes) # all cell types
length(unique(Major_CellTypes)) # all cell types

# rename rownames with gene
FPKM$gene <- Hmisc::capitalize(tolower(FPKM$gene))
if(length(FPKM$gene) == length(unique(FPKM$gene))){
        rownames(FPKM) <- FPKM$gene
}
FPKM <- FPKM[,-1]

# aggregate all same type of cells' expression
CellTypes <- data.frame(Major_CellTypes = Major_CellTypes)
rownames(CellTypes) <- colnames(FPKM)


FPKM_short <- cbind.data.frame(CellTypes,t(FPKM))
system.time(t_FPKM <- aggregate(. ~ Major_CellTypes, data = FPKM_short, FUN=mean))
rownames(t_FPKM) <- t_FPKM$Major_CellTypes
FPKM_short <- as.data.frame(t(t_FPKM[,-1]))
FPKM_short$genesymbol <- rownames(FPKM_short)
head(FPKM_short[,(ncol(FPKM_short)-3):ncol(FPKM_short)])
FPKM.summary <- FPKM_short
#====== 4.4 Identify Cell Types by Spearman correlation ==================================
Identify_Cell_Types_Spearman <- function(object, gendata,cluster_rows=F,
                                         cluster_cols = F,fontsize_row = 15,
                                         fontsize_col = 15,fontsize =20,
                                         title = ""){
        "
        Calculate Average Expression of each ident of seurat object,
        Calculate spearman correlation between above results and provided gendata dataset
        "
        if(class(object) != "seurat") {
                stop(paste("Error : ", object, " is not a seurat object"))
        }
        if(class(gendata) != "data.frame") {
                stop(paste("Error : ", gendata, " is not a data frame"))
        }
        object.AverageExp <- AverageExpression(object)
        object.AverageExp$genesymbol <- toupper(rownames(object.AverageExp))
        object.Exp <- object.AverageExp[,c(ncol(object.AverageExp),
                                       1:(ncol(object.AverageExp)-1))]# move last column to the first
        print("Merge genes expression:")
        table(gendata$genesymbol %in% object.Exp$genesymbol)
        # merge ===============
        object.Exp_gendata <- inner_join(object.Exp, gendata, by = "genesymbol")
        object.Exp_gendata <- object.Exp_gendata[order(object.Exp_gendata$genesymbol),]
        rownames(object.Exp_gendata) <- object.Exp_gendata$genesymbol
        object.Exp_gendata <- object.Exp_gendata[,-1]
        # Spearman correlation primary ==================
        c <- cor(object.Exp_gendata, method="spearman") # or naive_matrix
        diag(c) <-NA
        ident_num <- length(levels(object@ident))
        object_c_gendata <- c[(ident_num+1):nrow(c),1:ident_num]
        pheatmap(object_c_gendata,cex=.9,
                 cluster_rows= cluster_rows,
                 cluster_cols = cluster_cols,
                 fontsize_row = fontsize_row,
                 fontsize_col = fontsize_col,
                 fontsize = fontsize,
                 main = title)
        actual_cell_types <- apply(object_c_gendata, 2, which.max)
        rename_ident <- data.frame("old.ident.ids" = colnames(object_c_gendata),
                                   "new.cluster.ids" = rownames(object_c_gendata)[actual_cell_types])
        print(rename_ident)
        return(rename_ident)
}

# correlate with FPKM_short
rename_ident <- Identify_Cell_Types_Spearman(object = mouse_3eyes_cluster_3,
                           gendata = FPKM.summary,
                           cluster_rows=T,
                           cluster_cols = T,
                           fontsize_row = 8,
                           fontsize_col = 12,fontsize =20,
          title = "Spearman correlation: mouse_eyes vs pure cell types")
write.csv(rename_ident,file = "./output/rename_ident.csv")

#====== 4.5 Identify Cell Types by Spearman correlation using Mouse Cell Atlas===============================
lnames = load(file = "/Users/yah2014/Dropbox/Public/Olivier/R/scRNAseq-MouseAgedEyes/data/mca.Rda")
lnames
library(cowplot)
p1 <- DimPlot(object = mca, reduction.use = "FItSNE", no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = mca, reduction.use = "FItSNE", no.legend = TRUE, group.by = "Tissue", 
              do.return = TRUE, vector.friendly = TRUE, pt.size = 0.1) + ggtitle("Tissue") + 
        theme(plot.title = element_text(hjust = 0.5))
FeaturePlot(mca, c("S100a9", "Sftpc"), reduction.use = "FItSNE", dark.theme = TRUE, 
            pt.size = 0.1, vector.friendly = TRUE)
plot_grid(p1, p2)
