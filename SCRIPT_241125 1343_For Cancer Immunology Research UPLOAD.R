#############################################---CopyKat---#############################################################################################################
#############################################---CopyKat---#############################################################################################################
#############################################---CopyKat---#############################################################################################################

exp.rawdata <- as.matrix(dat@assays$RNA@counts)

library(devtools)
library(copykat)

copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="MDA3944_Merged", distance="euclidean", norm.cell.names="", n.cores=28,output.seg="FLASE")

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

head(pred.test)

head(CNA.test[ , 1:5])

my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

  heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
  
  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")

dim(pred.test)
dim(exp.rawdata)

copykat_calls <- rep(NA, nrow(dat@meta.data))
names(copykat_calls) <- rownames(dat@meta.data)
copykat_calls[pred.test$cell.names[which(pred.test$copykat.pred=="diploid")]] <- "diploid"
copykat_calls[pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]] <- "aneuploid"

head(copykat_calls)
table(copykat_calls)
sum(is.na(copykat_calls))

copykat_calls <- factor(copykat_calls)
dat@meta.data$copykat_calls <- copykat_calls

saveRDS(dat, "your_path.rds")

rm(list = ls())

#################################---QC---##############################################################################################################################
#################################---QC---##############################################################################################################################
#################################---QC---##############################################################################################################################

dat <- readRDS("your_path.rds")
dat

dat[["percent.MT"]] <- PercentageFeatureSet(dat, pattern = "^MT") ## if human, "MT"

rownames(dat)
range(dat[["percent.MT"]])## !

VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3, group.by = "orig.ident")
plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.MT", group.by = "orig.ident")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot1 + plot2

VlnPlot(dat, "nCount_RNA", group.by = "orig.ident")

mean(dat$nFeature_RNA)
sd(dat$nFeature_RNA)

mean(dat$nCount_RNA)
sd(dat$nCount_RNA)

mean(dat$percent.MT)
sd(dat$percent.MT)

dat <- subset(dat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.MT < 10)

saveRDS(dat, file = "your_path.rds")

#################---Cell Cycle Regression RUN ON CLUSTER---############################################################################################################
#################---Cell Cycle Regression RUN ON CLUSTER---############################################################################################################
#################---Cell Cycle Regression RUN ON CLUSTER---############################################################################################################

s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes 

options(future.globals.maxSize = +Inf)

dat <- NormalizeData(dat) 

dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(dat), 20) 

plot1 <- VariableFeaturePlot(dat) 
plot1 <- LabelPoints(plot = plot1, points = top10, repel = T) 
plot1 

all.genes <- rownames(dat) 
dat <- ScaleData(dat, features = all.genes) 

dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 

head(dat[[]]) 

RidgePlot(dat, features = c("PCNA", "TOP2A", "MCM6", "MKI67", "CDK6"), ncol = 2) 

dat <- RunPCA(dat, features = c(s.genes, g2m.genes)) 
DimPlot(dat) 

#############################################---Harmony---#############################################################################################################
#############################################---Harmony---#############################################################################################################
#############################################---Harmony---#############################################################################################################

dat <- dat %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)

harmony_embeddings <- Embeddings(dat, 'harmony')
harmony_embeddings[1:4, 1:4]

DimHeatmap(dat, dims = 1:30, cells = 500, balanced = TRUE)

dat <- JackStraw(dat, num.replicate = 100)
dat <- ScoreJackStraw(dat, dims = 1:20)
JackStrawPlot(dat, dims = 1:20)

ElbowPlot(dat)

dat7 <- dat %>% 
  RunUMAP(reduction = "harmony", dims = 1:5) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:5) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

DefaultAssay(dat) <- "RNA"

#######################################################USING CYTOTALK##################################################################################################
#######################################################USING CYTOTALK##################################################################################################
#######################################################USING CYTOTALK##################################################################################################

library(reticulate)
use_condaenv("r_reticulate_CytoTalk")
py_config()
conda_remove("r-reticulate")
conda_create(envname = "r_reticulate_CytoTalk", python_version = "3.7.3")  # Create a new Conda environment to facilitate the Python module installation.
use_condaenv("r_reticulate_CytoTalk")
conda_install(envname = "r_reticulate_CytoTalk", "pybind11")  # Install two necessary Python modules for correctly compiling and using the "pcst_fast" Python module.
conda_install(envname = "r_reticulate_CytoTalk", "numpy")
conda_install(envname = "r_reticulate_CytoTalk", "git+https://github.com/fraenkel-lab/pcst_fast.git", pip = TRUE) # To install the "pcst_fast" module.
py_config()
load(file = "scrna_cyto.rda")
lst_scrna <- scrna_cyto

type_a <- "13"
type_b <- "10"

results <- CytoTalk::run_cytotalk(lst_scrna, type_a, type_b, pcg = CytoTalk::pcg_human, lrp = CytoTalk::lrp_human, dir_out = "your_path")
str(results)
results$pathways
results$params

sce <- as.SingleCellExperiment(dat)
sce
lst_scrna <- CytoTalk::from_single_cell_experiment(sce)
dat$cell_type=Idents(dat)
metadata=dat@meta.data
metadata$cell_type
dat$cell_type

index=match(lst_scrna$cell_types, rownames(metadata))
lst_scrna$cell_types=metadata$cell_type[index]
lst_scrna$cell_types
lst_scrna

dat

##################################FOR SUBSETTING BASED ON METADATA
# Renaming based on clusters you defined with similar gene expression
dat <- RenameIdents(object = dat, '9' = '10', '22' = '13', '16' = '13', '5' = '13', '23' = '13', '16' = '13')

# Subset on a value in the object meta data
dat <- subset(x = dat, subset = MVP == "MVP")
NONMVP <- subset(x = dat, subset = MVP == "NONMVP")
dat <- subset(x = original, subset = MVP == "MVP")

# To make life easier
MVP <- NONMVP

sce <- as.SingleCellExperiment(MVP)
sce
lst_scrna <- CytoTalk::from_single_cell_experiment(sce)
MVP$cell_type=Idents(MVP)
metadata=MVP@meta.data
metadata$cell_type
MVP$cell_type

# View(lst_scrna)
index=match(lst_scrna$cell_types, rownames(metadata))
lst_scrna$cell_types=metadata$cell_type[index]
lst_scrna$cell_types
lst_scrna

library(reticulate)
use_condaenv("r_reticulate_CytoTalk")
conda_remove("r-reticulate")

type_a <- "PDGFRB Positive"
type_b <- "25"
results <- CytoTalk::run_cytotalk(lst_scrna, type_a, type_b, pcg = CytoTalk::pcg_human, lrp = CytoTalk::lrp_human, dir_out = "your_path")

#######################biomaRt & clusterProfiler#######################################################################################################################
#######################biomaRt & clusterProfiler#######################################################################################################################
#######################biomaRt & clusterProfiler#######################################################################################################################

library(clusterProfiler)
purrr::simplify
stats::filter
library(biomaRt)
library(tidyverse)

# Prepare the differentially expressed gene list
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$MVP, sep = "_")
immune.combined$GlobalCPSC <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "Candidate perivascular stromal cells_MVP", ident.2 = "Candidate perivascular stromal cells_NONMVP", verbose = FALSE)
head(b.interferon.response, n = 15)
df <- b.interferon.response
write.xlsx(x = b.interferon.response, rowNames = TRUE, file = "your_path.xlsx")
FeaturePlot(immune.combined, features = c("NNMT", "TNFAIP6", "CHI3L1", "COL1A1", "S100A4", "OLFML2B"), split.by = "MVP", max.cutoff = 3, 
            cols = c("grey", "red")) + theme(text = element_text(size = 20))
plots <- VlnPlot(immune.combined, features = c("NNMT", "TNFAIP6", "CHI3L1"), split.by = "MVP", group.by = "GlobalCPSC", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

###########################################VOLCANO PLOT MVP UPREGULATED################################################################################################
###########################################VOLCANO PLOT MVP UPREGULATED################################################################################################
###########################################VOLCANO PLOT MVP UPREGULATED################################################################################################

df <- read.csv("your_path.csv")
View(df)

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

df$diffexpressed <- 'NO'
df$diffexpressed[df$avg_log2FC > 0.6 & df$p_val_adj < 0.05] <- 'UP'
df$diffexpressed[df$avg_log2FC < -0.6 & df$p_val_adj < 0.05] <- 'DOWN'

top30degs <- head(df[order(df$p_val_adj), 'gene'], 15)
df$delabel <- ifelse(df$gene %in% top30degs, df$gene, NA)

#Labeling top 20 and top downregulated
df$delabel <- factor(df$gene, levels = c("your_genes"))

# -log10pvalue looks better_no title or
ggplot(data = df, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = delabel), show.legend = FALSE) +
  geom_vline(xintercept = c(-0.6, 0.6), col = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = c(0.05), col = 'gray', linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  coord_cartesian(ylim = c(0, 300), xlim = c(-2, 4)) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adjusted p-value")) +
  geom_text_repel(max.overlaps = Inf) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20)) +
  theme(text = element_text(size = 20))

library(ggrepel)

#############################Script for GO######################################################################################################
#############################Script for GO######################################################################################################
#############################Script for GO######################################################################################################

filter_df1 <- read.csv("your_path.csv", header=TRUE, dec = ".", fill = TRUE) #min.pct0.25 log2=0.25 only.pos = TRUE

nrow(filter_df1)
annotation = getBM(attributes = c('external_gene_name', 'entrezgene_id', 'description', 'gene_biotype'),
                   filters = ('external_gene_name'),
                   values = filter_df1$gene,
                   mart = ensembl106)

nrow(annotation)

annotated_df = left_join(filter_df1, annotation,
                         by=c('gene'='external_gene_name'))

dim(annotated_df)
dim(filter_df1)
dim(annotation)

anno_df2 = annotated_df[annotated_df$p_val_adj < 0.05,]
ent_gene = getBM(attributes = c('entrezgene_id'),
                 filters = ('external_gene_name'),
                 values = anno_df2$gene,
                 mart = ensembl106)
ent_gene = ent_gene$entrezgene_id
ent_gene

ego = enrichGO(gene = ent_gene,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
)
ego
barplot(ego, showCategory =20)
dotplot(ego, showCategory = 30) + theme(axis.text.y = element_text(size = 12)) + theme(text=element_text(size=21)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "BrBG")))
cnetplot(ego)
goplot(ego)

####################################################HALLMARK Gene Sets###########################################################################################################################
####################################################HALLMARK Gene Sets###########################################################################################################################
####################################################HALLMARK Gene Sets###########################################################################################################################

levels(dat@active.ident)

set.seed(123)  # Set seed for reproducibility
cell_ids <- WhichCells(dat, idents = c("your_cell_populations"))

cell_sample <- sample(cell_ids, size = length(cell_ids) * 0.2)

immune_cells <- subset(x = dat, cells = cell_sample)


library(Seurat)

BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
install.packages("XML")

library(GSEABase)
metabolic_gene_sets <- getGmt("your_path.gmt")
metabolic_list <- lapply(metabolic_gene_sets, geneIds)
names(metabolic_list) <- sapply(metabolic_gene_sets, setName)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

library(AUCell)
immune_rankings <- AUCell_buildRankings(immune_cells@assays$RNA@data, plotStats = TRUE)

immune_metabolic_AUC <- AUCell_calcAUC(metabolic_list, immune_rankings)

set.seed(333)
par(mfrow=c(1,1)) 
cells_assignment <- AUCell_exploreThresholds(immune_metabolic_AUC, plotHist=TRUE, assign=TRUE) 

immune_cells[["metabolic_AUC"]] <- CreateAssayObject(data = getAUC(immune_metabolic_AUC))

# Assuming `immune_cells` is your Seurat object with all relevant cell types
expr_matrix <- as.matrix(immune_cells@assays$RNA@data)  # Use normalized data
# Rank genes for each cell
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = TRUE)
# Calculate AUCell scores for each metabolic pathway in `metabolic_list`
metabolic_AUC <- AUCell_calcAUC(metabolic_list, cells_rankings)
auc_matrix <- as.matrix(getAUC(metabolic_AUC))
colnames(auc_matrix) <- colnames(immune_cells)
# Convert AUC scores to a format compatible with Seurat
immune_cells[["metabolic_scores"]] <- CreateAssayObject(data = auc_matrix)
colnames(auc_matrix)

DefaultAssay(immune_cells) <- "metabolic_scores"
rownames(immune_cells[["metabolic_scores"]])
VlnPlot(immune_cells, features = "HALLMARK-INFLAMMATORY-RESPONSE", group.by = "MSC", assay = "metabolic_scores")

# Calculating p-value
glycolysis_scores <- FetchData(immune_cells, vars = "HALLMARK-INFLAMMATORY-RESPONSE")
immune_cells$Glycolysis_Score <- glycolysis_scores$"HALLMARK-INFLAMMATORY-RESPONSE"
# Kruskal-Wallis test to see if there's a difference in glycolysis scores across cell types
kruskal.test(Glycolysis_Score ~ MSC, data = immune_cells@meta.data)
pairwise_results <- pairwise.wilcox.test(immune_cells$Glycolysis_Score, immune_cells$MSC, p.adjust.method = "BH")
pairwise_results

# ANOVA
# Extract metadata with glycolysis scores and cell types
data <- immune_cells@meta.data
# Perform ANOVA
anova_result <- aov(Glycolysis_Score ~ MSC, data = data)
summary(anova_result)
qqnorm(residuals(anova_result))
qqline(residuals(anova_result))
shapiro.test(residuals(anova_result))
#I don't have this library(car)
#I need the library leveneTest(Glycolysis_Score ~ MSC, data = data)
# Perform Tukey's HSD test for pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
tukey_result

install.packages("openxlsx")
library(openxlsx)
# Perform ANOVA and Tukey HSD test
anova_result <- aov(Glycolysis_Score ~ MSC, data = data)
tukey_result <- TukeyHSD(anova_result)

# Convert the Tukey HSD results to a data frame and include row names
tukey_df <- as.data.frame(tukey_result$MSC)
tukey_df$Comparison <- rownames(tukey_df)  # Add row names as a new column
rownames(tukey_df) <- NULL  # Remove row names from the data frame

# Create a new workbook and add the data
wb <- createWorkbook()
addWorksheet(wb, "Tukey HSD Results")
writeData(wb, "Tukey HSD Results", tukey_df)

# Save the workbook
saveWorkbook(wb, "Tukey_HSD_Results_Myeloid_Inflammatory Response.xlsx", overwrite = TRUE)

#############################################################CellChat########################################################################################################
#############################################################CellChat########################################################################################################
#############################################################CellChat########################################################################################################

library(devtools)
library(CellChat)
library(Seurat)
library(reticulate)

levels(dat@active.ident)
DimPlot(dat, reduction = "umap")
dat <- subset(x = dat, idents = c("your_cell_populations"))
saveRDS(dat, file = "your_path.rds")

data.input <- GetAssayData(dat, assay = "RNA", slot = "data")
# For older object
data.input <- GetAssayData(dat, assay = "RNA", layer = "data")

labels <- Idents(dat)
meta <- data.frame(group = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
saveRDS(cellchat, file = "your_path.rds")

CellChatDB.human <- CellChatDB.human

showDatabaseCategory(CellChatDB.human)

CellChatDB.use <- CellChatDB.human

# Set a subset of cellchatDB for cell-cell communications analysis
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
CellChatDB.use <- subsetDB(CellChatDB.human, search = "ECM-Receptor")
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Cell-Cell Contact")

# Update CellChatDB
interaction <- CellChatDB.human$interaction
complex <- CellChatDB.human$complex
cofactor <- CellChatDB.human$cofactor
geneInfo <- CellChatDB.human$geneInfo

write.csv(interaction, file = "CellChat_Interaction.csv")
write.csv(complex, file = "CellChat_Complex.csv")
write.csv(cofactor, file = "CellChat_Cofactor.csv")
write.csv(geneInfo, file = "CellChat_geneInfo.csv")

CellChatDB.human.updated <- list()
CellChatDB.human.updated$interaction <- interaction
CellChatDB.human.updated$complex <- complex
CellChatDB.human.updated$cofactor <- cofactor
CellChatDB.human.updated$geneInfo <- geneInfo

CellChatDB.use <- CellChatDB.human.updated

library(patchwork)
library(circlize)
options(stringsAsFactors = FALSE)

CellChatDB.use <- CellChatDB.human

cellchat@DB <- CellChatDB.use

# Subset data that is only in the cellchatDB
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Optional. The function projectData in CellChat is used to map inferred cell-cell communication data from one species to another by projecting data using a known
# protein-protein interaction (PPI) network, such as PPI.human.
cellchat <- projectData(cellchat, PPI.human)

# Only sse the projected data to compute communication probabilities
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")

dev.off

# Examine signaling sent from each cell type
mat <- cellchat@net$weight
par(mfrow = c(3,5), xpd = TRUE)
for (i in 1:nrow(mat)) {mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <-mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])}

saveRDS(cellchat, file = "your_path.rds")


cellchat@netP[["pathways"]]

extractEnrichedLR(cellchat, signaling = c(cellchat@netP[["pathways"]]),
                  geneLR.return = TRUE)

netAnalysis_contribution(cellchat,
                         signaling = c(cellchat@netP[["pathways"]]),
                         title = "Contribution of each LR pairs")

netAnalysis_contribution(cellchat,
                         signaling = c(cellchat@netP[["pathways"]][1:3]),
                         title = "Contribution of Ligand-Receptor Pairs for Top 3 Pathways")

extractEnrichedLR(cellchat, signaling = "TNF", geneLR.return = FALSE)
netAnalysis_contribution(cellchat, signaling = "TNF")

# You can plot any signaling pathways from those in your dataset
extractEnrichedLR(cellchat, signaling = "CCL", geneLR.return = FALSE)
netAnalysis_contribution(cellchat, signaling = "CCL")

extractEnrichedLR(cellchat, signaling = "COLLAGEN", geneLR.return = FALSE)
netAnalysis_contribution(cellchat, signaling = "COLLAGEN")

# Chord diagram
par(mfrow = c(1, 1), xpd=TRUE)
par(cex = 0.5)
netVisual_aggregate(cellchat, signaling = "TNF", layout = "chord")
netVisual_chord_cell(cellchat, signaling = "TNF")
netVisual_chord_gene(cellchat, signaling = "TNF")
netVisual_chord_gene(cellchat, signaling = "PDGF")

# Chord diagram: group cell clusters
levels(cellchat@idents)

group.cellType <- c(rep("",4), rep("DC", 4),rep("TC",4))
names(group.cellType)<- levels(cellchat@idents)
par(mfrow = c(1,1),xpd = TRUE)
par(cex = 0.5)
netVisual_chord_cell(cellchat,signaling = "TNF",
                     group= group.cellType,
                     title.name = paste0("TNF_", "signaling network"))

# Chord diagram: define source and target cell types
levels(cellchat@idents)
netVisual_chord_gene(cellchat, sources.use = c("your_cell_populations"),
                     signaling = "VEGF", legend.pos.y = 8)

#############################################################Crossreferencing HALLMARKwith CellChat########################################################################################################
#############################################################Crossreferencing HALLMARKwith CellChat########################################################################################################
#############################################################Crossreferencing HALLMARKwith CellChat########################################################################################################

library(GSEABase)

# Load the Hallmark glycolysis gene set
hallmark_glycolysis <- getGmt("your_path.gmt")

# Extract the list of glycolysis genes
glycolysis_genes <- geneIds(hallmark_glycolysis)[[1]]  # Assuming only one gene set in the .gmt

# Assuming `cellchat` is your CellChat object
library(CellChat)

# Extract interaction data from the CellChat object
interaction_data <- cellchat@DB$interaction

# Split receptor complexes into individual genes
interaction_data$receptor_split <- strsplit(interaction_data$receptor, "_")

# Combine ligand and individual receptor genes into a single column
interaction_data$all_genes <- mapply(function(ligand, receptors) {
  unique(c(ligand, receptors))
}, interaction_data$ligand, interaction_data$receptor_split)

# Filter interactions with glycolysis genes
glycolysis_interactions <- interaction_data[sapply(interaction_data$all_genes, function(genes) {
  any(genes %in% glycolysis_genes)
}), ]

# View the interactions
print(glycolysis_interactions)

# Count the number of interactions involving glycolysis genes per pathway
pathway_summary <- glycolysis_interactions %>%
  group_by(pathway_name) %>%
  summarise(
    InteractionCount = n(),
    GenesInvolved = unique(unlist(all_genes[all_genes %in% glycolysis_genes]))
  )

# View pathway-level summary
print(pathway_summary)

library(ggplot2)

ggplot(pathway_summary, aes(x = reorder(pathway_name, InteractionCount), y = InteractionCount)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Upregulated CellChat Pathways with HALLMARK Glycolysis Gene Set Gene Overlap",
       x = "Pathway",
       y = "Number of Interactions")

write.csv(glycolysis_interactions, "your_path.csv", row.names = FALSE)

# Get relevant ligands and receptors from glycolysis interactions
glycolysis_genes_involved <- unique(unlist(glycolysis_interactions$all_genes))

# Extract the expression matrix
expression_matrix <- cellchat@data.signaling

# Check the structure of the matrix
str(expression_matrix)

# Confirm overlap of glycolysis-related genes with the expression matrix
glycolysis_genes_involved <- unique(unlist(glycolysis_interactions$all_genes))
overlapping_genes <- intersect(rownames(expression_matrix), glycolysis_genes_involved)
print(overlapping_genes)

# Subset the expression matrix for overlapping glycolysis genes
expression_subset <- expression_matrix[overlapping_genes, , drop = FALSE]

# Inspect the subset
str(expression_subset)

# Extract metadata from CellChat
cell_metadata <- cellchat@meta
cell_metadata$idents <- cell_metadata$group  # Assuming 'group' corresponds to cell types

# Align metadata with columns of the expression matrix
cell_metadata <- cell_metadata[colnames(expression_subset), ]

library(dplyr)

# Combine expression data with metadata
expression_data <- as.data.frame(t(expression_subset))
expression_data$CellType <- cell_metadata$idents

# Summarize by cell type
cell_type_summary <- expression_data %>%
  group_by(CellType) %>%
  summarise(across(everything(), mean))

print(cell_type_summary)

colnames(cell_type_summary)

ggplot(cell_type_summary, aes(x = reorder(CellType, -BMP8A), y = BMP8A)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Average Expression of BMP8A by Cell Type",
       x = "Cell Type",
       y = "Expression of BMP8A")

ggplot(cell_type_summary, aes(x = reorder(CellType, -MIF), y = MIF)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Average Expression of MIF by Cell Type",
       x = "Cell Type",
       y = "Expression of MIF")

library(dplyr)

# Add a column for total expression
cell_type_summary <- cell_type_summary %>%
  mutate(TotalExpression = rowSums(select(., -CellType)))

# Plot total expression
ggplot(cell_type_summary, aes(x = reorder(CellType, -TotalExpression), y = TotalExpression)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Total Expression of All Fatty Acid Metabolism Genes by Cell Type",
       x = "Cell Type",
       y = "Total Expression")

library(pheatmap)

# Convert to matrix format
heatmap_matrix <- as.matrix(cell_type_summary[-1])  # Remove CellType column
rownames(heatmap_matrix) <- cell_type_summary$CellType

# Create heatmap
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,  # Cluster cell types
         cluster_cols = TRUE,  # Cluster genes
         main = "Expression of Fatty Acid Metabolism Genes across Cell Types")

write.csv(cell_type_summary, "your_path.csv", row.names = FALSE)
