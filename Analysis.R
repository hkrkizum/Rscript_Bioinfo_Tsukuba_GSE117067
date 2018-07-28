setwd("X:/Bioinfomatics/Class2018/Bioinfo_Tsukuba2018/GSE117067/Rscript")
library(TCC)
# library(Biostrings)
library(biomaRt)
library(gplots)
library(tidyverse)

# Data loarding ------------------------------------------------------------

Rawdata_list <- list.files("Rawdata/",full.names=T, include.dirs = F)
Sample_list <- c("WT_1","WT_2","WT_3","WT_4","KO_1","KO_2","KO_3","KO_4")

Rawdata_list
WT_1 <- readr::read_table2(Rawdata_list[1],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "WT_1"))

WT_2 <- readr::read_table2(Rawdata_list[2],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "WT_2"))
WT_3 <- readr::read_table2(Rawdata_list[3],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "WT_3"))
WT_4 <- readr::read_table2(Rawdata_list[4],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "WT_4"))
KO_1 <- readr::read_table2(Rawdata_list[5],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "KO_1"))
KO_2 <- readr::read_table2(Rawdata_list[6],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "KO_2"))
KO_3 <- readr::read_table2(Rawdata_list[7],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "KO_3"))
KO_4 <- readr::read_table2(Rawdata_list[8],
                           skip = 2, col_names = c("gene_id",
                                                   "Chr",
                                                   "Start",
                                                   "End",
                                                   "Strand",
                                                   "Length",
                                                   "KO_4"))


full_join(WT_1, WT_2) %>% 
    full_join(., WT_3) %>% 
    full_join(., WT_4) %>% 
    full_join(., KO_1) %>% 
    full_join(., KO_2) %>% 
    full_join(., KO_3) %>% 
    full_join(., KO_4) %>% 
    as.data.frame -> rawdata
rownames(rawdata) <- rawdata$gene_id

write.table(rawdata, "rawdata.txt", sep = "\t",row.names = F)

rawdata_paper <- read_delim("preetika_Oct.w_LimmaVoom_DEG_results.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)




# DEG Analysis by TCC  -----------------------------------------------------

DEG_rawdata <- rawdata %>% 
    dplyr::select(c(7:14))                                # select data
Group <- c("WT","WT","WT","WT","KO","KO","KO","KO")       # set Group

## set parameters -------------------------
in_f <- DEG_rawdata                        
out_f <- "df_result.txt"
param_method <- "edger"
param_FDR <- 0.1

## Excution -------------------------------
tcc <- new("TCC", in_f, Group)
tcc <- filterLowCountGenes(tcc)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                       iteration = 3, FDR = 0.05, floorPDEG = 0.05)
Normalized_rawdata <- getNormalizedData(tcc)
tcc <- estimateDE(tcc, test.method = param_method, FDR = param_FDR)
result <- getResult(tcc, sort = TRUE)


result$gene_id <- as.character(result$gene_id)
Normalized_rawdata <- data.frame(Normalized_rawdata) %>% 
    mutate(gene_id = rownames(Normalized_rawdata))
result %>% left_join(., rawdata) %>%
    dplyr::arrange(q.value) %>% 
    dplyr::select(c(1,8:12,2:7)) %>% 
    dplyr::left_join(., Normalized_rawdata) %>% 
    dplyr::select(c(1:6, 13:20, 7:10))-> df_result
write.table(df_result, out_f, sep = "\t", row.names = F)


# df_result %>% 
#     dplyr::filter(q.value < 0.05) %>% dim() 
# 
df_result %>%
    dplyr::filter(q.value < 0.1) %>% dim()
106 * 0.9
# 
# df_result %>% 
#     dplyr::filter(q.value < 0.2) %>% dim() 
# 199 * 0.8
# 
# df_result %>% 
#     dplyr::filter(q.value < 0.3) %>% dim() 
# 335 * 0.7

## make MAplot -----------------------------------------------
png("Result/MAplot.png", height = 3000, width = 3000, res = 300)
plot(tcc)
dev.off()
plot(tcc)

## Get gene list --------------------------------------------
df_result %>%  
    dplyr::filter(q.value < 0.1) %>% 
    dplyr::select(gene_id, m.value) -> DEG_list



# VennDiagram  ------------------------------------------------------------------
library(VennDiagram)

rawdata_paper %>% dplyr::filter(BH.qvalue < 0.25) %>% 
    dplyr::filter(logFC > log10(1.6) | logFC < -log10(1.6)) -> DEG_paper

Mydata <- DEG_list$gene_id
Paper_data <- DEG_paper$id

Venn_data <- list(MyData=Mydata, Paper=Paper_data)
venn.diagram(Venn_data, 
             filename = "Result/Venn_compare.tiff",
             imagetype = "tiff",
             height=3000, width=3000,
             resolution = 500,
             units = "px",
             fill=c(4,7),
             cat.pos=c(0,20),
             cat.dist=c(0.02,0.02),
             cat.cex=c(2,2),
             cex=c(2,2,2)
)

# Clustering Samples --------------------------------------------------------
head(Normalized_rawdata) 
rownames(Normalized_rawdata) <- Normalized_rawdata$gene_id
Cluster <- Normalized_rawdata[,-9]
head(Cluster)
tmp <- clusterSample(Cluster, dist.method = "spearman",
                     hclust.method = "average",
                     unique.pattern = TRUE)
png("Result/sample_cluster.png", height = 3000, width = 3000, res = 300)
plot(tmp)
dev.off()
plot(tmp)

# heat map ---------------------------------------------------------------------
df_result %>%
    dplyr::filter(q.value < 0.1) %>% 
    dplyr::select(gene_id, 
           WT_1, WT_2, WT_3, WT_4,
           KO_1, KO_2, KO_3, KO_4) -> df_heatmap

# Convert Gene ID ----------------------------------------------
db <- useMart("ensembl")
Mg <- useDataset("mmusculus_gene_ensembl", mart = db)

ensid <- df_heatmap$gene_id
res <- getBM(attributes = "external_gene_name",
             filters = "ensembl_gene_id",
             values = ensid,
             mart = Mg)
df_heatmap$gene_symbol <- res$external_gene_name

rownames(df_heatmap) <- df_heatmap$gene_symbol
head(df_heatmap)

df_heatmap <- df_heatmap[,c(-1,-10)]
df_heatmap

heatmap(as.matrix(df_heatmap),Colv = NA)


pdf(file = "Result/Cluster.pdf",
     width = 6, height = 18) # defaults to 7 x 7 inches

lm<-matrix(c(0,2,3,1,4,0),ncol=3)
heatmap.2(as.matrix(df_heatmap),
          Rowv = TRUE,
          Colv = FALSE,
          dendrogram = "none",
          col=greenred(75),
          scale="row",
          key=TRUE,
          keysize=1,
          symkey=TRUE,
          density.info="histogram",
          trace="none",
          cexRow=1,
          margin=c(4,8),
          main="Heat Map 2 (Z score Data)",
          lmat=lm,lwid=c(1,5,2),lhei=c(1,9)
)
dev.off() 


# GO for DAVID --------------------------------------------------------
df_result %>% dplyr::filter(q.value < param_FDR) %>% 
    getBM(attributes = c("ensembl_gene_id","external_gene_name"),
          filters = "ensembl_gene_id",
          values = .,
          mart = Mg) -> tmp
colnames(tmp)[1]  <-　"gene_id"

df_result %>% 
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::filter(m.value > 0) %>% 
    left_join(., tmp) %>% 
    dplyr::select(external_gene_name) -> DEG_list_up

write.table(DEG_list_up, file = "Result/DEG_list_up.txt", 
            row.names = F, col.names = F, quote = F) 

df_result %>% 
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::filter(m.value < 0) %>% 
    left_join(., tmp) %>% 
    dplyr::select(external_gene_name) -> DEG_list_down

write.table(DEG_list_down, file = "Result/DEG_list_down.txt", 
            row.names = F, col.names = F, quote = F) 


# GSEA Gene set enrichment analysis ----------------------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)

MSigDB_c5 <- read.gmt("X:/Bioinfomatics/Data/Annotation/c5.all.v6.2.symbols.gmt")
MSigDB_c2 <- read.gmt("X:/Bioinfomatics/Data/Annotation/c2.all.v6.2.symbols.gmt")
DisOnt <- read_delim("X:/Bioinfomatics/Data/Annotation/all_gene_disease_associations.tsv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
disease2gene <- DisOnt[, c("diseaseId", "geneId")]
disease2name <- DisOnt[, c("diseaseId", "diseaseName")]

# make Gene list --------------------------
df_result %>% dplyr::filter(q.value < param_FDR) %>% 
    getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene"),
          filters = "ensembl_gene_id",
          values = .,
          mart = Mg) -> tmp
colnames(tmp)[1]  <-　"gene_id"
df_result %>% 
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::arrange(desc(m.value)) %>% 
    left_join(., tmp) %>% 
    dplyr::select(gene_id,
                  entrezgene,
                  external_gene_name,
                  m.value) -> d

d %>% dplyr::filter(!is.na(entrezgene)) -> tmp
geneList_GSEA <- tmp$m.value
names(geneList_GSEA) <- as.character(toupper(tmp$entrezgene))
geneList_GSEA = sort(geneList_GSEA, decreasing = TRUE)       # for GSEA

d %>% dplyr::filter(m.value < 0) -> tmp
geneList_OA_down <- 
geneList_OA_down_symbol <- toupper(tmp$external_gene_name)          # for OA down, with Gene symbol
geneList_OA_down_entrezgene <- tmp$entrezgene                # for OA down, with Entrez gene ID


d %>% dplyr::filter(m.value < 0) -> tmp
geneList_OA_down <- toupper(tmp$external_gene_name)          # for OA down





## GSEA -----------------------------------
df_result %>% dplyr::filter(q.value < param_FDR) %>% 
    getBM(attributes = c("ensembl_gene_id","external_gene_name"),
          filters = "ensembl_gene_id",
          values = .,
          mart = Mg) -> tmp
colnames(tmp)[1]  <-　"gene_id"
df_result %>% 
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::arrange(desc(m.value)) %>% 
    left_join(., tmp) %>% 
    dplyr::select(external_gene_name, m.value) -> d
write.table(d[,1], "Result/DEG_list.txt", sep = "\t", row.names = F, quote = F, col.names = F) 

geneList <- d$m.value
names(geneList) <- as.character(toupper(d$external_gene_name))
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

res_GSEA <- GSEA(geneList = geneList,
            exponent = 1,
            nPerm        = 1000,
            minGSSize    = 10,
            maxGSSize    = 10000,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            verbose      = TRUE,
            # by = "DOSE",
            TERM2GENE = MSigDB_c5
            )


## gseGO ----------------------------------
df_result %>% dplyr::filter(q.value < param_FDR) %>% 
    getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene"),
          filters = "ensembl_gene_id",
          values = .,
          mart = Mg) -> tmp
colnames(tmp)[1]  <-　"gene_id"
df_result %>% 
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::arrange(desc(m.value)) %>% 
    left_join(., tmp) %>% 
    dplyr::select(gene_id,entrezgene,external_gene_name, m.value) -> d

d %>% dplyr::filter(!is.na(entrezgene)) -> d

geneList <- d$m.value
names(geneList) <- as.character(toupper(d$entrezgene))
geneList = sort(geneList, decreasing = TRUE)
geneList

res_gseGO <- gseGO(geneList = geneList,
                   ont = "BP",
                   OrgDb = "org.Mm.eg.db",
                   exponent = 1,
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose      = TRUE
                   )

res_gseGO <- gseGO(geneList = geneList,
                   ont = "CC",
                   OrgDb = "org.Mm.eg.db",
                   exponent = 1,
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose      = TRUE
                   )

res_gseGO <- gseGO(geneList = geneList,
                   ont = "MF",
                   OrgDb = "org.Mm.eg.db",
                   exponent = 1,
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose      = TRUE
                   )



## GO over-representation test ---------------------------------------------------------
df_result %>% dplyr::filter(q.value < param_FDR) %>% 
    dplyr::filter(m.value < 0) %>% 
    getBM(attributes = c("ensembl_gene_id","external_gene_name"),
          filters = "ensembl_gene_id",
          values = .,
          mart = Mg) -> tmp
colnames(tmp)[1]  <-　"gene_id"
head(tmp)
dim(tmp)
df_result %>% 
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::arrange(desc(m.value)) %>% 
    right_join(., tmp) %>% 
    dplyr::select(external_gene_name, m.value) -> d
write.table(d[,1], "Result/DEG_list.txt", sep = "\t", row.names = F, quote = F, col.names = F) 
geneList <- toupper(d$external_gene_name)
length(geneList)
head(geneList)

# rawdata %>% 
#     getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene"),
#           filters = "ensembl_gene_id",
#           values = .,
#           mart = Mg) -> tmp
# colnames(tmp)[1]  <-　"gene_id"
# write.table(tmp, "all_gene.txt", sep = "\t", quote = FALSE, col.names = FALSE)


universe <- read_delim("all_gene.txt", "\t",
                       escape_double = FALSE, col_names = FALSE,
                       trim_ws = TRUE)
universe <- toupper(tmp$external_gene_name)
universe    
res_enrich_MSigDB_c5 <- enricher(gene = geneList,
                                 minGSSize = 10,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.2,
                                 TERM2GENE= MSigDB_c5
                                 )
res_enrich_MSigDB_c2 <- enricher(gene = geneList,
                                 # universe = universe,
                                 pvalueCutoff = 0.05,
                                 minGSSize = 10,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.2,
                                 TERM2GENE= MSigDB_c2)

head(summary(res_enrich_MSigDB_c5))
head(summary(res_enrich_MSigDB_c2))
barplot(res_enrich_MSigDB_c5, showCategory=20)
barplot(res_enrich_MSigDB_c2, showCategory=20)


## KEGG over-representation test ---------------------------------------------------------
df_result %>% dplyr::filter(q.value < param_FDR) %>% 
    dplyr::filter(m.value < 0) %>% 
    getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene"),
          filters = "ensembl_gene_id",
          values = .,
          mart = Mg) -> tmp
colnames(tmp)[1]  <-　"gene_id"
df_result %>% 
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::arrange(desc(m.value)) %>% 
    left_join(., tmp) %>% 
    dplyr::select(gene_id,entrezgene,external_gene_name, m.value) -> d

d %>% dplyr::filter(!is.na(entrezgene)) -> d
geneList <- d$entrezgene

res_enrich_KEGG <- enrichKEGG(gene = geneList,
                              organism = 'mmu',
                              pvalueCutoff = 0.05)

head(res_enrich_KEGG)































# 
# # Disease Ontrogy enrichment analysis ------------------------------------------------
# df_result %>% dplyr::filter(q.value < param_FDR) %>% 
#     getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene"),
#           filters = "ensembl_gene_id",
#           values = .,
#           mart = Mg) -> tmp
# colnames(tmp)[1]  <-　"gene_id"
# df_result %>% 
#     dplyr::filter(q.value < param_FDR) %>% 
#     dplyr::arrange(desc(m.value)) %>% 
#     left_join(., tmp) %>% 
#     dplyr::select(gene_id,entrezgene,external_gene_name, m.value) -> d
# 
# d %>% dplyr::filter(!is.na(entrezgene)) -> d
# 
# deg <- as.character(d$entrezgene)
# head(geneList)
# str(geneList)
# read_delim("all_gene.txt", "\t", 
#            escape_double = FALSE, col_names = FALSE, 
#            trim_ws = TRUE) %>%
#     dplyr::filter(!is.na(X4)) -> universe
# 
# deg <- universe$X4[1:1000]
# deg
# res_DisOn <- enricher(gene = deg,
#                       pvalueCutoff = 0.00001,
#                       pAdjustMethod = "BH",
#                       qvalueCutoff = 0.00005,
#                       TERM2GENE = disease2gene,
#                       TERM2NAME = disease2name 
#                       )
# summary(res_DisOn)


