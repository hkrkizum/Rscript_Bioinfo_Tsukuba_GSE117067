setwd("X:/Bioinfomatics/Class2018/Bioinfo_Tsukuba2018/GSE117067/Rscript_Bioinfo_Tsukuba_GSE117067")
setwd("C:/Bioinfomatics/Class2018/Bioinfo_Tsukuba/Rscript")
library(TCC)
library(Biostrings)
library(biomaRt)
library(gplots)
library(VennDiagram)
library(clusterProfiler)
library(org.Mm.eg.db)
library(MeSH.db)
library(meshr)
library(fdrtool)
library(MeSH.Mmu.eg.db)
library(tidyverse)

## my function ---------------------------
grep_chr_1 <- function(x, y){   #ハイフンから前
    a <- str_locate(x, y)
    b <- substr(x, 1, a-1)
    return(b)
}

grep_chr_2 <- function(x, y){   #ハイフンから後
    a <- str_locate(x, y)
    b <- substr(x, a+1, nchar(x))
    return(b)
}

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

write.table(rawdata, "Result/rawdata.txt", sep = "\t",row.names = F)

rawdata_paper <- read_delim("preetika_Oct.w_LimmaVoom_DEG_results.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

head(rawdata)
head(rawdata_paper)

rawdata_paper %>%  # Check number of DEG in paper 
    dplyr::filter(BH.qvalue < 0.25) %>%
    arrange(desc(logFC)) %>% 
    dim()

# DEG Analysis by TCC  -----------------------------------------------------

DEG_rawdata <- rawdata %>% 
    dplyr::select(c(7:14))                                # select data
Group <- c("WT","WT","WT","WT","KO","KO","KO","KO")       # set Group

## set parameters -------------------------
in_f <- DEG_rawdata                        
out_f <- "Result/df_result.txt"
param_method <- "edger"
param_FDR <- 0.25

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

df_result %>%
  dplyr::select(-(2:6)) %>% 
  dplyr::filter(q.value < param_FDR) -> tmp
write.table(tmp, "Result/df_result_2.txt",sep = "\t", quote = F, row.names = F)

df_result %>%    # Check DEG list
    dplyr::filter(q.value < 0.25) %>% 
    dim()

## make MAplot -----------------------------------------------
out_f2 <- "Result/MAplot_rev2.png"                  #出力ファイル名を指定してout_f2に格納
param_fig <- c(3000, 3000, 300)               #ファイル出力時の横幅と縦幅を指定(単位はピクセル)
param_mar <- c(4, 4, 0, 0)             #下、左、上、右の順で余白を指定(単位は行)
param_col <- "magenta"                 #色を指定
param_cex <- 3                       #点の大きさを指定(2なら通常の2倍、0.5なら通常の0.5倍)
param_pch <- 15                        #点の形状を指定(詳細はこちらとか)

png(out_f2, 
    pointsize=13,
    width=param_fig[1],
    height=param_fig[2],
    res = param_fig[3])#出力ファイルの各種パラメータを指定
par(mar=param_mar)                     #余白を指定
plot(tcc, 
     FDR=param_FDR, 
     xlim=c(-2, 20), 
     ylim=c(-4.5, 4.5),
     cex=1.2, cex.lab=1.2,                #param_FDRで指定した閾値を満たすDEGをマゼンタ色にして描画
     cex.axis=1.2, main="",               #param_FDRで指定した閾値を満たすDEGをマゼンタ色にして描画
     xlab="A = (log2(G2) + log2(G1))/2",  #param_FDRで指定した閾値を満たすDEGをマゼンタ色にして描画
     ylab="M = log2(G2) - log2(G1)" )
legend("topright", 
       c(paste("DEG(FDR<", param_FDR, ")", sep=""),"non-DEG"),
       col=c("magenta", "black"),
       pch=20, cex=1.2)                   #凡例を作成
dev.off()     


# png("Result/MAplot.png", height = 3000, width = 3000, res = 300)
# plot(tcc)
# dev.off()
plot(tcc)

## Get gene list --------------------------------------------
df_result %>%  
    dplyr::filter(q.value < param_FDR) %>% 
    dplyr::select(gene_id, m.value) -> DEG_list

# VennDiagram  ------------------------------------------------------------------
rawdata_paper %>% dplyr::filter(BH.qvalue < 0.25) -> DEG_paper

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
             cat.pos=c(0, 25),
             cat.dist=c(0.02,0.05),
             cat.cex=c(2,2),
             cex=c(2,2,2)
)

## Schatter plot My data vs Paper -----------------------------------------------------------
colnames(rawdata_paper)[2:9] <- c("WT_1_p", "WT_2_p","WT_3_p","WT_4_p",
                                  "KO_1_p", "KO_2_p", "KO_3_p", "KO_4_p")
rawdata_paper$gene_id <- rawdata_paper$id
rawdata_paper
Normalized_rawdata %>% 
  inner_join(., rawdata_paper) %>%
  dplyr::select("gene_id", "WT_1","WT_2","WT_3","WT_4","KO_1","KO_2","KO_3","KO_4",
                "WT_1_p", "WT_2_p","WT_3_p","WT_4_p","KO_1_p", "KO_2_p", "KO_3_p", "KO_4_p") %>% 
  dplyr::mutate(mean_WT = (WT_1 + WT_2 + WT_3 + WT_4)/4) %>% 
  dplyr::mutate(mean_KO = (KO_1 + KO_2 + KO_3 + KO_4)/4) %>% 
  dplyr::mutate(mean_WT_p = (WT_1_p + WT_2_p + WT_3_p + WT_4_p)/4) %>% 
  dplyr::mutate(mean_KO_p = (KO_1_p + KO_2_p + KO_3_p + KO_4_p)/4) %>% 
  dplyr::select(gene_id, mean_WT,mean_KO, mean_WT_p, mean_KO_p) %>% 
  dplyr::arrange(desc(mean_WT)) -> df_schatter
df_schatter[,-1] <- round(df_schatter[,-1], 2)
head(df_schatter)

res <- lm(df_schatter[1000:dim(df_schatter)[1],]$mean_WT ~ df_schatter[1000:dim(df_schatter)[1],]$mean_WT_p)
summary(res)
res$coefficients
res <- lm(df_schatter[1000:dim(df_schatter)[1],]$mean_KO ~ df_schatter[1000:dim(df_schatter)[1],]$mean_KO_p)
summary(res)
res$coefficients

g <-ggplot(df_schatter, aes(x=mean_WT_p, y=mean_WT_p)) +
  geom_point() +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlab("FPRM in paper data") +
  ylab("FPRM in my analysis") +
  theme(panel.background=element_rect(fill="white"),
        axis.line=element_line(colour="black"),
        axis.title=element_text(size=10),
        axis.text=element_text(size=10)
  )
g
ggsave(filename = "Result/Schatter_plot_WT.tiff", plot = g,
       width = 128, height = 128, 
       units = ("mm"), dpi = 300)


g <-ggplot(df_schatter, aes(x=mean_KO_p, y=mean_KO)) +
  geom_point() +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  xlab("FPRM in paper data") +
  ylab("FPRM in my analysis") +
  theme(panel.background=element_rect(fill="white"),
        axis.line=element_line(colour="black"),
        axis.title=element_text(size=10),
        axis.text=element_text(size=10)
  )
g
ggsave(filename = "Result/Schatter_plot_KO.tiff", plot = g,
       width = 128, height = 128, 
       units = ("mm"), dpi = 300)

g <-ggplot(df_schatter[10:dim(df_schatter)[1],], aes(x=mean_WT_p, y=mean_WT_p)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_hline(yintercept=0) + geom_vline(xintercept=0) +
    xlab("FPRM in paper data") +
    ylab("FPRM in my analysis") +
    theme(panel.background=element_rect(fill="white"),
          axis.line=element_line(colour="black"),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10)
    )
g
ggsave(filename = "Result/Schatter_plot_WT_cut10.tiff", plot = g,
       width = 128, height = 128, 
       units = ("mm"), dpi = 300)


g <-ggplot(df_schatter[10:dim(df_schatter)[1],], aes(x=mean_KO_p, y=mean_KO)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_hline(yintercept=0) + geom_vline(xintercept=0) +
    xlab("FPRM in paper data") +
    ylab("FPRM in my analysis") +
    theme(panel.background=element_rect(fill="white"),
          axis.line=element_line(colour="black"),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10)
    )
g
ggsave(filename = "Result/Schatter_plot_KO_cut10.tiff", plot = g,
       width = 128, height = 128, 
       units = ("mm"), dpi = 300)



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
    dplyr::filter(q.value < param_FDR) %>% 
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
     width = 6, height = 25) # defaults to 7 x 7 inches

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
          lmat=lm,lwid=c(1,4,2),lhei=c(1,9)
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

DAVIT_Down_BP <- read_delim("Webresult/DAVIT_Down_BP.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

tmp1 <- grep_chr_1(DAVIT_Down_BP$Term, "~") 
tmp2 <- grep_chr_2(DAVIT_Down_BP$Term, "~") 
DAVIT_Down_BP$Term_id <- tmp1
DAVIT_Down_BP$Term_name <- tmp2

DAVIT_Down_BP <- DAVIT_Down_BP %>% 
    dplyr::select(Term_id, Term_name, 4:13) 
write.table(DAVIT_Down_BP, "Result/DAVIT_Down_BP.txt", sep = "\t",
            quote = F, row.names = F)


# Enrichment analysis ----------------------------------------------------------------------------------------------
MSigDB_c5 <- read.gmt("X:/Bioinfomatics/Data/Annotation/c5.all.v6.2.symbols.gmt")
MSigDB_c2 <- read.gmt("X:/Bioinfomatics/Data/Annotation/c2.all.v6.2.symbols.gmt")
DisOnt <- read_delim("X:/Bioinfomatics/Data/Annotation/all_gene_disease_associations.tsv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
disease2gene <- DisOnt[, c("diseaseId", "geneId")]
disease2name <- DisOnt[, c("diseaseId", "diseaseName")]

head(MSigDB_c5)
# make Backgrounf Gene list -----------------------------------------------------------------------------------

# df_result %>% 
#     getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene"),
#           filters = "ensembl_gene_id",
#           values = .$gene_id,
#           mart = Mg) -> tmp
# write.table(tmp, "Result/df_result_allgene.txt", 
#             sep = "\t", quote = FALSE, col.names = T, row.names = F)

df_result_allgene <- read_delim("Result/df_result_allgene.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(df_result_allgene)[1]  <-　"gene_id"
df_result %>% left_join(., df_result_allgene) -> df_result_allgene
universe_symbol <- toupper(df_result_allgene$external_gene_name)
universe_Ensambl <- toupper(df_result_allgene$gene_id)

df_result_allgene %>% dplyr::select(entrezgene, m.value) -> tmp
tmp <- na.omit(tmp)
tmp[!duplicated(tmp$entrezgene),] -> tmp
universe_entranz <- tmp$entrezgene

# for GSEA ----------------------------------------------------------------------------------------------------
df_result_allgene <- read_delim("Result/df_result_allgene.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(df_result_allgene)[1]  <-　"gene_id"
df_result %>% left_join(., df_result_allgene) -> df_result_allgene

df_result_allgene %>% dplyr::select(gene_id, external_gene_name, m.value) -> tmp
tmp[!duplicated(tmp$gene_id),] -> tmp
tmp[!duplicated(tmp$external_gene_name),] -> tmp
tmp[duplicated(tmp$external_gene_name),1]

geneList_GSEA <- tmp$m.value
names(geneList_GSEA) <- as.character(toupper(tmp$external_gene_name))
geneList_GSEA <- sort(geneList_GSEA, decreasing = TRUE)       

df_result_allgene %>% dplyr::select(entrezgene, m.value) -> tmp
tmp <- na.omit(tmp)
tmp[!duplicated(tmp$entrezgene),] -> tmp
geneList_gseGO <- tmp$m.value
names(geneList_gseGO) <- as.character(tmp$entrezgene)
geneList_gseGO <- sort(geneList_gseGO, decreasing = TRUE)       

# for Over expression Analysis  -----------------------------------------------------------------------------
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

d %>% dplyr::filter(m.value < 0) -> tmp
geneList_OA_down_Ensemble   <- tmp$gene_id                      # for OA down, with Ensemble
geneList_OA_down_symbol     <- toupper(tmp$external_gene_name)  # for OA down, with Gene symbol

tmp %>% 
    dplyr::select(entrezgene) %>%
    na.omit() %>% 
    .[!duplicated(.$entrezgene),] -> geneList_OA_down_entrezgene # for OA down, with Entrez gene ID

d %>% dplyr::filter(m.value > 0) -> tmp
geneList_OA_up_Ensemble   <- tmp$gene_id                        # for OA up, with Ensemble
geneList_OA_up_symbol     <- toupper(tmp$external_gene_name)    # for OA up, with Gene symbol

tmp %>% 
    dplyr::select(entrezgene) %>%
    na.omit() %>% 
    .[!duplicated(.$entrezgene),] -> geneList_OA_up_entrezgene # for OA up, with Entrez gene ID

## GSEA -----------------------------------------------------------------------------------------
res_GSEA <- GSEA(geneList = geneList_GSEA,
            exponent = 1,
            nPerm        = 1000,
            minGSSize    = 10,
            maxGSSize    = 500,
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            verbose      = TRUE,
            by = "fgsea",
            TERM2GENE = MSigDB_c5
            )
rownames(res_GSEA@result) <- 1:dim(res_GSEA@result)[1]
res_GSEA@result$ID <- grep_chr_2(res_GSEA@result$ID, "_") %>% 
    gsub("_", " ", .) %>% 
    tolower() -> res_GSEA@result$ID

res_GSEA@result[res_GSEA@result$core_enrichment %>% grep("DLG1",.),] %>% head()

gseaplot(res_GSEA, res_GSEA@result$ID[1])

res_GSEA@result %>% 
    dplyr::select(-Description) %>% 
    write.table(., "Result/GSEA.txt", sep = "\t",row.names = F,quote = F)
dim(res_GSEA@result)

## gseGO ----------------------------------
res_gseGO <- gseGO(geneList = geneList_gseGO,
                   ont = "BP",
                   OrgDb = "org.Mm.eg.db",
                   exponent = 1,
                   nPerm        = 1000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   verbose      = TRUE
                   )
res_gseGO@result$Description
png("Result/GSEAplot.png", height = 2000, width = 3000, res = 300)
gseaplot(res_gseGO, res_gseGO@result$ID[2])
dev.off()
write.table(res_gseGO@result, "Result/gesGO.txt", sep = "\t",row.names = F,quote = F)


res_gseGO_CC <- gseGO(geneList = geneList_gseGO,
                   ont = "CC",
                   OrgDb = "org.Mm.eg.db",
                   exponent = 1,
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   verbose      = TRUE
                   )

res_gseGO_MF <- gseGO(geneList = geneList_gseGO,
                   ont = "MF",
                   OrgDb = "org.Mm.eg.db",
                   exponent = 1,
                   nPerm        = 10000,
                   minGSSize    = 10,
                   maxGSSize    = 1000,
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   verbose      = TRUE
                   )

## GO over-representation test ---------------------------------------------------------
res_enrich_MSigDB_c5 <- enricher(gene = geneList_OA_down_symbol,
                                 universe = universe_symbol,
                                 minGSSize = 10,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.1,
                                 TERM2GENE= MSigDB_c5
                                 )

png("Result/barplot_res_enrich_MSigDB_c5.png", height = 3000, width = 3000, res = 300)
barplot(res_enrich_MSigDB_c5, showCategory=20)
dev.off()
barplot(res_enrich_MSigDB_c5, showCategory=20)
write.table(res_enrich_MSigDB_c5@result, "Result/res_enrich_MSigDB_c5.txt",
            sep = "\t", row.names = F, quote = F)

res_enrich_GO <- enrichGO(gene = geneList_OA_down_Ensemble,
                          universe = universe_Ensambl,
                          OrgDb    = org.Mm.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.1)
png("Result/barplot_res_enrich_GO.png", height = 2000, width = 3000, res = 300)
barplot(res_enrich_GO, showCategory=10)
dev.off()
barplot(res_enrich_GO, showCategory=20)

## KEGG over-representation test ---------------------------------------------------------
res_enrich_KEGG <- enrichKEGG(gene = geneList_OA_down_entrezgene,
                              organism = 'mmu',
                              minGSSize = 10,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.1)

res_enrich_KEGG <- enrichKEGG(gene = geneList_OA_up_entrezgene,
                              organism = 'mmu',
                              minGSSize = 10,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.1)

write.table(res_enrich_KEGG@result, "Result/res_enrich_KEGG.txt",
            sep = "\t", row.names = F, quote = F)
res_enrich_KEGG@result


browseKEGG(res_enrich_KEGG, res_enrich_KEGG@result$ID[1])
res_enrich_KEGG@result$geneID



## Meshr enrichent analysis --------------------------------------------------------------

# rawdata_paper %>%
#     getBM(attributes = c("ensembl_gene_id","external_gene_name", "entrezgene"),
#           filters = "ensembl_gene_id",
#           values = .$id,
#           mart = Mg) -> tmp
# write.table(tmp, "Result/df_result__paper_allgene.txt",
#             sep = "\t", quote = FALSE, col.names = T, row.names = F)
df_result_paper_allgene <- read_delim("Result/df_result__paper_allgene.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(df_result_paper_allgene)[1]  <-　"gene_id"
rawdata_paper %>% 
    dplyr::filter(BH.qvalue < 0.25) %>%  
    left_join(., df_result_paper_allgene) %>% 
    dplyr::select(gene_id,
                  entrezgene,
                  external_gene_name,
                  logFC) -> d
d %>% dplyr::filter(logFC < 0) -> tmp
geneList_Paper_down_Ensemble   <- tmp$gene_id                      # for OA down, with Ensemble
geneList_Paper_down_symbol     <- toupper(tmp$external_gene_name)  # for OA down, with Gene symbol
tmp %>% 
    dplyr::select(entrezgene) %>%
    na.omit() %>% 
    .[!duplicated(.$entrezgene),] -> tmp
geneList_Paper_down_entrezgene <- tmp$entrezgene # for OA down, with Gene symbol

df_result_paper_allgene %>% 
    na.omit() %>% 
    .[!duplicated(.$entrezgene),] -> tmp
universe_paper_entranz <- tmp$entrezgene

# Up Only -------------------------------------------------------
meshParams_up <- new("MeSHHyperGParams",   
                  geneIds = geneList_OA_up_entrezgene, 
                  universeGeneIds = universe_entranz,
                  annotation = "MeSH.Mmu.eg.db",
                  category = "F",
                  database = "gene2pubmed",
                  pvalueCutoff = 1,
                  pAdjust = "BH")
res_meshr_up <- meshHyperGTest(meshParams_up)
head(summary(res_meshr_up))
res_meshr_up@ORA %>% 
    dplyr::filter(BH < 0.25) -> tmp
unique(tmp$MESHTERM)

# Down Only -------------------------------------------------------
meshParams_down <- new("MeSHHyperGParams",   
                     geneIds = geneList_OA_down_entrezgene, 
                     universeGeneIds = universe_entranz,
                     annotation = "MeSH.Mmu.eg.db",
                     category = "F",
                     database = "gene2pubmed",
                     pvalueCutoff = 1,
                     pAdjust = "BH")
res_meshr_down <- meshHyperGTest(meshParams_down)
head(summary(res_meshr_down))
res_meshr_down@ORA %>% 
    dplyr::filter(BH < 0.25) -> tmp
unique(tmp$MESHTERM)

# All --------------------------------------------------------------
meshParams_all <- new("MeSHHyperGParams", 
                  geneIds = c(geneList_OA_up_entrezgene,geneList_OA_down_entrezgene), 
                  universeGeneIds = universe_entranz,
                  annotation = "MeSH.Mmu.eg.db",
                  category = "F",
                  database = "gene2pubmed",
                  pvalueCutoff = 1,
                  pAdjust = "BH")
meshParams_all <- meshHyperGTest(meshParams_all)
head(summary(meshParams_all))
meshParams_all@ORA %>% 
    dplyr::filter(BH < 0.3) -> tmp
unique(tmp$MESHTERM)

# All --------------------------------------------------------------
meshParams_paper <- new("MeSHHyperGParams", 
                  geneIds = geneList_Paper_down_entrezgene, 
                  universeGeneIds = universe_paper_entranz,
                  annotation = "MeSH.Mmu.eg.db",
                  category = "F",
                  database = "gene2pubmed",
                  pvalueCutoff = 1,
                  pAdjust = "BH")
meshR <- meshHyperGTest(meshParams_paper)
head(summary(meshR))
meshR@ORA %>% 
    dplyr::filter(BH < 0.3) -> tmp
unique(tmp$MESHTERM)


meshR <- meshHyperGTest(meshParams)
write.table(meshR@ORA, "Result/meshr.txt",
            sep = "\t", row.names = F, quote = F)


## hypergeometrix with database --------------------------------------------------------------
Schizo_gene_list <- read_delim("Webresult/Schizo_gene_list.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
head(Schizo_gene_list)

X <- c(geneList_OA_down_symbol, geneList_OA_up_symbol)
X <- data.frame("Gene" = X)
X$Gene <- as.character(X$Gene)
str(X)
N <- universe_symbol
N <- data.frame("Gene" = N)
N$Gene <- as.character(N$Gene)
str(N)
Schizo_gene_list
SZ <- Schizo_gene_list %>% 
    dplyr::select(Gene) %>% 
    data.frame()
colnames(SZ) <- "Gene"
str(SZ)

dim(N)
dim(SZ)
SZ %>% 
    inner_join(., N) %>% 
    dim() -> n

SZ %>% 
    inner_join(., X) %>% 
    dim() -> x

N <- dim(N)[1]
N
n <- n[1]
n
X <- dim(X)[1]
X
x <- x[1]
x
sum(dhyper(x = x:X,
           m = n,
           n = N-n,
           k = X
           ))

## hypergeometrix with same database --------------------------------------------------------
ASD_gene <- read_csv("Webresult/ASD_gene.txt", 
                     col_names = FALSE) %>% data.frame()
ASD_gene$X1 <- as.character(ASD_gene$X1)
X <- c(geneList_OA_down_symbol, geneList_OA_up_symbol)
X <- data.frame("Gene" = X)
X$Gene <- as.character(X$Gene)
str(X)
N <- universe_symbol
N <- data.frame("Gene" = N)
N$Gene <- as.character(N$Gene)
str(N)

colnames(ASD_gene) <- "Gene"
ASD_gene$Gene <- toupper(ASD_gene$Gene)
str(ASD_gene)

dim(N)
dim(ASD_gene)
ASD_gene %>% 
    inner_join(., N) %>% 
    dim() -> n
n
ASD_gene %>% 
    inner_join(., X) %>% 
    dim() -> x
x
N <- dim(N)[1]
N
n <- n[1]
n
X <- dim(X)[1]
X
x <- x[1]
x

sum(dhyper(x = x:X,
           m = n,
           n = N-n,
           k = X
))


x <- data.frame(
    cell   = c("A", "A", "B", "B", "C", "C"),
    sample = c("A1", "A2", "B1", "B2", "C1", "C2"),
    weight = c(0.32, 0.33, 0.21, 0.22, 0.37, 0.36)
)
g <- ggplot(x, aes(x = cell, y = weight, fill = sample))
g <- g + geom_bar(stat = "identity")
g <- g + scale_fill_nejm()
plot(g)
p <- ggplot(mtcars, aes(x = as.factor(gear), fill = as.factor(vs))) +
    geom_bar()
p


## Compare DAIV vs GSEA ------------------------------------------------------------------------------------------------
DAVIT_Down_BP %>% dplyr::filter(Benjamini < 0.4) -> DAVIT
dim(DAVIT)
res_GSEA@result %>% 
    dplyr::select(Description, qvalues) %>% 
    dplyr::filter(qvalues < 0.4)-> res_GSEA_cut
dim(res_GSEA_cut)
res_GSEA_cut$Description
res_GSEA_cut$Description <- grep_chr_2(res_GSEA_cut$Description, "_")
res_GSEA_cut$Description <- gsub("_", " ", res_GSEA_cut$Description)
res_GSEA_cut$Description <- tolower(res_GSEA_cut$Description)
head(res_GSEA_cut$Description)
colnames(DAVIT)[2] <- "Term"
colnames(res_GSEA_cut)[1] <- "Term"
inner_join(DAVIT, res_GSEA_cut) -> tmp
write.table(tmp, "Result/marged_DAVID_and_GSEA.txt",
            sep = "\t", quote = F, row.names = F)


a <- DAVIT$Term
b <- res_GSEA_cut$Term

Venn_data <- list(DAVIT=a, GSEA=b)
venn.diagram(Venn_data, 
             filename = "Result/Venn_DAVITvsGSEA.tiff",
             imagetype = "tiff",
             height=3000, width=3000,
             resolution = 500,
             units = "px",
             fill=c(4,7),
             cat.pos=c(0, 25),
             cat.dist=c(0.02,0.05),
             cat.cex=c(2,2),
             cex=c(2,2,2)
)

res_enrich_GO@result %>%
    dplyr::filter(qvalue < 0.4) -> res_enrich_GO_cut
colnames(res_enrich_GO_cut)[1] <- "Term_id"
inner_join(DAVIT, res_enrich_GO_cut) -> tmp
write.table(tmp, "Result/marged_DAVID_and_enrichGO.txt",
            sep = "\t", quote = F, row.names = F)

a <- DAVIT$Term_id
b <- res_enrich_GO_cut$Term_id

Venn_data <- list(DAVIT=a, enrichGO=b)
venn.diagram(Venn_data, 
             filename = "Result/Venn_DAVITvsenrichGO.tiff",
             imagetype = "tiff",
             height=3000, width=3000,
             resolution = 500,
             units = "px",
             fill=c(4,7),
             cat.pos=c(355, 25),
             cat.dist=c(0.1,0.04),
             cat.cex=c(2,2),
             cex=c(2,2,2)
)

res_enrich_MSigDB_c5@result %>%
    dplyr::filter(qvalue < 0.4) -> res_enrich_MSigDB_c5_cut
colnames(res_enrich_MSigDB_c5_cut)[1] <- "Term"
head(res_enrich_MSigDB_c5_cut$Term)
res_enrich_MSigDB_c5_cut$Term <- grep_chr_2(res_enrich_MSigDB_c5_cut$Term, "_") %>% 
    gsub("_", " ", .) %>% 
    tolower()
head(res_enrich_MSigDB_c5_cut$Term)

inner_join(DAVIT, res_enrich_MSigDB_c5_cut) -> tmp
tmp
write.table(tmp, "Result/marged_DAVID_and_MSigDB_c5.txt",
            sep = "\t", quote = F, row.names = F)

a <- DAVIT$Term
b <- res_enrich_MSigDB_c5_cut$Term

Venn_data <- list(DAVIT=a, enricher=b)
venn.diagram(Venn_data, 
             filename = "Result/Venn_DAVITvsenricher.tiff",
             imagetype = "tiff",
             height=3000, width=3000,
             resolution = 500,
             units = "px",
             fill=c(4,7),
             cat.pos=c(355, 25),
             cat.dist=c(0.1,0.04),
             cat.cex=c(2,2),
             cex=c(2,2,2)
)