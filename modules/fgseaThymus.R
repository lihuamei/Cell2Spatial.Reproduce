library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(ggplot2)

data <- openxlsx::read.xlsx("for.gesa.cell2spatial.xlsx", rowNames = T) %>% rename(SYMBOL = gene)
data <- subset(data, cluster == "Medulla")
gene <- data$SYMBOL
gene <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gene <- dplyr::distinct(gene, SYMBOL, .keep_all = TRUE)

data_all <- data %>%
    inner_join(gene, by = "SYMBOL")
dim(data_all)
head(data_all)

data_all_sort <- data_all %>%
    arrange(desc(avg_log2FC))
head(data_all_sort)

geneList <- data_all_sort$avg_log2FC
names(geneList) <- data_all_sort$ENTREZID
head(geneList)

kegg_gmt <- read.gmt("C:\\Users\\Administrator\\Downloads\\c5.all.v2023.1.Hs.entrez.gmt")
gsea <- GSEA(geneList, TERM2GENE = kegg_gmt)

openxlsx::write.xlsx(gsea, file = "medulla.fgsea.xlsx")


gseaplot2(gsea, c(4, 16))
ggsave("Medulla.gsea.pdf", width = 6, height = 4)


# =============================================

data <- openxlsx::read.xlsx("for.gesa.cell2spatial.xlsx", rowNames = T) %>% rename(SYMBOL = gene)
data <- subset(data, cluster == "Cortex")
gene <- data$SYMBOL
gene <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gene <- dplyr::distinct(gene, SYMBOL, .keep_all = TRUE)

data_all <- data %>%
    inner_join(gene, by = "SYMBOL")
dim(data_all)
head(data_all)

data_all_sort <- data_all %>%
    arrange(desc(avg_log2FC))
head(data_all_sort)

geneList <- data_all_sort$avg_log2FC
names(geneList) <- data_all_sort$ENTREZID
head(geneList)

kegg_gmt <- read.gmt("C:\\Users\\Administrator\\Downloads\\c5.all.v2023.1.Hs.entrez.gmt")
gsea <- GSEA(geneList, TERM2GENE = kegg_gmt)

openxlsx::write.xlsx(gsea, file = "cortex.fgsea.xlsx")


gseaplot2(gsea, c(302, 643))
ggsave("Cortex.gsea.pdf", width = 6, height = 4)
