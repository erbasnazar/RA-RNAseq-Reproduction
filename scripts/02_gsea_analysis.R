# =============================================================================
# 02_gsea_analysis.R
#
# Romatoid Artritte İskelet Kası Transkriptomik Analizi
# GEO: GSE283881 -- Gen Seti Zenginleştirme Analizi (GSEA), Hallmark gen setleri
#
# Ön koşul: 01_differential_expression.R dosyasının çalıştırılmış olması
# (bu script `resOrdered` nesnesini kullanır)
# =============================================================================

library(DESeq2)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# -----------------------------------------------------------------------------
# 1. Ensembl ID'lerin gen sembollerine (SYMBOL) eşlenmesi
# -----------------------------------------------------------------------------
res_df <- as.data.frame(resOrdered)
res_df$ENSEMBL <- rownames(res_df)

gene_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df$ENSEMBL,
  columns = c("SYMBOL"),
  keytype = "ENSEMBL"
)

geneList_df <- merge(
  res_df,
  gene_df,
  by = "ENSEMBL"
)

# -----------------------------------------------------------------------------
# 2. GSEA için sıralı (ranked) gen listesinin oluşturulması
#    Aynı sembole eşlenen genlerde en yüksek test istatistiği (stat) seçilir
# -----------------------------------------------------------------------------
tmp <- tapply(
  geneList_df$stat,
  geneList_df$SYMBOL,
  max
)

tmp <- sort(tmp, decreasing = TRUE)

geneList_symbol2 <- setNames(
  as.numeric(tmp),
  names(tmp)
)

# -----------------------------------------------------------------------------
# 3. MSigDB Hallmark gen setlerinin (H koleksiyonu) indirilmesi
# -----------------------------------------------------------------------------
hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hallmark_genes <- hallmark[, c("gs_name", "gene_symbol")]

# -----------------------------------------------------------------------------
# 4. GSEA analizinin çalıştırılması
# -----------------------------------------------------------------------------
gsea_result <- GSEA(
  geneList = geneList_symbol2,
  TERM2GENE = hallmark_genes,
  pvalueCutoff = 1,
  verbose = TRUE
)

gsea_df <- as.data.frame(gsea_result)
head(gsea_df)
dim(gsea_df)

# -----------------------------------------------------------------------------
# 5. Romatoid artrit/immün yanıt ile ilişkili yolakların filtrelenmesi
# -----------------------------------------------------------------------------
important <- gsea_df[
  grep(
    "INTERFERON_GAMMA|INTERFERON_ALPHA|TNFA|INFLAMMATORY|IL6_JAK_STAT3|MTORC1|HYPOXIA|UNFOLDED",
    gsea_df$Description
  ),
]

important[, c("Description", "NES", "p.adjust")]

# -----------------------------------------------------------------------------
# 6. Sıralı gen listesinin .rnk formatında dışa aktarılması
#    (GSEA Desktop veya başka araçlarla çapraz doğrulama için)
# -----------------------------------------------------------------------------
rnk <- data.frame(
  Gene = names(geneList_symbol2),
  Score = geneList_symbol2
)

write.table(
  rnk,
  file = "RA_IFN_vs_RA_NT.rnk",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
