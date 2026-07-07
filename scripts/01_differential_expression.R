# =============================================================================
# 01_differential_expression.R
#
# Romatoid Artritte İskelet Kası Transkriptomik Analizi
# GEO: GSE283881 -- IFN (interferon-stimulated) vs NT (non-treated) karşılaştırması
#
# Amaç: Ham count matrisinden DESeq2 ile diferansiyel gen ifadesi analizi
# =============================================================================

library(DESeq2)

# -----------------------------------------------------------------------------
# 1. Count matrisinin okunması
#    Satırlar: genler (Ensembl ID) | Sütunlar: örnekler
# -----------------------------------------------------------------------------
counts <- read.delim(
  file.choose(),
  row.names = 1,
  check.names = FALSE
)

# -----------------------------------------------------------------------------
# 2. Örnek gruplarının tanımlanması (12 IFN + 12 NT örnek)
# -----------------------------------------------------------------------------
condition <- factor(c(
  rep("IFN", 12),
  rep("NT", 12)
))

colData <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

# -----------------------------------------------------------------------------
# 3. DESeq2 nesnesinin oluşturulması
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# -----------------------------------------------------------------------------
# 4. Düşük ekspresyonlu genlerin filtrelenmesi
#    (tüm örnekler toplamında en az 10 okuma sayısı olan genler tutulur)
# -----------------------------------------------------------------------------
dds <- dds[rowSums(counts(dds)) >= 10, ]

# -----------------------------------------------------------------------------
# 5. Diferansiyel ekspresyon analizinin çalıştırılması
# -----------------------------------------------------------------------------
dds <- DESeq(dds)

# -----------------------------------------------------------------------------
# 6. Sonuçların alınması ve düzenlenmesi (adjusted p-value'ya göre sıralı)
# -----------------------------------------------------------------------------
res <- results(dds)
resOrdered <- res[order(res$padj), ]

summary(res)
head(as.data.frame(resOrdered))

# Bir sonraki script (02_gsea_analysis.R) `resOrdered` nesnesini kullanır.
