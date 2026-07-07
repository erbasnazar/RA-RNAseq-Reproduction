RA-RNAseq-Reproduction
Reproduction of RNA-seq analyses from:
Oliver et al. (2025)
Communications Biology
Data: GEO GSE283881 — IFN (interferon-stimulated) vs NT (non-treated) samples, 12 vs 12
Project Goal
To reproduce and interpret transcriptomic analyses performed on rheumatoid arthritis skeletal muscle myobundles.
Scripts
scripts/01_differential_expression.R — Differential gene expression analysis using DESeq2
scripts/02_gsea_analysis.R — Gene symbol mapping, ranked gene list, GSEA with MSigDB Hallmark gene sets
Analyses
PCA
Heatmap
Volcano Plot
GSEA
Tools
R
DESeq2
clusterProfiler
msigdbr
org.Hs.eg.db
EnhancedVolcano
pheatmap
Main Findings
IFN-γ is the dominant transcriptional regulator.
RA-specific enrichment of Hypoxia, MTORC1 and UPR pathways.
JAK/STAT signaling plays a central role in muscle dysfunction.
