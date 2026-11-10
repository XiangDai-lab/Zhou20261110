# 加载必要的 R 包
library(Seurat)
library(scop)
library(data.table)
library(Matrix)
library(ggplot2)
library(scTenifoldKnk)
library(ggrepel)
library(igraph)

# 设置工作目录到你存放数据的路径
setwd("/mnt/disk3/xulutao/zy.qlxa")

# ==========================================
# 1. 读取并高度纯化细胞注释 (仅保留 Tumor 细胞)
# ==========================================
cat("正在读取细胞注释文件并提取纯 Tumor 细胞...\n")
meta_data <- fread("GSE181294_scRNAseq.ano.csv.gz", data.table = FALSE)

# 修正 fread 表头问题
if (meta_data[1, "V2"] == "cells" || meta_data[1, "V2"] == "celltype") {
  meta_data <- meta_data[-1, ]
}
rownames(meta_data) <- meta_data$V1

target_tissue_col <- "V3"          # 样本名/组织来源 
celltype_col      <- "V2"          # 细胞类型

# 【核心更新】：不仅提取 -T- (肿瘤样本)，且严格限定 V2 为 "Tumor" 细胞
tumor_meta <- meta_data[grepl("-T-", meta_data[[target_tissue_col]]) & meta_data[[celltype_col]] == "Tumor", ]
tumor_barcodes <- rownames(tumor_meta)
tumor_samples <- unique(tumor_meta[[target_tissue_col]])

cat(sprintf("经过深度过滤，共锁定 %d 个纯肿瘤 (Tumor) 细胞，来自 %d 个肿瘤样本。\n", length(tumor_barcodes), length(tumor_samples)))

# ==========================================
# 2. 读取原始数据 (RAW data) 并合并
# ==========================================
raw_data_dir <- "RAW_data"
cat("正在逐个读取包含 Tumor 细胞的样本表达矩阵...\n")

srt_list <- list()
all_files <- list.files(raw_data_dir, full.names = TRUE)

for (sample_id in tumor_samples) {
  file_path <- all_files[grepl(paste0(sample_id, "\\.count\\.csv\\.gz"), all_files)]
  
  if (length(file_path) > 0) {
    cat(sprintf("  - 读取并处理样本: %s ...\n", sample_id))
    mat <- fread(file_path[1], data.table = FALSE)
    rownames(mat) <- mat[, 1]
    mat <- mat[, -1]
    
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
    
    if (!grepl(paste0("^", sample_id), colnames(mat)[1])) {
      colnames(mat) <- paste0(sample_id, "_", colnames(mat))
    }
    
    srt_list[[sample_id]] <- CreateSeuratObject(counts = mat)
  }
}

cat("正在合并所有 Tumor 细胞数据...\n")
if (length(srt_list) > 1) {
  srt <- merge(srt_list[[1]], y = srt_list[2:length(srt_list)])
  # 修复 Seurat V5 特性
  srt <- JoinLayers(srt)
} else if (length(srt_list) == 1) {
  srt <- srt_list[[1]]
} else {
  stop("未找到任何矩阵文件！")
}

rm(srt_list, mat)
gc()

# ==========================================
# 3. 提取目标 Tumor 细胞并添加 Metadata
# ==========================================
cat("将矩阵精确切割至仅剩 Tumor 细胞...\n")
valid_cells <- intersect(colnames(srt), tumor_barcodes)
srt <- subset(srt, cells = valid_cells)

# 提取分级信息 (High_Grade vs Low_Grade) 并写入 metadata
tumor_meta$Tumor_Grade <- ifelse(grepl("-HG", tumor_meta[[target_tissue_col]]), "High_Grade", 
                                 ifelse(grepl("-LG", tumor_meta[[target_tissue_col]]), "Low_Grade", "Other"))

srt <- AddMetaData(srt, metadata = tumor_meta[valid_cells, ])

# ==========================================
# 4. Scop 标准化与降维 (纯 Tumor 细胞群体)
# ==========================================
cat("正在对纯 Tumor 细胞运行 scop 核心分析流程 (质控、归一化、降维)...\n")
srt <- standard_scop(srt)

# ==========================================
# 5. 可视化 1：纯 Tumor 细胞的 UMAP 图
# ==========================================
cat("\n【图 1】正在生成纯 Tumor 细胞群的 UMAP 图像...\n")
# 因为所有细胞都是 Tumor，所以按肿瘤分级(Tumor_Grade)来着色，观察异质性
p_tumor_umap <- CellDimPlot(
  srt,
  group.by = "Tumor_Grade",  
  reduction = "UMAP",       
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
) + ggplot2::labs(title = "UMAP of Pure Tumor Cells (Colored by Grade)")
print(p_tumor_umap)
Sys.sleep(3)


# ==========================================
# 6. scTenifoldKnk: Tumor 特异性虚拟敲除分析
# ==========================================
cat("\n==================================================\n")
cat("【图 2, 3, 4】开始执行 Tumor 特异性 OR51E2 虚拟敲除分析 (scTenifoldKnk)...\n")

target_gene <- "OR51E2"

srt <- JoinLayers(srt)
counts_matrix <- GetAssayData(srt, layer = "counts", assay = "RNA")

cat("正在过滤低表达基因以加速 Tumor 网络构建...\n")
gene_expr_pct <- rowSums(counts_matrix > 0) / ncol(counts_matrix)

# 找出在至少 5% 细胞中表达的基因
keep_genes <- rownames(counts_matrix)[gene_expr_pct > 0.05]

# 强制保留我们要敲除的目标基因
if (!(target_gene %in% keep_genes)) {
  cat(sprintf(" -> [注意] %s 表达占比不足 5%%，强制将其保留以进行 KO 分析。\n", target_gene))
  keep_genes <- unique(c(keep_genes, target_gene))
} else {
  cat(sprintf(" -> [好消息] %s 在 Tumor 细胞中表达丰富，已自动保留。\n", target_gene))
}

filtered_counts <- counts_matrix[keep_genes, ]

if (target_gene %in% rownames(filtered_counts)) {
  cat("正在构建 Tumor 特异性调控网络并执行虚拟敲除...\n")
  cat("(60核心满载运算中，构建 10 个子网络，请稍候...)\n")
  
  set.seed(2026)
  
  ko_res <- scTenifoldKnk(countMatrix = filtered_counts, 
                          gKO = target_gene, 
                          nc_nNet = 10,      
                          nc_nCells = 500,   
                          nCores = 60)       
  
  dr_df <- ko_res$diffRegulation
  dr_df <- dr_df[order(dr_df$distance, decreasing = TRUE), ]
  top_genes <- head(dr_df$gene[dr_df$gene != target_gene & dr_df$p.adj < 0.05], 20)
  
  # --- KO图1：虚拟敲除显著性排序散点图 ---
  cat("绘制散点图...\n")
  dr_df$Significance <- ifelse(dr_df$p.adj < 0.05, "Significant", "Not Significant")
  dr_df$Significance[dr_df$gene == target_gene] <- "Knocked Out"
  dr_df$Label <- ifelse(dr_df$gene %in% c(target_gene, head(top_genes, 10)), dr_df$gene, "")
  
  p_ko_scatter <- ggplot(dr_df, aes(x = rank(-distance), y = distance, color = Significance)) +
    geom_point(alpha = 0.7) +
    geom_text_repel(aes(label = Label), size = 4, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c("Knocked Out" = "#E41A1C", "Significant" = "#FF7F00", "Not Significant" = "#999999")) +
    scale_y_log10() + 
    theme_classic() +
    labs(title = paste("Virtual Knockout of", target_gene, "in Tumor Cells"),
         subtitle = "Manifold Alignment Distance Plot",
         x = "Gene Rank (Descending by Distance)",
         y = "Distance (Perturbation Score, log10 scale)")
  print(p_ko_scatter)
  Sys.sleep(3)
  
  # --- KO图2：Top 扰动基因条形图 ---
  cat("绘制条形图...\n")
  top_df <- dr_df[dr_df$gene %in% top_genes, ]
  p_ko_bar <- ggplot(top_df, aes(x = reorder(gene, distance), y = distance, fill = -log10(p.adj))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = paste("Top 20 Genes Perturbed by", target_gene, "KO (Tumor Only)"),
         x = "Genes", y = "Perturbation Distance")
  print(p_ko_bar)
  Sys.sleep(3)
  
  # --- KO图3：基因互作网络图 ---
  cat("绘制网络图...\n")
  if ("tensorNetworks" %in% names(ko_res)) {
    wt_mat <- as.matrix(ko_res$tensorNetworks$WT)
    nodes_to_plot <- c(target_gene, head(top_genes, 15)) 
    nodes_to_plot <- intersect(nodes_to_plot, rownames(wt_mat))
    
    if (length(nodes_to_plot) > 1) {
      sub_mat <- wt_mat[nodes_to_plot, nodes_to_plot]
      diag(sub_mat) <- 0
      
      # 有向图避免报错
      g <- graph_from_adjacency_matrix(sub_mat, mode = "directed", weighted = TRUE)
      
      # 红色正调控，蓝色负调控
      E(g)$color <- ifelse(E(g)$weight > 0, "tomato", "steelblue")
      E(g)$weight <- abs(E(g)$weight)
      
      # 过滤噪音边
      if(length(E(g)$weight) > 0) {
        threshold <- quantile(E(g)$weight, 0.7)
        g <- delete_edges(g, which(E(g)$weight < threshold))
      }
      
      V(g)$color <- ifelse(V(g)$name == target_gene, "red", "lightblue")
      V(g)$size <- ifelse(V(g)$name == target_gene, 25, 15)
      E(g)$width <- E(g)$weight * 10
      
      plot(g, 
           main = paste("Ego Network around", target_gene, "in Tumor Cells"),
           vertex.label.cex = 0.8,
           vertex.label.color = "black",
           edge.arrow.size = 0.4,
           layout = layout_with_fr)
      cat("\n[恭喜] 肿瘤特异性的 scTenifoldKnk 虚拟敲除分析全部完成！\n")
    } else {
      cat("可用于绘制网络的显著节点不足。\n")
    }
  }
} else {
  cat(sprintf("[致命错误] 目标基因 %s 在提取的矩阵中完全找不到。\n", target_gene))
}