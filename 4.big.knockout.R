# 加载必要的 R 包
library(Seurat)
library(scop)
library(data.table)
library(Matrix)
library(ggplot2) # 用于微调图表细节
library(scTenifoldKnk) # 用于单细胞虚拟敲除
library(ggrepel) # 用于散点图标签防重叠
library(igraph)  # 用于绘制互作网络图

# 设置工作目录到你存放数据的路径
setwd("/mnt/disk3/xulutao/zy.qlxa")

# ==========================================
# 1. 读取并处理作者提供的细胞注释 (Metadata)
# ==========================================
cat("正在读取细胞注释文件...\n")
meta_data <- fread("GSE181294_scRNAseq.ano.csv.gz", data.table = FALSE)

# 修正 fread 没有识别到正确表头的问题：
# 你的输出显示第一行是列名变量（cells, sample），第二行才是真正数据。
if (meta_data[1, "V2"] == "cells" || meta_data[1, "V2"] == "celltype") {
  meta_data <- meta_data[-1, ]
}
# 设置行名为 V1 列（Cell Barcode）
rownames(meta_data) <- meta_data$V1

# 根据测试输出确认的变量名：
target_tissue_col <- "V3"          # V3 列对应样本名/组织来源 (例如 SCG-PCA21-N-LG)
celltype_col      <- "V2"          # V2 列对应细胞类型 (例如 Mast cells)

# 提取前列腺肿瘤细胞 (-T- 代表 Tumor 样本)
tumor_meta <- meta_data[grepl("-T-", meta_data[[target_tissue_col]]), ]
tumor_barcodes <- rownames(tumor_meta)
tumor_samples <- unique(tumor_meta[[target_tissue_col]])
cat(sprintf("筛选出 %d 个前列腺肿瘤细胞，来自 %d 个肿瘤样本。\n", length(tumor_barcodes), length(tumor_samples)))

# ==========================================
# 2. 读取原始数据 (RAW data) 并合并
# ==========================================
raw_data_dir <- "RAW_data"
cat("正在逐个读取前列腺肿瘤样本的表达矩阵并合并 (这可能需要几分钟)...\n")

srt_list <- list()
all_files <- list.files(raw_data_dir, full.names = TRUE)

for (sample_id in tumor_samples) {
  # 匹配对应的 count.csv.gz 文件 (例如: GSM5494348_SCG-PCA21-T-LG.count.csv.gz)
  file_path <- all_files[grepl(paste0(sample_id, "\\.count\\.csv\\.gz"), all_files)]
  
  if (length(file_path) > 0) {
    cat(sprintf("  - 读取样本: %s ...\n", sample_id))
    mat <- fread(file_path[1], data.table = FALSE)
    # 第一列是基因名，设置为行名
    rownames(mat) <- mat[, 1]
    mat <- mat[, -1]
    
    # 转换为稀疏矩阵以节省内存，这是处理 scRNA-seq 的推荐做法
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
    
    # 确保列名(barcode)带有样本前缀，和 metadata 中的格式保持一致
    if (!grepl(paste0("^", sample_id), colnames(mat)[1])) {
      colnames(mat) <- paste0(sample_id, "_", colnames(mat))
    }
    
    srt_list[[sample_id]] <- CreateSeuratObject(counts = mat)
  } else {
    cat(sprintf("  [警告] 找不到样本 %s 的文件，已跳过。\n", sample_id))
  }
}

cat("正在合并所有肿瘤样本数据...\n")
if (length(srt_list) > 1) {
  srt <- merge(srt_list[[1]], y = srt_list[2:length(srt_list)])
  # 【修复 Seurat V5 特性】：合并后将不同样本的独立 layers (counts.1, counts.2...) 统一组合为单一 counts layer
  srt <- JoinLayers(srt)
} else if (length(srt_list) == 1) {
  srt <- srt_list[[1]]
} else {
  stop("未在 RAW_data 中找到任何肿瘤样本对应的表达矩阵，请检查解压路径！")
}

# 清理内存
rm(srt_list, mat)
gc()

# ==========================================
# 3. 过滤数据并添加 metadata
# ==========================================
cat("添加细胞注释并提取目标肿瘤细胞...\n")
# 交集获取在表达矩阵和 metadata 中都存在的细胞（双保险）
valid_cells <- intersect(colnames(srt), tumor_barcodes)
srt <- subset(srt, cells = valid_cells)

# 仅向 Seurat 对象中加入肿瘤相关的 metadata
srt <- AddMetaData(srt, metadata = tumor_meta[valid_cells, ])

# ==========================================
# 4. 使用 scop 包的标准 pipeline 进行降维聚类
# ==========================================
cat("正在运行 scop 核心分析流程 (包含质控、归一化、降维、UMAP)...\n")
# standard_scop() 是 scop 包提供的一键式单细胞处理流程
srt <- standard_scop(srt)

# ==========================================
# 5. 可视化：绘制并打印全局细胞 UMAP 图
# ==========================================
cat("正在生成全细胞 UMAP 图像...\n")
p_cell <- CellDimPlot(
  srt,
  group.by = celltype_col,  
  reduction = "UMAP",       
  label = TRUE,             
  label_insitu = TRUE,      
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
print(p_cell)


# ==========================================
# 6. 特异性基因 (OR51E2) 的表达模式探究
# ==========================================
cat("\n==================================================\n")
cat("开始探索 OR51E2 基因的表达情况...\n")

target_gene <- "OR51E2"

if (target_gene %in% rownames(srt)) {
  
  srt$Tumor_Grade <- ifelse(grepl("-HG", srt[[target_tissue_col]]), "High_Grade", 
                            ifelse(grepl("-LG", srt[[target_tissue_col]]), "Low_Grade", "Other"))
  
  cat("1. 正在生成 [全局基因表达 UMAP 图]...\n")
  p_gene_overall <- FeatureDimPlot(srt, features = target_gene, reduction = "UMAP", theme_use = "theme_blank")
  print(p_gene_overall)
  Sys.sleep(2)
  
  cat("2. 正在生成 [不同肿瘤分级的拆分表达 UMAP 图]...\n")
  p_gene_split <- FeatureDimPlot(srt, features = target_gene, reduction = "UMAP", split.by = "Tumor_Grade", theme_use = "theme_blank")
  print(p_gene_split)
  Sys.sleep(2)
  
  cat("3. 正在生成 [细胞类型的定量小提琴图]...\n")
  p_vln <- FeatureStatPlot(srt, group.by = celltype_col, bg.by = "Tumor_Grade", stat.by = target_gene, add_box = TRUE) + 
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_vln)
  Sys.sleep(2)
  
} else {
  cat(sprintf("[警告] 在表达矩阵中未找到基因: %s。\n", target_gene))
}

# ==========================================
# 7. scTenifoldKnk 虚拟敲除分析 (Virtual KO)
# ==========================================
cat("\n==================================================\n")
cat("7. 开始执行目标基因的单细胞虚拟敲除分析 (scTenifoldKnk)...\n")

if (target_gene %in% rownames(srt)) {
  # 提取表达矩阵: 修复 Seurat V5 的参数更新 (slot -> layer)
  # 双重保险：提取前强制检查并合并所有 layer
  srt <- JoinLayers(srt)
  counts_matrix <- GetAssayData(srt, layer = "counts", assay = "RNA")
  
  # 【性能优化】过滤极低表达基因，以加速张量网络计算
  cat("正在过滤低表达基因以加速构建网络...\n")
  gene_expr_pct <- rowSums(counts_matrix > 0) / ncol(counts_matrix)
  
  # 找出在至少 5% 细胞中表达的基因
  keep_genes <- rownames(counts_matrix)[gene_expr_pct > 0.05]
  
  # 【核心修复】：强制保留我们要敲除的目标基因 (即使它的表达率低于 5%)
  if (!(target_gene %in% keep_genes)) {
    cat(sprintf(" -> [注意] %s 表达细胞占比不足 5%%，已强制将其保留在矩阵中以进行 KO 分析。\n", target_gene))
    keep_genes <- unique(c(keep_genes, target_gene))
  }
  
  # 进行矩阵子集提取
  filtered_counts <- counts_matrix[keep_genes, ]
  
  if (target_gene %in% rownames(filtered_counts)) {
    cat("构建单细胞调控张量网络并执行虚拟敲除...\n")
    cat("(算力全开模式: 60核心并行，构建 10 个子网络获得高精度结果，这可能需要一段时间...)\n")
    
    # 设定随机种子以保证结果可重复
    set.seed(2026)
    
    # 运行虚拟敲除
    # 充分利用 64核 512G 服务器算力
    ko_res <- scTenifoldKnk(countMatrix = filtered_counts, 
                            gKO = target_gene, 
                            nc_nNet = 10,      # 提高到 10 个子网络，获取更高精度和稳定性的结果
                            nc_nCells = 500,   # 维持默认 500 足够代表性
                            nCores = 60)       # 启用 60 个线程（留4个核心给系统防止无响应）
    
    # 提取差异调节(受扰动)结果
    dr_df <- ko_res$diffRegulation
    
    # 获取受影响最显著的 Top 20 靶基因
    dr_df <- dr_df[order(dr_df$distance, decreasing = TRUE), ]
    top_genes <- head(dr_df$gene[dr_df$gene != target_gene & dr_df$p.adj < 0.05], 20)
    
    # --- KO图1：虚拟敲除显著性排序散点图 (Manhattan-like plot) ---
    cat("绘制受虚拟敲除影响的基因排序散点图...\n")
    dr_df$Significance <- ifelse(dr_df$p.adj < 0.05, "Significant", "Not Significant")
    dr_df$Significance[dr_df$gene == target_gene] <- "Knocked Out"
    dr_df$Label <- ifelse(dr_df$gene %in% c(target_gene, head(top_genes, 10)), dr_df$gene, "")
    
    p_ko_scatter <- ggplot(dr_df, aes(x = rank(-distance), y = distance, color = Significance)) +
      geom_point(alpha = 0.7) +
      geom_text_repel(aes(label = Label), size = 4, max.overlaps = 20, show.legend = FALSE) +
      scale_color_manual(values = c("Knocked Out" = "#E41A1C", "Significant" = "#FF7F00", "Not Significant" = "#999999")) +
      theme_classic() +
      labs(title = paste("Virtual Knockout of", target_gene),
           subtitle = "Manifold Alignment Distance Plot",
           x = "Gene Rank (Descending by Distance)",
           y = "Distance (Perturbation Score)")
    print(p_ko_scatter)
    Sys.sleep(2)
    
    # --- KO图2：Top 扰动基因条形图 ---
    cat("绘制 Top 下游靶标基因条形图...\n")
    top_df <- dr_df[dr_df$gene %in% top_genes, ]
    p_ko_bar <- ggplot(top_df, aes(x = reorder(gene, distance), y = distance, fill = -log10(p.adj))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_viridis_c() +
      theme_minimal() +
      labs(title = paste("Top 20 Genes Perturbed by", target_gene, "KO"),
           x = "Genes", y = "Perturbation Distance")
    print(p_ko_bar)
    Sys.sleep(2)
    
    # --- KO图3：基因互作网络图 (Ego-Network) ---
    cat("绘制目标基因互作网络图...\n")
    if ("tensorNetworks" %in% names(ko_res)) {
      # 提取野生型(WT)网络的降维/平滑后张量矩阵作为邻接矩阵
      wt_mat <- as.matrix(ko_res$tensorNetworks$WT)
      # 选取被敲除基因及其 Top 受影响基因组成子网络
      nodes_to_plot <- c(target_gene, head(top_genes, 15)) 
      
      # 防御性编程：确保这些基因确实存在于网络矩阵中
      nodes_to_plot <- intersect(nodes_to_plot, rownames(wt_mat))
      
      if (length(nodes_to_plot) > 1) {
        sub_mat <- wt_mat[nodes_to_plot, nodes_to_plot]
        # 去除自环（对角线设为0）
        diag(sub_mat) <- 0
        
        # 将矩阵转化为 igraph 对象
        g <- graph_from_adjacency_matrix(sub_mat, mode = "undirected", weighted = TRUE)
        
        # 简单修饰一下网络图的外观
        V(g)$color <- ifelse(V(g)$name == target_gene, "red", "lightblue")
        V(g)$size <- ifelse(V(g)$name == target_gene, 25, 15)
        # 边的粗细和权重成正比
        E(g)$width <- E(g)$weight * 5
        
        plot(g, 
             main = paste("Ego Network around", target_gene),
             vertex.label.cex = 0.8,
             vertex.label.color = "black",
             edge.color = "gray80",
             layout = layout_with_fr) # Fruchterman-Reingold 弹簧布局
        cat("scTenifoldKnk 虚拟敲除分析全部完成！\n")
      } else {
        cat("可用于绘制网络的显著节点不足。\n")
      }
    }
    
  } else {
    cat(sprintf("基因 %s 在过滤后的矩阵中缺失 (可能是表达量极低)。为加速计算，该基因被剔除了，无法进行 KO 分析。\n", target_gene))
  }
}