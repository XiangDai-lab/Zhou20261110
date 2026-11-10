# 加载必要的 R 包
library(Seurat)
library(scop)
library(data.table)
library(Matrix)

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
# 它会自动计算 PCA 和 UMAP 降维
srt <- standard_scop(srt)

# ==========================================
# 5. 可视化：绘制并打印 UMAP 图
# ==========================================
cat("正在生成 UMAP 图像...\n")
# 使用 scop 包特有的 CellDimPlot 函数绘制带有作者细胞注释的 UMAP 图
p <- CellDimPlot(
  srt,
  group.by = celltype_col,  # 映射我们在上面定义的细胞类型列 (V2)
  reduction = "UMAP",       # standard_scop 默认生成的降维名称
  label = TRUE,             # 在图中显示细胞类型的标签
  label_insitu = TRUE,      # 标签位于对应细胞群的中心
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

# 打印图像到 R 的 Plot 区域
print(p)
cat("分析完成！图像已输出到 Plot 区域。\n")