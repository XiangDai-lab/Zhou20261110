# 加载必要的 R 包
library(Seurat)
library(scop)
library(data.table)

# 设置工作目录到你存放数据的路径
setwd("/mnt/disk3/xulutao/zy.qlxa")

# ==========================================
# 1. 读取并处理作者提供的细胞注释 (Metadata)
# ==========================================
cat("正在读取细胞注释文件...\n")
# 读取注释文件，假设第一列是 Cell Barcode
meta_data <- fread("GSE181294_scRNAseq.ano.csv.gz", data.table = FALSE)
rownames(meta_data) <- meta_data[, 1] 

# ---------------------------------------------------------
# 【测试模块】：确认未知的列名变量
# ---------------------------------------------------------
cat("\n==================================================\n")
cat("【测试输出】查看 meta_data 的列名：\n")
print(colnames(meta_data))

cat("\n【测试输出】查看 meta_data 的前 5 行数据：\n")
print(head(meta_data, 5))
cat("==================================================\n")

# !!! 测试专用阻断 !!!
# 运行到这里会暂停，方便你查看控制台输出的列名。
# 确认后，请修改下方的 target_tissue_col 等参数，然后将下面这行 stop() 注释掉。
stop("【测试暂停】请查看控制台输出的列名和数据。确认无误并修改下方变量后，注释此行再重新运行！")
# ---------------------------------------------------------

# 【需要修改】：根据上面的测试输出，修改为正确的列名和字段
target_tissue_col <- "Tissue"          # <-- 请替换为真实的表示组织来源的列名 (例如：'tissue_type' 或 'disease')
target_tumor_name <- "Prostate Tumor"  # <-- 请替换为真实的表示前列腺肿瘤的字段名
celltype_col      <- "CellType"        # <-- 请替换为真实的作者注释细胞类型的列名 (例如：'cell_type' 或 'annotation')

# 仅保留前列腺肿瘤细胞的 metadata
tumor_meta <- meta_data[meta_data[[target_tissue_col]] == target_tumor_name, ]
tumor_barcodes <- rownames(tumor_meta)
cat(sprintf("筛选出 %d 个前列腺肿瘤细胞。\n", length(tumor_barcodes)))


# ==========================================
# 2. 读取原始数据 (RAW data)
# ==========================================
raw_data_dir <- "RAW_data"

cat("正在读取表达矩阵...\n")
# 因为你解压出来的是包含各个样本的独立的 .count.csv.gz 和 .mtx.gz，
# 建议通过 Read10X 或者直接分别读取合并（如果后续报错，我们再调整这里的读取方式）
# 暂时保留常规的 10x 读取方法，如果作者是标准 10x 格式通常能识别：
tryCatch({
  counts_matrix <- Read10X(data.dir = raw_data_dir) 
}, error = function(e) {
  stop("Read10X 读取失败，RAW data 可能非标准 10X 输出结构，请先确认文件读取逻辑。")
})

# 创建 Seurat 对象，并同时添加所有的 metadata
srt <- CreateSeuratObject(counts = counts_matrix, meta.data = meta_data)

# ==========================================
# 3. 筛选数据：仅保留前列腺肿瘤细胞
# ==========================================
cat("正在提取前列腺肿瘤细胞...\n")
# 根据前面提取的肿瘤细胞 barcodes 进行取子集
srt <- subset(srt, cells = tumor_barcodes)

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
  group.by = celltype_col,  # 使用我们在上面定义的细胞类型列
  reduction = "UMAP",       # standard_scop 默认生成的降维名称通常为 UMAP 或者 StandardUMAP2D
  label = TRUE,             # 在图中显示细胞类型的标签
  label_insitu = TRUE,      # 标签位于对应细胞群的中心
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

# 打印图像到 R 的 Plot 区域
print(p)
cat("分析完成！图像已输出到 Plot 区域。\n")