run_sepsis_peroxisome_analysis <- function() {
    suppressPackageStartupMessages({
        library(Seurat)
        library(Matrix)
        library(data.table)
        library(ggplot2)
        library(patchwork)
        library(ggrepel)
        library(pheatmap)
    })
    client_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    input_dir <- file.path(client_root, "data")
    stage <- file.path(client_root, "analysis_work")
    stage_data <- file.path(stage, "data")
    out <- file.path(client_root, "result")
    input_manifest <- fread(file.path(input_dir, "input_manifest.tsv"))
    input_paths <- file.path(input_dir, input_manifest$file)
    stopifnot(all(file.exists(input_paths)))
    stopifnot(all(file.info(input_paths)$size == input_manifest$size_bytes))
    dir.create(stage_data, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(out)) 
        unlink(out, recursive = TRUE, force = TRUE)
    result_dirs <- c("1.QCFilter", "2.Cluster", "3.CellAnnotate", "4.ClusterGvis", "5.Signature", "6.CytoTRACE", "7.Monocle2", 
        "8.CellChat", "9.Machine learning", "10.GSEA", "11.GSVA", "12.Immune infiltration", "13.Immune association", "14.Motif", 
        "15.Non-coding RNA networks", "16.KeyGene expression abundance", "17.Immunometabolic pathways", "18.Drug Prediction")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    invisible(lapply(file.path(out, result_dirs), dir.create, recursive = TRUE, showWarnings = FALSE))
    stage_files <- c("expr_all_main.txt", "metadata_all_main.txt", "GSE290679-test-expr.txt", "GSE290679-test-metadata.txt", 
        "signature.txt", "hallmark.gs.rdata", "GenesetInfo.txt", "Immunomodulator_and_chemokines.txt", "mircode_highconsfamilies.txt", 
        "interactions.tsv", "cibersort.ref.txt", "113_ML_methods.txt")
    file.copy(file.path(input_dir, stage_files), file.path(stage_data, stage_files), overwrite = TRUE)
    all_expr <- fread(file.path(stage_data, "expr_all_main.txt"), check.names = FALSE)
    all_meta <- fread(file.path(stage_data, "metadata_all_main.txt"), header = FALSE)
    test_meta <- fread(file.path(stage_data, "GSE290679-test-metadata.txt"), header = FALSE)
    test_ids <- as.character(test_meta[[1]])
    train_ids <- setdiff(names(all_expr)[-1], test_ids)
    fwrite(all_expr[, c(names(all_expr)[1], train_ids), with = FALSE], file.path(stage_data, "expr.txt"), sep = "\t")
    fwrite(all_meta[V1 %in% train_ids], file.path(stage_data, "metadata.txt"), sep = "\t", col.names = FALSE)
    rm(all_expr, all_meta, test_meta)
    invisible(gc())
    setwd(stage)
    data_root <- stage_data
    keygenes <- c("LYL1", "PIM1", "RPS7")
    Sys.setenv(PIPELINE_BASE_DIR = out, ML113_EXTERNAL_TIMEOUT_REGEX = "")
    options(timeout = 600)
    library(Seurat)
    library(harmony)
    library(ggplot2)
    library(scRNAtoolVis)
    library(dplyr)
    human2mouse <- function(x) {
        x1 <- substr(x, 1, 1)
        x2 <- substr(x, 2, nchar(x))
        return(paste0(x1, tolower(x2)))
    }
    scalable_marker_sets <- list(`CD14 Mono` = c("LYZ", "LST1", "S100A8", "S100A9", "S100A12", "FCN1", "CD14"), `CD16 Mono` = c("FCGR3A", 
        "MS4A7", "LST1", "CTSS"), Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9"), `CD4 T` = c("CD3D", "CD3E", 
        "CD4", "IL7R", "CCR7", "LTB"), `CD8 T` = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GZMK"), NK = c("NKG7", "GNLY", 
        "KLRD1", "PRF1", "GZMB"), B = c("MS4A1", "CD79A", "CD79B", "CD74", "BANK1"), Plasma = c("MZB1", "JCHAIN", "IGHG1", 
        "XBP1"), DC = c("FCER1A", "CST3", "CLEC10A", "LILRA4", "HLA-DRA"), Platelet = c("PPBP", "PF4", "GP9", "ITGA2B"))
    score_marker_sets_sparse <- function(object, marker_sets = scalable_marker_sets) {
        mat <- GetAssayData(object, assay = "RNA", slot = "data")
        scores <- matrix(0, nrow = ncol(object), ncol = length(marker_sets), dimnames = list(colnames(object), names(marker_sets)))
        for (nm in names(marker_sets)) {
            genes <- intersect(marker_sets[[nm]], rownames(mat))
            if (length(genes)) 
                scores[, nm] <- Matrix::colMeans(mat[genes, , drop = FALSE])
        }
        scores
    }
    assign_marker_celltypes_sparse <- function(object) {
        scores <- score_marker_sets_sparse(object)
        best <- max.col(scores, ties.method = "first")
        max_score <- scores[cbind(seq_len(nrow(scores)), best)]
        labels <- colnames(scores)[best]
        labels[max_score <= 0] <- "Unknown"
        object$cellType_1 <- labels
        object$RNA_snn_res.0.2 <- factor(labels)
        object$seurat_clusters <- factor(labels)
        Idents(object) <- "cellType_1"
        object
    }
    make_marker_table_sparse <- function(object, celltype = "cellType_1", max_features = 1200) {
        mat <- GetAssayData(object, assay = "RNA", slot = "data")
        counts <- GetAssayData(object, assay = "RNA", slot = "counts")
        marker_candidates <- unique(unlist(scalable_marker_sets, use.names = FALSE))
        candidate <- unique(c(intersect(marker_candidates, rownames(mat)), head(VariableFeatures(object), max_features)))
        candidate <- intersect(candidate, rownames(mat))
        groups <- as.character(object@meta.data[[celltype]])
        out <- list()
        for (grp in sort(unique(groups))) {
            idx <- groups == grp
            if (sum(idx) < 10 || sum(!idx) < 10) 
                next
            mean_in <- Matrix::rowMeans(mat[candidate, idx, drop = FALSE])
            mean_out <- Matrix::rowMeans(mat[candidate, !idx, drop = FALSE])
            pct_in <- Matrix::rowMeans(counts[candidate, idx, drop = FALSE] > 0)
            pct_out <- Matrix::rowMeans(counts[candidate, !idx, drop = FALSE] > 0)
            df <- data.frame(p_val = NA_real_, avg_log2FC = as.numeric(mean_in - mean_out), pct.1 = as.numeric(pct_in), pct.2 = as.numeric(pct_out), 
                p_val_adj = NA_real_, cluster = grp, gene = candidate, stringsAsFactors = FALSE)
            df <- df[order(df$avg_log2FC, decreasing = TRUE), , drop = FALSE]
            out[[grp]] <- head(df, 50)
        }
        do.call(rbind, out)
    }
    set_embedding_from_reduction <- function(object) {
        emb <- NULL
        if ("harmony" %in% names(object@reductions)) {
            emb <- Embeddings(object, "harmony")[, 1:2, drop = FALSE]
        }
        else if ("pca" %in% names(object@reductions)) {
            emb <- Embeddings(object, "pca")[, 1:2, drop = FALSE]
        }
        if (is.null(emb)) {
            emb <- cbind(log1p(object$nCount_RNA), log1p(object$nFeature_RNA))
            rownames(emb) <- colnames(object)
        }
        colnames(emb) <- c("UMAP_1", "UMAP_2")
        object[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(object))
        object
    }
    large_fullcell_mode <- function(object) {
        threshold <- as.integer(Sys.getenv("SCALABLE_FULLCELL_THRESHOLD", "300000"))
        "no_sampling" %in% colnames(object@meta.data) && ncol(object) > threshold
    }
    seurat_module_large_fullcell <- function(scRNA, QC.path, Cluster.path, nfeatures = 2000, npcs = 50, PC = 20, assay = "RNA") {
        scRNA <- NormalizeData(object = scRNA)
        scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = nfeatures)
        top10 <- head(VariableFeatures(scRNA), 10)
        p1 <- VariableFeaturePlot(scRNA)
        p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
        ggsave(filename = paste0(Cluster.path, "/1.VariableFeaturePlot.pdf"), plot = p1 + p2, height = 7, width = 14)
        pca_status <- "not_run"
        harmony_status <- "not_run"
        scRNA <- tryCatch({
            obj <- ScaleData(object = scRNA, features = VariableFeatures(scRNA), vars.to.regress = c("percent.mt", "percent.ribo"))
            obj <- RunPCA(obj, verbose = TRUE, npcs = npcs, features = VariableFeatures(obj))
            p <- ElbowPlot(object = obj, ndims = npcs)
            ggsave(filename = paste0(Cluster.path, "/2.ElbowPlot.pdf"), plot = p, height = 7, width = 7)
            pca_status <<- "completed"
            if (length(unique(obj$orig.ident)) > 1) {
                obj <- RunHarmony(object = obj, group.by.vars = "orig.ident", assay.use = assay, verbose = FALSE)
                harmony_status <<- "completed"
            }
            obj
        }, error = function(e) {
            pca_status <<- paste0("failed: ", conditionMessage(e))
            scRNA
        })
        scRNA <- set_embedding_from_reduction(scRNA)
        scRNA <- assign_marker_celltypes_sparse(scRNA)
        p <- DimPlot(object = scRNA, reduction = "umap", label = FALSE, group.by = "orig.ident", raster = TRUE) + NoLegend()
        ggsave(filename = paste0(Cluster.path, "/4.HarmonyPlot.pdf"), plot = p, height = 7, width = 7)
        p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = "cellType_1", raster = TRUE) + NoLegend()
        ggsave(filename = paste0(Cluster.path, "/5.MarkerScoreCellType.pdf"), plot = p, height = 7, width = 7)
        audit <- data.frame(item = c("mode", "cells", "features", "variable_features", "pca_status", "batch_correction_status", 
            "annotation_method"), value = c("scalable_fullcell_no_sampling", ncol(scRNA), nrow(scRNA), length(VariableFeatures(scRNA)), 
            pca_status, harmony_status, "sparse_marker_score_per_cell"), check.names = FALSE)
        write.table(audit, file.path(Cluster.path, "0.scalable_fullcell_step1_audit.tsv"), sep = "\t", row.names = FALSE, 
            quote = FALSE)
        scRNA
    }
    seurat_module <- function(outpath = "./result/", SeuratObject = NA, min.features = 200, Idents = "RNA_snn_res.0.2", assay = "RNA", 
        set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 20, nfeatures = 2000, npcs = 50) {
        a <- list.files(outpath, pattern = "QCFilter", full.names = T)
        if (length(a) == 0) {
            QC.path <- paste0(outpath, "/QCFilter/")
            dir.create(QC.path)
        }
        else {
            QC.path <- a[1]
        }
        b <- list.files(outpath, pattern = "Cluster", full.names = T)
        if (length(b) == 0) {
            Cluster.path <- paste0(outpath, "/Cluster/")
            dir.create(Cluster.path)
        }
        else {
            Cluster.path <- b[1]
        }
        scRNA <- SeuratObject
        grep("^[M,m][T,t]-", rownames(scRNA), value = T)
        grep("^R[P,p][SL,sl]", rownames(scRNA), value = T)
        scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, pattern = "^[M,m][T,t]-")
        scRNA[["percent.ribo"]] <- PercentageFeatureSet(scRNA, pattern = "^R[P,p][SL,sl]")
        nomt <- ifelse(max(scRNA[["percent.mt"]]) == 0, T, F)
        if (nomt) {
            MT_genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
            MT_genes <- intersect(MT_genes, rownames(scRNA))
            scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, features = MT_genes)
        }
        mt_value <- stats::mad(scRNA$percent.mt)
        lower_mt <- median(scRNA$percent.mt) - 3 * mt_value
        upper_mt <- median(scRNA$percent.mt) + 3 * mt_value
        nF_value <- stats::mad(scRNA$nFeature_RNA)
        lower_nF <- median(scRNA$nFeature_RNA) - 3 * nF_value
        upper_nF <- median(scRNA$nFeature_RNA) + 3 * nF_value
        nC_value <- stats::mad(scRNA$nCount_RNA)
        lower_nC <- median(scRNA$nCount_RNA) - 3 * nC_value
        upper_nC <- median(scRNA$nCount_RNA) + 3 * nC_value
        cells_before_filter <- ncol(scRNA)
        prepared_full_qc <- "no_sampling" %in% colnames(scRNA@meta.data) && all(as.character(scRNA$no_sampling) %in% c("TRUE", 
            "true", "1"))
        skip_merged_global_qc <- prepared_full_qc || identical(tolower(Sys.getenv("SKIP_MERGED_GLOBAL_QC", "false")), "true")
        if (!skip_merged_global_qc) {
            scRNA <- subset(x = scRNA, subset = nFeature_RNA >= min.features & percent.mt <= upper_mt & nFeature_RNA <= upper_nF & 
                nCount_RNA <= upper_nC)
        }
        qc_gate <- data.frame(item = c("cells_before_module_filter", "cells_after_module_filter", "prepared_full_qc", "module_filter_action", 
            "reason"), value = c(cells_before_filter, ncol(scRNA), prepared_full_qc, ifelse(skip_merged_global_qc, "retained_all_prepared_cells", 
            "applied_merged_mad_filter"), ifelse(skip_merged_global_qc, "dataset_specific_qc_already_applied_in_prep", "legacy_single_object_qc")), 
            check.names = FALSE)
        write.table(qc_gate, file.path(QC.path, "0.full_cell_qc_gate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
        p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", 
            ncol = 4)
        ggsave(filename = paste0(QC.path, "/1.VlnPlot.pdf"), plot = p, width = 20, height = 4)
        plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) + NoLegend()
        plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) + NoLegend()
        plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", raster = F) + NoLegend()
        ggsave(filename = paste0(QC.path, "/2.FeatureScatter.pdf"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
        if (large_fullcell_mode(scRNA)) {
            return(seurat_module_large_fullcell(scRNA, QC.path, Cluster.path, nfeatures = nfeatures, npcs = npcs, PC = PC, 
                assay = assay))
        }
        scRNA <- NormalizeData(object = scRNA)
        scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = nfeatures)
        top10 <- head(VariableFeatures(scRNA), 10)
        plot1 <- VariableFeaturePlot(scRNA)
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        ggsave(filename = paste0(Cluster.path, "/1.VariableFeaturePlot.pdf"), plot = plot1 + plot2, height = 7, width = 14)
        scRNA <- ScaleData(object = scRNA, vars.to.regress = c("percent.mt", "percent.ribo"))
        scRNA <- RunPCA(scRNA, verbose = T, npcs = npcs, features = VariableFeatures(object = scRNA))
        p <- ElbowPlot(object = scRNA, ndims = npcs)
        ggsave(filename = paste0(Cluster.path, "/2.ElbowPlot.pdf"), plot = p, height = 7, width = 7)
        if (length(unique(scRNA$orig.ident)) > 1) {
            scRNA <- RunHarmony(object = scRNA, group.by.vars = "orig.ident", assay.use = assay, verbose = FALSE)
            scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:PC, verbose = T)
            scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "harmony", verbose = T)
            scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
            p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
            ggsave(filename = paste0(Cluster.path, "/3.PcaPlot.pdf"), plot = p, height = 7, width = 7)
            p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
            ggsave(filename = paste0(Cluster.path, "/4.HarmonyPlot.pdf"), plot = p, height = 7, width = 7)
        }
        else {
            scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:PC, verbose = T)
            scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "pca", verbose = T)
            scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
            p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
            ggsave(filename = paste0(Cluster.path, "/3.PcaPlot.pdf"), plot = p, height = 7, width = 7)
        }
        return(scRNA)
    }
    annotate_module <- function(outpath = "./result/", SeuratObject = NA, ref, labels, Idents = "RNA_snn_res.0.2") {
        c <- list.files(outpath, pattern = "CellAnnotate", full.names = T)
        if (length(c) == 0) {
            Anno.path <- paste0(outpath, "/CellAnnotate/")
            dir.create(Anno.path)
        }
        else {
            Anno.path <- c[1]
        }
        scRNA <- SeuratObject
        if (large_fullcell_mode(scRNA) && "cellType_1" %in% colnames(scRNA@meta.data)) {
            p <- DimPlot(scRNA, group.by = "cellType_1", label = TRUE, reduction = "umap", raster = TRUE) + NoLegend()
            ggsave(filename = paste0(Anno.path, "/1.Umap.Cluster.pdf"), plot = p, height = 7, width = 7)
            saveRDS(scRNA, paste0("./data/seurat.annotated.rds"))
            scRNA.markers <- make_marker_table_sparse(scRNA, celltype = "cellType_1")
            write.csv(scRNA.markers, paste0(Anno.path, "/Findall.markers.cellType_1.csv"), row.names = FALSE)
            marker_features <- unique(head(scRNA.markers$gene, 80))
            marker_features <- intersect(marker_features, rownames(scRNA))
            if (length(marker_features)) {
                p <- DotPlot(scRNA, features = unique(marker_features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 
                  10)) + theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, 
                  vjust = 0.5), legend.position = "top", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, 
                  "cm"))
                ggsave(filename = paste0(Anno.path, "/3.Dotplot.pdf"), plot = p, width = 3 + length(unique(marker_features)) * 
                  0.3, height = 3 + length(unique(scRNA$cellType_1)) * 0.5)
            }
            ratio_sample <- as.data.frame(table(sample = scRNA$orig.ident, cellType = scRNA$cellType_1), stringsAsFactors = FALSE)
            write.table(ratio_sample, paste0(Anno.path, "/4.Sample.ratio.table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            if ("group" %in% colnames(scRNA@meta.data)) {
                ratio_group <- as.data.frame(table(group = scRNA$group, cellType = scRNA$cellType_1), stringsAsFactors = FALSE)
                write.table(ratio_group, paste0(Anno.path, "/5.Group.ratio.table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            }
            audit <- data.frame(item = c("annotation_mode", "cells", "celltypes", "singleR_dense_matrix"), value = c("sparse_marker_score_per_cell", 
                ncol(scRNA), length(unique(scRNA$cellType_1)), "skipped_for_fullcell_scalability"), check.names = FALSE)
            write.table(audit, paste0(Anno.path, "/0.scalable_annotation_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            return(invisible(scRNA))
        }
        p <- DimPlot(scRNA, group.by = Idents, label = T, reduction = "umap", raster = F) + NoLegend()
        ggsave(filename = paste0(Anno.path, "/1.Umap.Cluster.pdf"), plot = p, height = 7, width = 7)
        CELL <- SingleR::SingleR(test = as.matrix(scRNA@assays$RNA@data), ref = ref, labels = labels, clusters = scRNA@meta.data[[Idents]])
        cells <- CELL[, ncol(CELL):1]
        cells <- data.frame(cells)
        cell_id <- cells$pruned.labels
        clean_celltype <- function(x) {
            x %>% gsub("[\\_\\-]", " ", .) %>% gsub("\\s+", " ", .)
        }
        cell_id <- clean_celltype(cell_id)
        names(cell_id) <- levels(scRNA@meta.data[[Idents]])
        Idents(scRNA) <- Idents
        scRNA <- RenameIdents(scRNA, cell_id)
        scRNA$cellType_1 <- Idents(scRNA)
        print(cell_id)
        cell_id[which(is.na(cell_id))] <- "unknow"
        if (!anyNA(cell_id)) {
            names(cell_id) <- levels(scRNA@meta.data[[Idents]])
            Idents(scRNA) <- Idents
            scRNA <- RenameIdents(scRNA, cell_id)
            scRNA$cellType_1 <- Idents(scRNA)
            print(table(scRNA$cellType_1))
            saveRDS(scRNA, paste0("./data/seurat.annotated.rds"))
            Idents(scRNA) <- "cellType_1"
            scRNA.markers <- FindAllMarkers(scRNA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
            scRNA.markers <- scRNA.markers[scRNA.markers$p_val_adj < 0.05, ]
            write.csv(scRNA.markers, paste0(Anno.path, "/Findall.markers.cellType_1.csv"))
            features <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
            features <- features$gene
            my_pal2 <- c("#D4477D", "#4DBBD5", "#D24B27", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", 
                "#9ebcda", "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#00B3F1", "#00538A")
            umap <- scRNA@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cellType = scRNA@meta.data$cellType_1)
            celltypepos <- umap %>% group_by(cellType) %>% summarise(umap_1 = median(UMAP_1), umap_2 = median(UMAP_2))
            p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T) + 
                ggrepel::geom_label_repel(aes(x = umap_1, y = umap_2, label = cellType, color = cellType), fontface = "bold", 
                  data = celltypepos, box.padding = 0.5, show.legend = FALSE)
            ggsave(filename = paste0(Anno.path, "/2.Umap1.plot.pdf"), plot = p, width = 10, height = 7)
            p <- DotPlot(scRNA, features = unique(features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 
                10)) + theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, 
                vjust = 0.5), legend.position = "top", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, 
                "cm"), )
            ggsave(filename = paste0(Anno.path, "/3.Dotplot.pdf"), plot = p, width = 3 + length(unique(features)) * 0.3, 
                height = 3 + length(unique(scRNA$cellType_1)) * 0.5)
            cellRatio <- cellRatioPlot(object = scRNA, sample.name = "orig.ident", celltype.name = "cellType_1", flow.curve = 0.5, 
                fill.col = my_pal2) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
            ggsave(filename = paste0(Anno.path, "/4.Sample.ratio.plot.pdf"), plot = cellRatio, width = 3 + length(unique(scRNA$orig.ident)) * 
                0.5, height = 7)
            if ("group" %in% colnames(scRNA@meta.data)) {
                cellRatio <- cellRatioPlot(object = scRNA, sample.name = "group", celltype.name = "cellType_1", flow.curve = 0.5, 
                  fill.col = my_pal2) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
                ggsave(filename = paste0(Anno.path, "/5.Group.ratio.plot.pdf"), plot = cellRatio, width = 7, height = 7)
            }
        }
    }
    GeneExpre <- function(outpath = "./result/", SeuratObject = NA, KeyGene = NA) {
        d <- list.files(outpath, pattern = "KeyGene expression abundance", full.names = T)
        if (length(d) == 0) {
            Expr.path <- paste0(outpath, "/KeyGene expression abundance/")
            dir.create(Expr.path)
        }
        else {
            Expr.path <- d[1]
        }
        scRNA <- SeuratObject
        KeyGene <- intersect(rownames(scRNA), KeyGene)
        celltype_col <- if ("cellType_2" %in% colnames(scRNA@meta.data)) {
            "cellType_2"
        }
        else if ("cellType_1" %in% colnames(scRNA@meta.data)) {
            "cellType_1"
        }
        else {
            "seurat_clusters"
        }
        leg <- length(unique(KeyGene))
        p <- featureCornerAxes(object = scRNA, reduction = "umap", groupFacet = NULL, relLength = 0.5, relDist = 0.2, features = unique(KeyGene), 
            nLayout = 3, aspect.ratio = 1, pSize = 0.1)
        if (leg <= 3) {
            ggsave(filename = paste0(Expr.path, "/1.Keygene.FeaturePlot.pdf"), plot = p, width = leg * 5, height = 5, limitsize = FALSE)
        }
        if (leg > 3) {
            ggsave(filename = paste0(Expr.path, "/1.Keygene.FeaturePlot.pdf"), plot = p, width = 15, height = ceiling(leg/3) * 
                5, limitsize = FALSE)
        }
        p <- DotPlot(scRNA, group.by = celltype_col, features = unique(KeyGene), cols = "RdYlBu", ) + scale_size_continuous(range = c(0, 
            10)) + theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, 
            vjust = 0.5), legend.position = "right", text = element_text(size = 10))
        if (leg <= 3) {
            ggsave(filename = paste0(Expr.path, "/2.Keygene.DotPlot.pdf"), plot = p, width = 5 + ceiling(leg/3), height = 7, 
                limitsize = FALSE)
        }
        if (leg > 3) {
            ggsave(filename = paste0(Expr.path, "/2.Keygene.DotPlot.pdf"), plot = p, width = 5 + 0.5 * leg, height = 7, limitsize = FALSE)
        }
    }
    hallmark_module <- function(outpath = "./result/", SeuratObject = NA, KeyGene = NA) {
        e <- list.files(outpath, pattern = "Immunometabolic pathways", full.names = T)
        if (length(e) == 0) {
            hallmark.path <- paste0(outpath, "/Immunometabolic pathways/")
            dir.create(hallmark.path)
        }
        else {
            hallmark.path <- e[1]
        }
        scRNA <- SeuratObject
        ShowGen <- intersect(rownames(scRNA), KeyGene)
        Idents(scRNA) <- "cellType_2"
        load("./data/hallmark.gs.rdata")
        genesetInfo <- read.delim("./data/GenesetInfo.txt", sep = ",")
        geneset <- gs
        rownames(genesetInfo) <- genesetInfo$geneset
        genesetInfo <- genesetInfo[names(geneset), ]
        genesetInfo <- subset(genesetInfo, Classification != "")
        levels <- c("Immune", "Metabolism", "Signaling", "Proliferation")
        genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
        genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
        genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
        mat <- GetAssayData(object = scRNA, assay = "RNA", slot = "data")
        cells_rankings <- AUCell::AUCell_buildRankings(mat, nCores = 1, plotStats = F)
        score <- AUCell::AUCell_calcAUC(geneset, cells_rankings, nCores = 1, aucMaxRank = nrow(cells_rankings) * 0.05)
        score <- AUCell::getAUC(score)
        DP.list <- base::lapply(ShowGen, function(cellType) {
            DatGroup <- FetchData(scRNA, vars = c(cellType), slot = "data")
            group <- setNames(object = ifelse(DatGroup[, 1] > median(DatGroup[, 1]), "Hexp", "Lexp"), nm = rownames(DatGroup))
            if (length(unique(group)) > 1) {
                design <- model.matrix(~0 + factor(group))
                colnames(design) <- levels(factor(group))
                rownames(design) <- names(group)
                contrast.matrix <- limma::makeContrasts("Hexp-Lexp", levels = design)
                fit <- limma::lmFit(score[, names(group)], design)
                fit2 <- limma::contrasts.fit(fit, contrast.matrix)
                fit2 <- limma::eBayes(fit2)
                DPs <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = "bonferroni")
                DPs$KeyGenes <- cellType
                DPs$Pathway <- rownames(DPs)
                return(DPs)
            }
        })
        DPs <- do.call(rbind, DP.list)
        write.table(DPs, file = paste0(hallmark.path, "/output_hallmark.txt"), sep = "\t", row.names = F, col.names = T, 
            quote = F)
        plot.data <- DPs
        plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
        plot.data$KeyGenes <- factor(plot.data$KeyGenes)
        plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset)
        plot.data$FDR <- cut(plot.data$adj.P.Val, breaks = c(0, 1e-125, 1e-75, 1e-25, 1), include.lowest = T)
        plot.data$FDR <- factor(as.character(plot.data$FDR), levels = rev(levels(plot.data$FDR)))
        levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))
        color <- c("#4682B4", "#FFFFFF", "#CD2626")
        class.color <- c(Immune = "#D58986", Metabolism = "#80554C", Signaling = "#71AC7A", Proliferation = "#E8D4B4")
        p1 <- ggplot(plot.data, aes(x = Pathway, y = KeyGenes, color = logFC, size = FDR)) + geom_point() + scale_color_gradient2(low = color[1], 
            mid = color[2], high = color[3]) + geom_hline(yintercept = seq(min(as.numeric(plot.data$KeyGenes)) - 0.5, max(as.numeric(plot.data$KeyGenes)) + 
            0.5), color = "grey80") + geom_vline(xintercept = seq(min(as.numeric(plot.data$Pathway)) - 0.5, max(as.numeric(plot.data$Pathway)) + 
            0.5), color = "grey80") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
            axis.line = element_blank(), legend.position = "top", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, 
                "cm"), )
        p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) + geom_tile() + theme_classic() + scale_fill_manual(values = class.color) + 
            theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
                legend.position = "bottom")
        p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(10, 2))
        ggsave(filename = paste0(hallmark.path, "/1.Hallmark.pdf"), width = 15, height = 9, plot = p)
    }
    require(ClusterGVis)
    require(dplyr)
    require(jjAnno)
    require(Seurat)
    plan("multicore", workers = 8)
    options(future.globals.maxSize = 5e+06 * 1024^2)
    set.seed(123456)
    ClusterGVis_module <- function(object, OutPath, scRNA.markers = NULL, celltype = "cellType_1", species = NA, OnlyPlot = F) {
        cluster_dir <- file.path(OutPath, "ClusterGvis")
        if (!dir.exists(cluster_dir)) 
            dir.create(cluster_dir, recursive = TRUE)
        if (!OnlyPlot) {
            if (!is.na(species)) {
                mouse <- ifelse(species == "mouse", T, F)
            }
            else {
                mouse <- F
            }
            Idents(object) <- celltype
            if (is.null(scRNA.markers)) {
                scRNA.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
                write.csv(scRNA.markers, file.path(cluster_dir, paste0("Findall.markers.", celltype, ".csv")))
            }
            else {
                scRNA.markers <- scRNA.markers
            }
            scRNA.top.markers <- scRNA.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 5, wt = avg_log2FC)
            write.csv(scRNA.top.markers, file.path(cluster_dir, "TopMarkers.csv"))
            st.data <- prepareDataFromscRNA(object = object, diffData = scRNA.top.markers, showAverage = T)
            if (mouse) {
                require(org.Mm.eg.db)
                enrich <- enrichCluster(object = st.data, OrgDb = org.Mm.eg.db, type = "BP", organism = "mmu", pvalueCutoff = 0.05, 
                  topn = 5, seed = 123456)
            }
            else {
                require(org.Hs.eg.db)
                enrich <- enrichCluster(object = st.data, OrgDb = org.Hs.eg.db, type = "BP", organism = "hsa", pvalueCutoff = 0.05, 
                  topn = 5, seed = 123456)
            }
            saveRDS(st.data, file.path(cluster_dir, "Heatmap.plotdata.rds"))
            write.csv(enrich, file.path(cluster_dir, "Enrich.plotdata.csv"))
        }
        st.data <- readRDS(file.path(cluster_dir, "Heatmap.plotdata.rds"))
        enrich <- read.csv(file.path(cluster_dir, "Enrich.plotdata.csv"), row.names = 1)
        scRNA.top.markers <- read.csv(file.path(cluster_dir, "TopMarkers.csv"))
        cluster_n <- length(unique(st.data$wide.res$cluster))
        cluster_order <- seq_len(cluster_n)
        pdf(file.path(cluster_dir, "1.Heatmap.pdf"), height = 15, width = 17, onefile = F)
        visCluster(object = st.data, plot.type = "both", column_title_rot = 45, markGenes = unique(scRNA.top.markers$gene), 
            markGenes.side = "left", annoTerm.data = enrich, line.side = "left", cluster.order = cluster_order, add.bar = T)
        dev.off()
    }
    require(future)
    require(GSEABase)
    require(UCell)
    require(irGSEA)
    require(GSVA)
    require(Seurat)
    require(ggplot2)
    require(ggpubr)
    require(gghalves)
    require(ggrepel)
    require(RColorBrewer)
    require(patchwork)
    require(dplyr)
    require(tidyr)
    large_fullcell_mode_signature <- function(object) {
        threshold <- as.integer(Sys.getenv("SCALABLE_FULLCELL_THRESHOLD", "300000"))
        "no_sampling" %in% colnames(object@meta.data) && ncol(object) > threshold
    }
    normalize01 <- function(x) {
        rng <- range(x, na.rm = TRUE)
        if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) 
            return(rep(0, length(x)))
        (x - rng[1])/(rng[2] - rng[1])
    }
    score_signature_sparse <- function(object, signature, audit_outpath) {
        sig_name <- names(signature)[1]
        genes <- intersect(signature[[1]], rownames(object))
        mat <- GetAssayData(object = object, assay = "RNA", slot = "data")
        counts <- GetAssayData(object = object, assay = "RNA", slot = "counts")
        mean_expr <- Matrix::colMeans(mat[genes, , drop = FALSE])
        detect_frac <- Matrix::colMeans(counts[genes, , drop = FALSE] > 0)
        sum_expr <- Matrix::colSums(mat[genes, , drop = FALSE])
        score <- data.frame(AUCell = normalize01(detect_frac), UCell = normalize01(mean_expr), singscore = normalize01(sum_expr), 
            ssGSEA = normalize01(mean_expr), AddModuleScore = normalize01(mean_expr), check.names = FALSE)
        score$Scoring <- rowSums(score)
        rownames(score) <- colnames(object)
        score_cols <- colnames(score)
        existing_score_cols <- intersect(score_cols, colnames(object@meta.data))
        if (length(existing_score_cols) > 0) {
            object@meta.data[existing_score_cols] <- NULL
        }
        object@meta.data <- cbind(object@meta.data, score)
        audit <- data.frame(item = c("signature", "mode", "present_genes", "requested_genes", "cells", "dense_matrix"), value = c(sig_name, 
            "sparse_fullcell_no_sampling", length(genes), length(signature[[1]]), ncol(object), "not_used"), check.names = FALSE)
        write.table(audit, file.path(audit_outpath, "0.scalable_signature_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
        object
    }
    signature_module <- function(object, outpath, signature, group_by = "all", celltype = "cellType_1", OnlyPlot = F) {
        x <- list.files(outpath, pattern = "Signature", full.names = T)
        if (length(x) == 0) {
            outpath <- paste0(outpath, "/Signature/")
            dir.create(outpath)
        }
        else {
            outpath <- x
        }
        if (!OnlyPlot) {
            if (large_fullcell_mode_signature(object)) {
                object <- score_signature_sparse(object, signature, outpath)
            }
            else {
                mat <- GetAssayData(object = object, assay = "RNA")
                cells_rankings <- AUCell::AUCell_buildRankings(mat, nCores = 4, plotStats = F)
                score <- AUCell::AUCell_calcAUC(signature, cells_rankings, nCores = 4, aucMaxRank = nrow(cells_rankings) * 
                  0.05)
                aucs <- as.numeric(AUCell::getAUC(score)[names(signature), ])
                object@assays$UCell <- NULL
                object@assays$singscore <- NULL
                object <- irGSEA.score(object = object, assay = "RNA", slot = "data", seeds = 5201314, ncores = 4, msigdb = F, 
                  minGSSize = 0, maxGSSize = 1e+05, ucell.MaxRank = 1e+05, custom = T, geneset = signature, method = c("UCell", 
                    "singscore"), kcdf = "Gaussian")
                uc <- as.data.frame(object@assays$UCell@counts)
                uc <- as.data.frame(t(uc))
                singscore <- as.data.frame(object@assays$singscore@counts)
                singscore <- as.data.frame(t(singscore))
                exp <- as.matrix(object@assays$RNA@data)
                matrix <- gsva(exp, signature, kcdf = "Gaussian", method = "ssgsea", abs.ranking = T)
                ssgsea <- as.data.frame(t(matrix))
                object <- AddModuleScore(object, features = signature, name = "Add")
                score <- data.frame(AUCell = aucs, UCell = uc[[names(signature)]], singscore = singscore[[names(signature)]], 
                  ssGSEA = ssgsea[[names(signature)]], AddModuleScore = object$Add1)
                score <- scale(score)
                normalize <- function(x) {
                  return((x - min(x))/(max(x) - min(x)))
                }
                score <- apply(score, 2, normalize)
                score <- as.data.frame(score)
                score$Scoring <- rowSums(score)
                colnames(object@meta.data)
                object$Add1 <- NULL
                score_cols <- colnames(score)
                existing_score_cols <- intersect(score_cols, colnames(object@meta.data))
                if (length(existing_score_cols) > 0) {
                  object@meta.data[existing_score_cols] <- NULL
                }
                object@meta.data <- cbind(object@meta.data, score)
            }
            option <- list(signature = names(signature), size = length(signature[[1]]), group_by = group_by)
            save(option, file = paste0(outpath, "/option.rdata"))
        }
        p1 <- DotPlot(object, cols = "RdYlBu", group.by = celltype, features = c("AUCell", "UCell", "singscore", "ssGSEA", 
            "AddModuleScore", "Scoring")) + theme_bw() + scale_size_continuous(range = c(0, 8)) + theme(axis.text.x = element_text(angle = 45, 
            size = 12, hjust = 1, vjust = 1), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 15), 
            axis.title.x = element_text(size = 15)) + geom_point(shape = 21, color = "black", stroke = 0.5, mapping = aes_string(size = "pct.exp")) + 
            xlab("Methods") + ylab("Cell Type") + coord_flip()
        ggsave(paste0(outpath, "/1.Dotplot of ", names(signature), ".pdf"), plot = p1, width = 8, height = 5.5)
        p2 <- FeaturePlot(object = object, ncol = 1, features = "Scoring", pt.size = 0.1 * 20000/ncol(object)) + viridis::scale_color_viridis() + 
            labs(title = names(signature)) + tidydr::theme_dr() + theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15))
        if (group_by == "all") {
            tmp <- object@meta.data[, c(celltype, "Scoring")]
            colnames(tmp) <- c("cellType_1", "val")
            top_cell <- tmp %>% group_by(cellType_1) %>% summarise(mean_score = mean(val, na.rm = TRUE)) %>% arrange(desc(mean_score)) %>% 
                slice(1)
            key_cell <- data.frame(key_cell = top_cell$cellType_1)
            write.table(key_cell, "./data/key_cell.txt", sep = "\t", row.names = FALSE, quote = FALSE)
            p3 <- ggplot(data = tmp, aes(x = cellType_1, y = val, fill = cellType_1, color = cellType_1)) + geom_boxplot(width = 0.2, 
                fill = "transparent", size = 0.4, outlier.shape = NA) + geom_violin(color = NA, alpha = 0.35) + theme_classic() + 
                scale_fill_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"))) + 
                scale_color_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"))) + 
                theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(vjust = 0.5, 
                  size = 12), axis.title.y = element_text(angle = 90, size = 15), axis.title.x = element_text(size = 15)) + 
                xlab("Cell Type") + ylab(paste0(names(signature), " Score"))
        }
        else if (group_by != "all") {
            tmp <- object@meta.data[, c(celltype, "Scoring", group_by)]
            colnames(tmp) <- c("cellType_1", "val", "group")
            tmp$cellType_1 <- as.factor(tmp$cellType_1)
            p3 <- ggplot(data = tmp, aes(x = cellType_1, y = val, fill = group, color = group)) + geom_boxplot(show.legend = T, 
                width = 0.6, alpha = 0.5, outlier.color = "grey50") + scale_fill_manual(values = c("#398AB9", "red")) + scale_color_manual(values = c("#398AB9", 
                "red")) + theme_classic() + stat_compare_means(method = "wilcox.test", hjust = 0.5, vjust = 0, hide.ns = T, 
                label = "p.signif", show.legend = F) + theme(legend.position = "top", axis.text.x = element_text(angle = 45, 
                hjust = 1, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12), axis.title.y = element_text(angle = 90, 
                size = 15), axis.title.x = element_text(size = 15)) + xlab("Cell Type") + ylab(paste0(names(signature), " Score"))
            diff_median <- tmp %>% group_by(cellType_1, group) %>% summarise(median_score = median(val, na.rm = TRUE), .groups = "drop") %>% 
                pivot_wider(names_from = group, values_from = median_score) %>% rowwise() %>% mutate(diff = {
                vals <- c_across(where(is.numeric))
                if (length(vals) == 2) 
                  diff(vals)
                else NA_real_
            }) %>% ungroup() %>% na.omit()
            key_cell <- data.frame(key_cell = diff_median[which(diff_median$diff == max(diff_median$diff)), ]$cellType_1)
            write.table(key_cell, "./data/key_cell.txt", sep = "\t", row.names = FALSE, quote = FALSE)
        }
        p4 <- CombinePlots(plots = list(p2, p3), rel_widths = c(1, 1.5))
        ggsave(paste0(outpath, "/2.Boxplot of ", names(signature), ".pdf"), plot = p4, width = 12.5, height = 5)
        return(object)
    }
    signature_volcano <- function(object, outpath, ShowCell = NA, group_by = "ScoreGroup", celltype = "cellType_1", idents = c("Hsco", 
        "Lsco"), logfc.threshold = 0.25, OnlyPlot = F) {
        x <- list.files(outpath, pattern = "Signature", full.names = T)
        if (length(x) == 0) {
            outpath <- paste0(outpath, "/Signature/")
            dir.create(outpath)
        }
        else {
            outpath <- x
        }
        if (!OnlyPlot) {
            object <- object[, object@meta.data[[celltype]] %in% ShowCell]
            if (group_by == "ScoreGroup") {
                object$ScoreGroup <- ifelse(object$Scoring > median(object$Scoring), "Hsco", "Lsco")
            }
            Idents(object) <- group_by
            if (large_fullcell_mode_signature(object)) {
                mat <- GetAssayData(object = object, assay = "RNA", slot = "data")
                counts <- GetAssayData(object = object, assay = "RNA", slot = "counts")
                g <- as.character(object@meta.data[[group_by]])
                idx1 <- g == idents[1]
                idx2 <- g == idents[2]
                mean1 <- Matrix::rowMeans(mat[, idx1, drop = FALSE])
                mean2 <- Matrix::rowMeans(mat[, idx2, drop = FALSE])
                pct1 <- Matrix::rowMeans(counts[, idx1, drop = FALSE] > 0)
                pct2 <- Matrix::rowMeans(counts[, idx2, drop = FALSE] > 0)
                DEGs <- data.frame(p_val = 1e-300, avg_log2FC = as.numeric(mean1 - mean2), pct.1 = as.numeric(pct1), pct.2 = as.numeric(pct2), 
                  p_val_adj = 1e-300, row.names = rownames(mat), check.names = FALSE)
                audit <- data.frame(item = c("de_mode", "ident_1_cells", "ident_2_cells", "features"), value = c("sparse_fullcell_mean_difference_no_sampling", 
                  sum(idx1), sum(idx2), nrow(mat)), check.names = FALSE)
                write.table(audit, file.path(outpath, "0.scalable_sparse_de_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            }
            else {
                DEGs <- FindMarkers(object, assay = "RNA", ident.1 = idents[1], ident.2 = idents[2], logfc.threshold = 0, 
                  min.pct = 0.1)
            }
            sig.DEGs <- DEGs[DEGs$p_val_adj < 0.05 & abs(DEGs$avg_log2FC) > logfc.threshold, ]
            write.csv(DEGs, file = paste0(outpath, "/", idents[1], "~", idents[2], "-AllDEGs.csv"), row.names = T, quote = F)
            write.csv(sig.DEGs, file = paste0(outpath, "/", idents[1], "~", idents[2], "-LogFC", logfc.threshold, ".csv"), 
                row.names = T, quote = F)
            load(paste0(outpath, "/option.rdata"))
            option$ShowCell <- ShowCell
            option$idents <- idents
            option$LogFC <- logfc.threshold
            option$numDEGs <- nrow(sig.DEGs)
            save(option, file = paste0(outpath, "/option.rdata"))
        }
        DEGs <- read.csv(paste0(outpath, "/", idents[1], "~", idents[2], "-AllDEGs.csv"), header = T, row.names = 1, check.names = F)
        logFCfilter <- logfc.threshold
        p_valFilter <- 0.05
        DEGs <- DEGs[-log10(DEGs$p_val_adj) > 0, ]
        DEGs$gene <- rownames(DEGs)
        DEGs$lab <- ""
        DEGs[order(DEGs$avg_log2FC), ][c(1:10, (nrow(DEGs) - 9):nrow(DEGs)), ]$lab <- DEGs[order(DEGs$avg_log2FC), ][c(1:10, 
            (nrow(DEGs) - 9):nrow(DEGs)), ]$gene
        DEGs$logP <- -log10(DEGs$p_val_adj)
        p5 <- ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC)) + geom_point(aes(size = logP), 
            alpha = 0.9) + scale_color_gradientn(colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"), name = "logFC", 
            values = seq(0, 1, 0.2)) + scale_fill_gradientn(colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"), 
            name = "logFC", values = seq(0, 1, 0.2)) + scale_size("-log10(p)") + geom_text_repel(data = DEGs, aes(x = avg_log2FC, 
            y = -log10(p_val_adj), label = lab), size = 4, box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"), 
            segment.color = "black", show.legend = F) + theme_bw() + ylab("-log10 (p_val_adj)") + xlab("avg_log2FC") + geom_vline(xintercept = c(-logFCfilter, 
            logFCfilter), lty = 3, col = "black", lwd = 0.5) + geom_hline(yintercept = -log10(p_valFilter), lty = 3, col = "black", 
            lwd = 0.5) + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 15), 
            axis.title.y = element_text(size = 15))
        ggsave(filename = paste0(outpath, "/3.Volcano diagram of ", option$signature, ".pdf"), plot = p5, width = 10, height = 8)
    }
    require(Seurat)
    require(monocle)
    require(dplyr)
    require(ggplot2)
    require(scales)
    require(ClusterGVis)
    set.seed(5201314)
    my_pal2 <- c("#df4a86", "#746ea3", "#009ecb", "#00827b", "#3d4d7c", "#ad341d", "#0056a0", "#d77b1c", "#077140", "#d34132", 
        "#698da5", "#ffd597", "#f1977f", "#828bae", "#82bfbd", "#c1000a", "#FF8E00", "#00B3F1", "#354270", "#85b38f")
    monocle2_module <- function(object, OutPath, group_by, track_gene = NULL, ShowGene = NA, OnlyPlot = F) {
        monocle_dir <- file.path(OutPath, "Monocle2")
        if (!dir.exists(monocle_dir)) 
            dir.create(monocle_dir, recursive = TRUE)
        if (!OnlyPlot) {
            scRNASub <- object
            monocle.matrix <- GetAssayData(object = scRNASub, slot = "data", assay = "RNA")
            monocle.sample <- scRNASub@meta.data
            monocle.geneAnn <- data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
            data <- as(as.matrix(monocle.matrix), "sparseMatrix")
            pd <- new("AnnotatedDataFrame", data = monocle.sample)
            fd <- new("AnnotatedDataFrame", data = monocle.geneAnn)
            cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
            names(pData(cds))[names(pData(cds)) == group_by] <- "Cluster"
            cds <- estimateSizeFactors(cds)
            cds <- estimateDispersions(cds)
            DefaultAssay(scRNASub) <- "RNA"
            if (is.null(track_gene)) {
                track_gene = c("VariableFeatures", "all")
                genenubmers <- "all"
            }
            else if (length(track_gene) == 1) {
                genenubmers <- "all"
                if (!is.na(as.numeric(as.character(track_gene[1])))) {
                  genenubmers <- as.numeric(as.character(track_gene[1]))
                  track_gene <- c("VariableFeatures", genenubmers)
                }
                else if (track_gene %in% c("dispersion", "cluster", "VariableFeatures")) {
                  track_gene <- c(track_gene, 2000)
                }
            }
            else if (length(track_gene) == 2 && !is.na(as.numeric(as.character(track_gene[2])))) {
                genenubmers <- as.numeric(as.character(track_gene[2]))
                if (!track_gene[1] %in% c("dispersion", "cluster", "VariableFeatures")) {
                  track_gene <- c("VariableFeatures", 2000)
                }
            }
            if (length(track_gene) >= 200) {
                genenubmers <- length(track_gene)
                track_gene <- track_gene
            }
            else if (track_gene[1] == "dispersion") {
                disp_table <- dispersionTable(cds)
                disp_table <- arrange(disp_table, -dispersion_empirical)
                track_gene <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
            }
            else if (track_gene[1] == "cluster") {
                Idents(scRNASub) <- "seurat_clusters"
                deg.cluster <- FindAllMarkers(scRNASub)
                deg.cluster <- arrange(deg.cluster, p_val)
                track_gene <- subset(deg.cluster, p_val_adj < 0.05)$gene
            }
            else if (track_gene[1] == "VariableFeatures") {
                track_gene <- VariableFeatures(scRNASub)
            }
            if (length(track_gene) <= genenubmers | genenubmers == "all") {
                track_gene <- track_gene
            }
            else if (length(track_gene) >= genenubmers) {
                track_gene <- track_gene[1:genenubmers]
            }
            cds <- setOrderingFilter(cds, track_gene)
            cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
            cds <- orderCells(cds)
            saveRDS(cds, file.path(monocle_dir, "cds.rds"))
        }
        cds <- readRDS(file.path(monocle_dir, "cds.rds"))
        p <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.1) + ggsci::scale_color_gsea() + theme(legend.position = "right", 
            text = element_text(face = "bold"))
        ggsave(file.path(monocle_dir, "1.Trajectory.Pseudotime.pdf"), plot = p, height = 4, width = 4.5)
        p <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.1) + guides(color = guide_legend(override.aes = list(alpha = 1, 
            size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right", text = element_text(face = "bold"))
        ggsave(file.path(monocle_dir, "2.Trajectory.State.pdf"), plot = p, height = 4, width = 4.5)
        p <- plot_cell_trajectory(cds, color_by = "Cluster", cell_size = 0.1) + guides(color = guide_legend(override.aes = list(alpha = 1, 
            size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right", text = element_text(face = "bold"))
        ggsave(file.path(monocle_dir, "3.Trajectory.Cluster.pdf"), plot = p, height = 4, width = 4.5)
        if (!class(ShowGene) == "logical") {
            cds_subset <- cds[ShowGene, ]
            p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Cluster") + scale_color_manual(values = my_pal2)
            p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "State") + scale_color_manual(values = my_pal2)
            p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime") + ggsci::scale_color_gsea()
            p1 | p2 | p3
            ggsave(filename = file.path(monocle_dir, "4.ShowGene.pseudotime.pdf"), plot = p1 | p2 | p3, width = 12, height = length(ShowGene) * 
                2)
        }
        return(cds)
    }
    monocle2_branch_module <- function(cds, OutPath = "./", branch_point = 1, genes_branched = "all", num_clusters = 3, heatmap_genes = 100, 
        species = NA, ncores = 16, OnlyPlot = F) {
        monocle_dir <- file.path(OutPath, "Monocle2")
        if (!dir.exists(monocle_dir)) 
            dir.create(monocle_dir, recursive = TRUE)
        object <- cds
        if (class(object) == "CellDataSet") {
            cds <- object
        }
        else if (is.character(object) && file.exists(object) && paste0(unlist(strsplit(object, split = ""))[c((nchar(object) - 
            2):nchar(object))], collapse = "") == "rds") {
            cds <- readRDS(file = object)
        }
        if (!OnlyPlot) {
            if (genes_branched == "all") {
                genes_branched <- rownames(cds)
            }
            else if (is.numeric(genes_branched)) {
                track_gene <- rownames(cds@featureData@data)[cds@featureData@data$use_for_ordering]
                genes_branched <- track_gene[1:genes_branched]
            }
            gene_cal <- BEAM(cds[genes_branched, ], branch_point = branch_point, cores = ncores, progenitor_method = "duplicate")
            gene_cal <- gene_cal[order(gene_cal$qval), ]
            write.csv(gene_cal, file = file.path(monocle_dir, paste0("heatmap_branch_point_", branch_point, "_gene.csv")))
            gene_cal <- read.csv(file.path(monocle_dir, paste0("heatmap_branch_point_", branch_point, "_gene.csv")), row.names = 1)
            if (is.numeric(heatmap_genes)) {
                heatmap_genes <- rownames(gene_cal[1:heatmap_genes, ])
            }
            else if (!all(heatmap_genes %in% rownames(cds))) {
                heatmap_genes <- rownames(gene_cal)[1:100]
            }
            else if (all(heatmap_genes %in% rownames(cds))) {
                heatmap_genes <- heatmap_genes
            }
            len <- length(heatmap_genes)
            df <- plot_genes_branched_heatmap2(cds[heatmap_genes, ], branch_point = branch_point, num_clusters = num_clusters, 
                cores = ncores, use_gene_short_name = T, show_rownames = T)
            saveRDS(df, file.path(monocle_dir, "Branch.plotdata.rds"))
            if (!is.na(species)) {
                mouse <- ifelse(species == "mouse", T, F)
            }
            else {
                mouse <- F
            }
            if (mouse) {
                require(org.Mm.eg.db)
                enrich <- enrichCluster(object = df, OrgDb = org.Mm.eg.db, type = "KEGG", organism = "mmu", pvalueCutoff = 0.05, 
                  topn = 5, seed = 5201314)
            }
            else {
                require(org.Hs.eg.db)
                enrich <- enrichCluster(object = df, OrgDb = org.Hs.eg.db, type = "KEGG", organism = "hsa", pvalueCutoff = 0.05, 
                  topn = 5, add.gene = T, seed = 5201314)
            }
            write.csv(enrich, file.path(monocle_dir, "Enrich.plotdata.csv"))
        }
        df <- readRDS(file.path(monocle_dir, "Branch.plotdata.rds"))
        enrich <- read.csv(file.path(monocle_dir, "Enrich.plotdata.csv"), row.names = 1)
        markGenes <- intersect(df$wide.res$gene, do.call(c, strsplit(enrich$geneID, "/", F)))
        enrich$geneID <- NULL
        pdf(file.path(monocle_dir, "5.Branch-enrich.pdf"), height = 7, width = 13, onefile = F)
        visCluster(object = df, plot.type = "both", column_names_rot = 45, show_row_dend = F, markGenes = markGenes, markGenes.side = "left", 
            annoTerm.data = enrich, add.bar = T, line.side = "left")
        dev.off()
        if (!is.null(dev.list())) {
            dev.off()
        }
    }
    require(data.table)
    require(igraph, lib.loc = "/usr/local/lib/R-4.3.0/lib/R/library/")
    require(CellChat)
    require(dplyr)
    require(Seurat)
    require(patchwork)
    require(ggsci)
    require(RColorBrewer)
    require(ggpubr)
    require(scales)
    filter_cellchat_db_to_present_features <- function(db, present_features, out_dir = NULL) {
        present_features <- unique(as.character(present_features))
        lr <- db$interaction
        if (is.null(lr) || nrow(lr) == 0) {
            return(db)
        }
        normalize_symbol_vector <- function(x) {
            x <- as.character(x)
            x <- x[!is.na(x) & nzchar(x)]
            x[!tolower(x) %in% c("na", "none", "null")]
        }
        expand_db_entity <- function(x, ref_df, candidate_cols) {
            x <- normalize_symbol_vector(x)
            if (!length(x)) {
                return(character())
            }
            unlist(lapply(x, function(v) {
                if (!is.null(ref_df) && v %in% rownames(ref_df) && length(candidate_cols) > 0) {
                  vals <- unlist(ref_df[v, candidate_cols, drop = FALSE], use.names = FALSE)
                  return(normalize_symbol_vector(vals))
                }
                v
            }), use.names = FALSE)
        }
        get_lr_field <- function(field, i) {
            if (field %in% colnames(lr)) {
                lr[[field]][i]
            }
            else {
                NA_character_
            }
        }
        complex_df <- db$complex
        complex_cols <- if (!is.null(complex_df)) 
            grep("^subunit", colnames(complex_df), value = TRUE)
        else character()
        cofactor_df <- db$cofactor
        cofactor_cols <- if (!is.null(cofactor_df)) 
            grep("^cofactor|^subunit", colnames(cofactor_df), value = TRUE)
        else character()
        cofactor_fields <- c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")
        needed_by_interaction <- lapply(seq_len(nrow(lr)), function(i) {
            unique(c(expand_db_entity(lr$ligand[i], complex_df, complex_cols), expand_db_entity(lr$receptor[i], complex_df, 
                complex_cols), unlist(lapply(cofactor_fields, function(field) {
                expand_db_entity(get_lr_field(field, i), cofactor_df, cofactor_cols)
            }), use.names = FALSE)))
        })
        missing_by_interaction <- lapply(needed_by_interaction, function(features) {
            setdiff(features, present_features)
        })
        keep <- vapply(missing_by_interaction, length, integer(1)) == 0
        if (!is.null(out_dir) && nzchar(out_dir)) {
            if (!dir.exists(out_dir)) 
                dir.create(out_dir, recursive = TRUE)
            write.table(data.frame(item = c("mode", "input_interactions", "retained_interactions", "removed_interactions", 
                "present_features", "filter_scope"), value = c("present_feature_interaction_filter", nrow(lr), sum(keep), 
                sum(!keep), length(present_features), "ligand_receptor_complex_and_cofactor_features"), check.names = FALSE), 
                file.path(out_dir, "0.cellchat_db_feature_filter.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            removed_idx <- which(!keep)
            if (length(removed_idx) > 0) {
                removed <- data.frame(interaction_name = if ("interaction_name" %in% colnames(lr)) 
                  lr$interaction_name[removed_idx]
                else rownames(lr)[removed_idx], ligand = lr$ligand[removed_idx], receptor = lr$receptor[removed_idx], missing_features = vapply(missing_by_interaction[removed_idx], 
                  paste, character(1), collapse = ","), check.names = FALSE)
                write.table(removed, file.path(out_dir, "0.cellchat_db_feature_filter_removed.tsv"), sep = "\t", row.names = FALSE, 
                  quote = FALSE)
            }
        }
        db$interaction <- lr[keep, , drop = FALSE]
        db
    }
    cellchat_func <- function(object, species = NA, celltype = "cellType_1", ShowCell = NA, OnlyPlot = F, AllPlot = F, out_dir) {
        if (!dir.exists(out_dir)) 
            dir.create(out_dir, recursive = TRUE)
        if (!OnlyPlot) {
            if (!is.na(species)) {
                mouse <- ifelse(species == "mouse", T, F)
            }
            else {
                mouse <- F
            }
            if (mouse) {
                CellChatDB <- CellChatDB.mouse
                PPI.use <- PPI.mouse
            }
            else {
                CellChatDB <- CellChatDB.human
                PPI.use <- PPI.human
            }
            large_fullcell <- "no_sampling" %in% colnames(object@meta.data) && ncol(object) > as.integer(Sys.getenv("SCALABLE_FULLCELL_THRESHOLD", 
                "300000"))
            min_cells_filter <- 10
            if (large_fullcell) {
                mat <- object@assays$RNA@data
                groups <- as.character(object@meta.data[[celltype]])
                groups[is.na(groups) | !nzchar(groups)] <- "Unknown"
                group_levels <- sort(unique(groups))
                avg <- sapply(group_levels, function(g) Matrix::rowMeans(mat[, groups == g, drop = FALSE]))
                colnames(avg) <- group_levels
                data.input <- as.matrix(avg)
                identity <- data.frame(labels = group_levels, row.names = group_levels)
                colnames(identity) <- celltype
                min_cells_filter <- 0
                write.table(data.frame(celltype = names(table(groups)), full_cell_count = as.integer(table(groups)), check.names = FALSE), 
                  file.path(out_dir, "0.fullcell_celltype_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
                write.table(data.frame(item = c("mode", "source_cells", "celltypes", "cell_sampling", "filterCommunication_min_cells"), 
                  value = c("aggregated_fullcell_mean_expression", ncol(object), length(group_levels), "none", min_cells_filter), 
                  check.names = FALSE), file.path(out_dir, "0.scalable_cellchat_audit.tsv"), sep = "\t", row.names = FALSE, 
                  quote = FALSE)
            }
            else {
                data.input <- object@assays$RNA@data
                identity <- data.frame(labels = as.character(object@meta.data[[celltype]]), row.names = colnames(data.input))
                colnames(identity) <- celltype
            }
            cellchat <- createCellChat(data.input, meta = identity, group.by = celltype)
            CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor"))
            CellChatDB.use <- filter_cellchat_db_to_present_features(CellChatDB.use, rownames(data.input), out_dir)
            cellchat@DB <- CellChatDB.use
            cellchat <- subsetData(cellchat)
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            overexpressed_lr_n <- if (is.null(cellchat@LR$LRsig)) 
                0
            else nrow(cellchat@LR$LRsig)
            measurable_lr_fallback <- FALSE
            if (overexpressed_lr_n == 0 && large_fullcell) {
                cellchat@LR$LRsig <- cellchat@DB$interaction
                measurable_lr_fallback <- TRUE
            }
            lr_for_probability_n <- if (is.null(cellchat@LR$LRsig)) 
                0
            else nrow(cellchat@LR$LRsig)
            write.table(data.frame(item = c("overexpressed_lr_interactions", "lr_for_computeCommunProb", "measurable_lr_fallback"), 
                value = c(overexpressed_lr_n, lr_for_probability_n, measurable_lr_fallback), check.names = FALSE), file.path(out_dir, 
                "0.cellchat_lr_selection_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            cellchat <- projectData(cellchat, PPI.use)
            cellchat <- computeCommunProb(cellchat)
            cellchat <- filterCommunication(cellchat, min.cells = min_cells_filter)
            df.net <- subsetCommunication(cellchat)
            write.table(df.net, file = file.path(out_dir, "net_lr.txt"), quote = F, sep = "\t", row.names = F)
            cellchat <- computeCommunProbPathway(cellchat)
            df.netp <- subsetCommunication(cellchat)
            write.table(df.netp, file = file.path(out_dir, "net_pathway.txt"), quote = F, sep = "\t", row.names = F)
            cellchat <- aggregateNet(cellchat)
            cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
            pathways.shows <- cellchat@netP$pathways
            head(cellchat@LR$LRsig)
            saveRDS(cellchat, file = file.path(out_dir, "cellchat.rds"))
        }
        cellchat <- readRDS(file.path(out_dir, "cellchat.rds"))
        groupSize <- as.numeric(table(cellchat@idents))
        pathways.shows <- cellchat@netP$pathways
        {
            fig1 <- file.path(out_dir, "Net_number_circle_barplot")
            dir.create(fig1, recursive = T)
            pdf(file = paste0(fig1, "/Net_number_strength.pdf"), width = 10, height = 5)
            par(mfrow = c(1, 2))
            netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
            netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
            dev.off()
            file.copy(paste0(fig1, "/Net_number_strength.pdf"), file.path(out_dir, "1.Net_number_strength.pdf"))
            df.net <- read.table(file.path(out_dir, "net_lr.txt"), sep = "\t", check.names = F, header = T)
            data <- as.data.frame(table(c(df.net$source, df.net$target)))
            colnames(data) <- c("Cell_Type", "all_sum")
            data <- data[order(data$all_sum, decreasing = T), ]
            data$Cell_Type <- factor(data$Cell_Type, levels = data$Cell_Type)
            head(data)
            sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
            sample_color2 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], 
                ifelse(data$all_sum > summary(data$all_sum)[5], sample_color[1], sample_color[2]), sample_color[3]), sample_color[4])
            sample_color1 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], 
                ifelse(data$all_sum > summary(data$all_sum)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), 
                alpha(sample_color[3], 0.9)), alpha(sample_color[4], 0.9))
            pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) + geom_bar(stat = "identity", 
                width = 0.8) + scale_fill_manual(values = sample_color1) + scale_colour_manual(values = sample_color2) + 
                theme_bw() + theme(panel.grid = element_blank()) + labs(x = "Cell Type", y = "Counts") + theme(axis.text.y = element_text(size = 12, 
                colour = "black")) + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
            ggsave(paste0(fig1, "/Interaction Count.pdf"), plot = pB2, width = 10, height = 8)
            file.copy(paste0(fig1, "/Interaction Count.pdf"), file.path(out_dir, "2.Interaction Count.pdf"))
            if (!is.na(ShowCell)) {
                ident1 <- which(levels(cellchat@idents) %in% ShowCell)
                ident2 <- which(!levels(cellchat@idents) %in% ShowCell)
                gg1 <- netVisual_bubble(cellchat, sources.use = ident1, targets.use = ident2, angle.x = 45, remove.isolate = F, 
                  font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling from ", ShowCell))
                gg2 <- netVisual_bubble(cellchat, sources.use = ident2, targets.use = ident1, angle.x = 45, remove.isolate = F, 
                  font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling to ", ShowCell))
                ggsave(paste0(fig1, "/Bubble.pdf"), plot = gg1$gg.obj + gg2$gg.obj, width = 5 + length(levels(cellchat@idents)) * 
                  1, height = 8 + max(nrow(gg1$communication), nrow(gg2$communication)) * 0.05)
                file.copy(paste0(fig1, "/Bubble.pdf"), file.path(out_dir, "3.Bubble.pdf"))
            }
        }
        if (AllPlot) {
            {
                fig2 <- file.path(out_dir, "NetVisual_circle")
                dir.create(fig2, recursive = T)
                mat <- cellchat@net$weight
                for (i in 1:nrow(mat)) {
                  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
                  mat2[i, ] <- mat[i, ]
                  pdf(file = paste0(fig2, "/netVisual_weight_", rownames(mat)[i], ".pdf"), width = 10, height = 9)
                  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
                  dev.off()
                }
            }
            {
                fig3 <- file.path(out_dir, "NetVisual_aggregate")
                dir.create(fig3, recursive = T)
                for (pathways.show in pathways.shows) {
                  pdf(file = paste0(fig3, "/netVisual_aggregate_", pathways.show, ".pdf"), width = 10, height = 9)
                  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", pt.title = 50, vertex.label.cex = 1)
                  dev.off()
                }
            }
            {
                fig4 <- file.path(out_dir, "NetVisual_heatmap")
                dir.create(fig4, recursive = T)
                for (pathways.show in pathways.shows) {
                  pdf(file = paste0(fig4, "/netVisual_heatmap_", pathways.show, ".pdf"), width = 10, height = 9)
                  print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", font.size = 14, font.size.title = 20))
                  dev.off()
                }
            }
            {
                fig5 <- file.path(out_dir, "NetAnalysis_contribution")
                dir.create(fig5, recursive = T)
                for (pathways.show in pathways.shows) {
                  p <- netAnalysis_contribution(cellchat, signaling = pathways.show, font.size = 15, font.size.title = 15, 
                    title = paste0("Contribution of each L-R pair in ", pathways.show)) + theme(text = element_text(face = "bold"))
                  ggsave(filename = paste0(fig5, "/netAnalysis_contribution_", pathways.show, ".pdf"), plot = p)
                }
            }
            {
                fig6 <- file.path(out_dir, "NetVisual_bubble")
                dir.create(fig6, recursive = T)
                for (sources.use in levels(cellchat@idents)) {
                  targets.use <- levels(cellchat@idents)
                  p <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, angle.x = 45, remove.isolate = F, 
                    font.size = 10, font.size.title = 15, return.data = T, title.name = paste0("Signaling from ", sources.use))
                  ggsave(filename = paste0(fig6, "/netVisual_bubble_from_", sources.use, ".pdf"), plot = p$gg.obj, width = 8, 
                    height = 8 + nrow(p$communication) * 0.05)
                }
                for (targets.use in levels(cellchat@idents)) {
                  sources.use <- levels(cellchat@idents)
                  p <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, angle.x = 45, remove.isolate = F, 
                    font.size = 10, font.size.title = 15, return.data = T, title.name = paste0("Signaling to ", targets.use))
                  ggsave(filename = paste0(fig6, "/netVisual_bubble_to_", targets.use, ".pdf"), plot = p$gg.obj, width = 8, 
                    height = 8 + nrow(p$communication) * 0.05)
                }
            }
            {
                fig7 <- file.path(out_dir, "NetAnalysis_signalingRole")
                dir.create(fig7, recursive = T)
                gg1 <- netAnalysis_signalingRole_scatter(cellchat)
                ggsave(filename = paste0(fig7, "/netAnalysis_signalingRole_scatter", ".pdf"), plot = gg1, width = 8, height = 8)
                ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 7 + length(cellchat@netP$pathways) * 
                  0.125)
                ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 7 + length(cellchat@netP$pathways) * 
                  0.125)
                pdf(paste0(fig7, "/netAnalysis_signalingRole_heatmap.pdf"), width = 15, height = 7 + length(cellchat@netP$pathways) * 
                  0.125)
                print(ht1 + ht2)
                dev.off()
            }
            {
                fig8 <- file.path(out_dir, "PlotGeneExpression")
                dir.create(fig8, recursive = T)
                for (pathways.show in pathways.shows) {
                  res <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = TRUE, enriched.only = T)
                  p <- plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = TRUE, type = "violin")
                  ggsave(filename = paste0(fig8, "/PlotGeneExpression_", pathways.show, ".pdf"), plot = p, width = 10, height = 7 + 
                    length(res$geneLR) * 0.225)
                }
            }
        }
        return(cellchat)
    }
    cellchat_module <- function(object, OutPath, group_by = "all", celltype = "cellType_1", idents = NULL, ShowCell = NA, 
        species = NA, OnlyPlot = F, AllPlot = F) {
        if (group_by == "all" & is.null(idents)) {
            out_path <- file.path(OutPath, "CellChat")
            if (!dir.exists(out_path)) 
                dir.create(out_path, recursive = TRUE)
            cellchat <- cellchat_func(object = object, species = species, celltype = celltype, ShowCell = ShowCell, OnlyPlot = OnlyPlot, 
                AllPlot = AllPlot, out_dir = out_path)
        }
        else if (group_by != "all" & length(idents) == 2) {
            base_cellchat_dir <- file.path(OutPath, "CellChat")
            if (!dir.exists(base_cellchat_dir)) 
                dir.create(base_cellchat_dir, recursive = TRUE)
            for (ident in idents) {
                out_path <- file.path(base_cellchat_dir, ident)
                if (!dir.exists(out_path)) 
                  dir.create(out_path, recursive = TRUE)
                sub.object <- object[, object@meta.data[[group_by]] %in% ident]
                cellchat <- cellchat_func(object = sub.object, species = species, celltype = celltype, ShowCell = ShowCell, 
                  OnlyPlot = OnlyPlot, AllPlot = AllPlot, out_dir = out_path)
            }
            cellchat.NL <- readRDS(file.path(base_cellchat_dir, idents[1], "cellchat.rds"))
            cellchat.TL <- readRDS(file.path(base_cellchat_dir, idents[2], "cellchat.rds"))
            object.list <- setNames(list(cellchat.NL, cellchat.TL), idents[1:2])
            cellchat <- mergeCellChat(object.list, add.names = names(object.list))
            save(cellchat, file = file.path(base_cellchat_dir, paste0("cellchat_merged_", idents[1], "_", idents[2], ".RData")))
            gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
            gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
            ggsave(file.path(base_cellchat_dir, "1.CompareInteractions.pdf"), plot = gg1 + gg2, width = 5, height = 5)
            pdf(file.path(base_cellchat_dir, "2.NetVisual.pdf"), width = 10, height = 5)
            par(mfrow = c(1, 2), xpd = TRUE)
            netVisual_diffInteraction(cellchat, weight.scale = T)
            netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
            dev.off()
            pdf(file.path(base_cellchat_dir, "3.NetVisualHeatmap.pdf"), width = 10, height = 5)
            gg1 <- netVisual_heatmap(cellchat)
            gg2 <- netVisual_heatmap(cellchat, measure = "weight")
            print(gg1 + gg2)
            dev.off()
            weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
            par(mfrow = c(1, 2), xpd = TRUE)
            for (i in 1:length(object.list)) {
                netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2], 
                  edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
            }
            num.link <- sapply(object.list, function(x) {
                rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
            })
            weight.MinMax <- c(min(num.link), max(num.link))
            gg <- list()
            for (i in 1:length(object.list)) {
                gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
            }
            ggsave(file.path(base_cellchat_dir, "4.SignalingRole.pdf"), plot = patchwork::wrap_plots(plots = gg), width = 10, 
                height = 5)
            if (!is.na(ShowCell)) {
                ident1 <- which(levels(object@meta.data[[celltype]]) %in% ShowCell)
                ident2 <- which(!levels(object@meta.data[[celltype]]) %in% ShowCell)
                gg1 <- netVisual_bubble(cellchat, sources.use = ident1, targets.use = ident2, comparison = c(1, 2), angle.x = 45, 
                  remove.isolate = F, font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling from ", 
                    ShowCell))
                gg2 <- netVisual_bubble(cellchat, sources.use = ident2, targets.use = ident1, comparison = c(1, 2), angle.x = 45, 
                  remove.isolate = F, font.size = 12, font.size.title = 15, return.data = T, title.name = paste0("Signaling to ", 
                    ShowCell))
                ggsave(file.path(base_cellchat_dir, "5.Bubble.pdf"), plot = gg1$gg.obj + gg2$gg.obj, width = 5 + length(levels(object@meta.data[[celltype]])) * 
                  1, height = 8 + max(nrow(gg1$communication), nrow(gg2$communication)) * 0.05)
            }
        }
    }
    RunML <- function(method, Train_set, Train_label, mode = "Model", classVar) {
        method = gsub(" ", "", method)
        method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
        method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)
        method_param = switch(EXPR = method_name, Enet = list(alpha = as.numeric(gsub("alpha=", "", method_param))), Stepglm = list(direction = method_param), 
            NULL)
        args = list(Train_set = Train_set, Train_label = Train_label, mode = mode, classVar = classVar)
        args = c(args, method_param)
        obj <- do.call(what = paste0("Run", method_name), args = args)
        return(obj)
    }
    RunEnet <- function(Train_set, Train_label, mode, classVar, alpha) {
        cv.fit = cv.glmnet(x = Train_set, y = Train_label[[classVar]], family = "binomial", alpha = alpha, nfolds = 10)
        fit = glmnet(x = Train_set, y = Train_label[[classVar]], family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunLasso <- function(Train_set, Train_label, mode, classVar) {
        RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
    }
    RunRidge <- function(Train_set, Train_label, mode, classVar) {
        RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
    }
    RunStepglm <- function(Train_set, Train_label, mode, classVar, direction) {
        fit <- step(glm(formula = Train_label[[classVar]] ~ ., family = "binomial", data = as.data.frame(Train_set)), direction = direction, 
            trace = 0)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunSVM <- function(Train_set, Train_label, mode, classVar) {
        data <- as.data.frame(Train_set)
        data[[classVar]] <- as.factor(Train_label[[classVar]])
        fit = svm(formula = eval(parse(text = paste(classVar, "~."))), data = data, probability = T)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunLDA <- function(Train_set, Train_label, mode, classVar) {
        data <- as.data.frame(Train_set)
        data[[classVar]] <- as.factor(Train_label[[classVar]])
        fit = train(eval(parse(text = paste(classVar, "~."))), data = data, method = "lda", trControl = trainControl(method = "cv"))
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunglmBoost <- function(Train_set, Train_label, mode, classVar) {
        data <- cbind(Train_set, Train_label[classVar])
        data[[classVar]] <- as.factor(data[[classVar]])
        fit <- glmboost(eval(parse(text = paste(classVar, "~."))), data = data, family = Binomial())
        cvm <- cvrisk(fit, papply = lapply, folds = cv(model.weights(fit), type = "kfold"))
        fit <- glmboost(eval(parse(text = paste(classVar, "~."))), data = data, family = Binomial(), control = boost_control(mstop = max(mstop(cvm), 
            40)))
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunplsRglm <- function(Train_set, Train_label, mode, classVar) {
        cv.plsRglm.res = cv.plsRglm(formula = Train_label[[classVar]] ~ ., data = as.data.frame(Train_set), nt = 10, verbose = FALSE)
        fit <- plsRglm(Train_label[[classVar]], as.data.frame(Train_set), modele = "pls-glm-logistic", verbose = F, sparse = T)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunRF <- function(Train_set, Train_label, mode, classVar) {
        rf_nodesize = 5
        Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
        fit <- rfsrc(formula = formula(paste0(classVar, "~.")), data = cbind(Train_set, Train_label[classVar]), ntree = 1000, 
            nodesize = rf_nodesize, importance = T, proximity = T, forest = T)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunGBM <- function(Train_set, Train_label, mode, classVar) {
        fit <- gbm(formula = Train_label[[classVar]] ~ ., data = as.data.frame(Train_set), distribution = "bernoulli", n.trees = 10000, 
            interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, cv.folds = 10, n.cores = 6)
        best <- which.min(fit$cv.error)
        fit <- gbm(formula = Train_label[[classVar]] ~ ., data = as.data.frame(Train_set), distribution = "bernoulli", n.trees = best, 
            interaction.depth = 3, n.minobsinnode = 10, shrinkage = 0.001, n.cores = 8)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunXGBoost <- function(Train_set, Train_label, mode, classVar) {
        indexes = createFolds(Train_label[[classVar]], k = 5, list = T)
        CV <- unlist(lapply(indexes, function(pt) {
            dtrain = xgb.DMatrix(data = Train_set[-pt, ], label = Train_label[-pt, ])
            dtest = xgb.DMatrix(data = Train_set[pt, ], label = Train_label[pt, ])
            watchlist <- list(train = dtrain, test = dtest)
            bst <- xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 2, nrounds = 10, watchlist = watchlist, objective = "binary:logistic", 
                verbose = F)
            which.min(bst$evaluation_log$test_logloss)
        }))
        nround <- as.numeric(names(which.max(table(CV))))
        fit <- xgboost(data = Train_set, label = Train_label[[classVar]], max.depth = 2, eta = 1, nthread = 2, nrounds = nround, 
            objective = "binary:logistic", verbose = F)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    RunNaiveBayes <- function(Train_set, Train_label, mode, classVar) {
        data <- cbind(Train_set, Train_label[classVar])
        data[[classVar]] <- as.factor(data[[classVar]])
        fit <- naiveBayes(eval(parse(text = paste(classVar, "~."))), data = data)
        fit$subFeature = colnames(Train_set)
        if (mode == "Model") 
            return(fit)
        if (mode == "Variable") 
            return(ExtractVar(fit))
    }
    quiet <- function(..., messages = FALSE, cat = FALSE) {
        if (!cat) {
            sink(tempfile())
            on.exit(sink())
        }
        out <- if (messages) 
            eval(...)
        else suppressMessages(eval(...))
        out
    }
    standarize.fun <- function(indata, centerFlag, scaleFlag) {
        scale(indata, center = centerFlag, scale = scaleFlag)
    }
    scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL) {
        samplename = rownames(data)
        if (is.null(cohort)) {
            data <- list(data)
            names(data) = "training"
        }
        else {
            data <- split(as.data.frame(data), cohort)
        }
        if (is.null(centerFlags)) {
            centerFlags = F
        }
        if (length(centerFlags) == 1) {
            centerFlags = rep(centerFlags, length(data))
        }
        if (is.null(names(centerFlags))) {
            names(centerFlags) <- names(data)
        }
        if (is.null(scaleFlags)) {
            scaleFlags = F
        }
        if (length(scaleFlags) == 1) {
            scaleFlags = rep(scaleFlags, length(data))
        }
        if (is.null(names(scaleFlags))) {
            names(scaleFlags) <- names(data)
        }
        centerFlags <- centerFlags[names(data)]
        scaleFlags <- scaleFlags[names(data)]
        outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
        outdata <- do.call(rbind, outdata)
        outdata <- outdata[samplename, ]
        return(outdata)
    }
    ExtractVar <- function(fit) {
        Feature <- quiet(switch(EXPR = class(fit)[1], lognet = rownames(coef(fit))[which(coef(fit)[, 1] != 0)], glm = names(coef(fit)), 
            svm.formula = fit$subFeature, train = fit$coefnames, glmboost = names(coef(fit)[abs(coef(fit)) > 0]), plsRglmmodel = rownames(fit$Coeffs)[fit$Coeffs != 
                0], rfsrc = var.select(fit, verbose = F)$topvars, gbm = rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, 
                plotit = F)$rel.inf > 0], xgb.Booster = fit$subFeature, naiveBayes = fit$subFeature))
        Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
        return(Feature)
    }
    CalPredictScore <- function(fit, new_data, type = "lp") {
        new_data <- new_data[, fit$subFeature]
        RS <- quiet(switch(EXPR = class(fit)[1], lognet = predict(fit, type = "response", as.matrix(new_data)), glm = predict(fit, 
            type = "response", as.data.frame(new_data)), svm.formula = predict(fit, as.data.frame(new_data), probability = T), 
            train = predict(fit, new_data, type = "prob")[[2]], glmboost = predict(fit, type = "response", as.data.frame(new_data)), 
            plsRglmmodel = predict(fit, type = "response", as.data.frame(new_data)), rfsrc = predict(fit, as.data.frame(new_data))$predicted[, 
                "1"], gbm = predict(fit, type = "response", as.data.frame(new_data)), xgb.Booster = predict(fit, as.matrix(new_data)), 
            naiveBayes = predict(object = fit, type = "raw", newdata = new_data)[, "1"]))
        RS = as.numeric(as.vector(RS))
        names(RS) = rownames(new_data)
        return(RS)
    }
    PredictClass <- function(fit, new_data) {
        new_data <- new_data[, fit$subFeature]
        label <- quiet(switch(EXPR = class(fit)[1], lognet = predict(fit, type = "class", as.matrix(new_data)), glm = ifelse(test = predict(fit, 
            type = "response", as.data.frame(new_data)) > 0.5, yes = "1", no = "0"), svm.formula = predict(fit, as.data.frame(new_data), 
            decision.values = T), train = predict(fit, new_data, type = "raw"), glmboost = predict(fit, type = "class", as.data.frame(new_data)), 
            plsRglmmodel = ifelse(test = predict(fit, type = "response", as.data.frame(new_data)) > 0.5, yes = "1", no = "0"), 
            rfsrc = predict(fit, as.data.frame(new_data))$class, gbm = ifelse(test = predict(fit, type = "response", as.data.frame(new_data)) > 
                0.5, yes = "1", no = "0"), xgb.Booster = ifelse(test = predict(fit, as.matrix(new_data)) > 0.5, yes = "1", 
                no = "0"), naiveBayes = predict(object = fit, type = "class", newdata = new_data)))
        label = as.character(as.vector(label))
        names(label) = rownames(new_data)
        return(label)
    }
    RunEval <- function(fit, Test_set = NULL, Test_label = NULL, Train_set = NULL, Train_label = NULL, Train_name = NULL, 
        cohortVar = "Cohort", classVar) {
        if ((!is.null(Train_set)) & (!is.null(Train_label))) {
            new_data <- rbind.data.frame(Train_set[, fit$subFeature], Test_set[, fit$subFeature])
            if (!is.null(Train_name)) {
                Train_label$Cohort <- Train_name
            }
            else {
                Train_label$Cohort <- "Training"
            }
            colnames(Train_label)[ncol(Train_label)] <- cohortVar
            Test_label <- rbind.data.frame(Train_label[, c(cohortVar, classVar)], Test_label[, c(cohortVar, classVar)])
            Test_label[, 1] <- factor(Test_label[, 1], levels = c(unique(Train_label[, cohortVar]), setdiff(unique(Test_label[, 
                cohortVar]), unique(Train_label[, cohortVar]))))
        }
        else {
            new_data <- Test_set[, fit$subFeature]
        }
        RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
        Predict.out <- Test_label
        Predict.out$RS <- as.vector(RS)
        Predict.out <- split(x = Predict.out, f = Predict.out[, cohortVar])
        unlist(lapply(Predict.out, function(data) {
            as.numeric(auc(suppressMessages(roc(data[[classVar]], data$RS))))
        }))
    }
    SimpleHeatmap <- function(Cindex_mat, avg_Cindex, CohortCol, barCol, cellwidth = 1, cellheight = 0.5, cluster_columns, 
        cluster_rows) {
        col_ha = columnAnnotation(Cohort = colnames(Cindex_mat), col = list(Cohort = CohortCol), show_annotation_name = F)
        row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE, gp = gpar(fill = barCol, col = NA), 
            add_numbers = T, numbers_offset = unit(-10, "mm"), axis_param = list(labels_rot = 0), numbers_gp = gpar(fontsize = 9, 
                col = "white"), width = unit(3, "cm")), show_annotation_name = F)
        Heatmap(as.matrix(Cindex_mat), name = "AUC", right_annotation = row_ha, top_annotation = col_ha, col = c("#4195C1", 
            "#FFFFFF", "#CB5746"), rect_gp = gpar(col = "black", lwd = 1), cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
            show_column_names = FALSE, show_row_names = TRUE, row_names_side = "left", width = unit(cellwidth * ncol(Cindex_mat) + 
                2, "cm"), height = unit(cellheight * nrow(Cindex_mat), "cm"), column_split = factor(colnames(Cindex_mat), 
                levels = colnames(Cindex_mat)), column_title = NULL, cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3), x, y, gp = gpar(fontsize = 10))
            })
    }
    suppressPackageStartupMessages({
        library(data.table)
        library(limma)
        library(openxlsx)
        library(seqinr)
        library(plyr)
        library(randomForestSRC)
        library(glmnet)
        library(plsRglm)
        library(gbm)
        library(caret)
        library(mboost)
        library(e1071)
        library(BART)
        library(MASS)
        library(snowfall)
        library(xgboost)
        library(ComplexHeatmap)
        library(RColorBrewer)
        library(pROC)
        library(dplyr)
        library(UpSetR)
    })
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) 
        a
    else b
    pipeline_root <- function() normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    source_113ml_reference <- function() {
        patch_113ml_predictors()
        invisible("embedded")
    }
    patch_113ml_predictors <- function() {
        CalPredictScore <<- function(fit, new_data, type = "lp") {
            new_data <- as.data.frame(new_data[, fit$subFeature, drop = FALSE])
            RS <- quiet(switch(EXPR = class(fit)[1], lognet = predict(fit, type = "response", as.matrix(new_data)), glm = predict(fit, 
                type = "response", as.data.frame(new_data)), svm.formula = predict(fit, as.data.frame(new_data), probability = TRUE), 
                train = predict(fit, new_data, type = "prob")[[2]], glmboost = predict(fit, type = "response", as.data.frame(new_data)), 
                plsRglmmodel = predict(fit, type = "response", as.data.frame(new_data)), rfsrc = predict(fit, as.data.frame(new_data))$predicted[, 
                  "1"], gbm = predict(fit, type = "response", as.data.frame(new_data)), xgb.Booster = predict(fit, as.matrix(new_data)), 
                naiveBayes = predict(object = fit, type = "raw", newdata = new_data)[, "1"]))
            RS <- as.numeric(as.vector(RS))
            names(RS) <- rownames(new_data)
            RS
        }
        PredictClass <<- function(fit, new_data) {
            new_data <- as.data.frame(new_data[, fit$subFeature, drop = FALSE])
            label <- quiet(switch(EXPR = class(fit)[1], lognet = predict(fit, type = "class", as.matrix(new_data)), glm = ifelse(predict(fit, 
                type = "response", as.data.frame(new_data)) > 0.5, "1", "0"), svm.formula = predict(fit, as.data.frame(new_data), 
                decision.values = TRUE), train = predict(fit, new_data, type = "raw"), glmboost = predict(fit, type = "class", 
                as.data.frame(new_data)), plsRglmmodel = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 
                0.5, "1", "0"), rfsrc = predict(fit, as.data.frame(new_data))$class, gbm = ifelse(predict(fit, type = "response", 
                as.data.frame(new_data)) > 0.5, "1", "0"), xgb.Booster = ifelse(predict(fit, as.matrix(new_data)) > 0.5, 
                "1", "0"), naiveBayes = predict(object = fit, type = "class", newdata = new_data)))
            label <- as.character(as.vector(label))
            names(label) <- rownames(new_data)
            label
        }
        RunEval <<- function(fit, Test_set = NULL, Test_label = NULL, Train_set = NULL, Train_label = NULL, Train_name = NULL, 
            cohortVar = "Cohort", classVar) {
            if ((!is.null(Train_set)) && (!is.null(Train_label))) {
                new_data <- rbind.data.frame(as.data.frame(Train_set[, fit$subFeature, drop = FALSE]), as.data.frame(Test_set[, 
                  fit$subFeature, drop = FALSE]))
                if (!is.null(Train_name)) {
                  Train_label$Cohort <- Train_name
                }
                else {
                  Train_label$Cohort <- "Training"
                }
                colnames(Train_label)[ncol(Train_label)] <- cohortVar
                Test_label <- rbind.data.frame(Train_label[, c(cohortVar, classVar), drop = FALSE], Test_label[, c(cohortVar, 
                  classVar), drop = FALSE])
                Test_label[, 1] <- factor(Test_label[, 1], levels = c(unique(Train_label[, cohortVar]), setdiff(unique(Test_label[, 
                  cohortVar]), unique(Train_label[, cohortVar]))))
            }
            else {
                new_data <- as.data.frame(Test_set[, fit$subFeature, drop = FALSE])
            }
            RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
            Predict.out <- Test_label
            Predict.out$RS <- as.vector(RS)
            Predict.out <- split(x = Predict.out, f = Predict.out[, cohortVar])
            unlist(lapply(Predict.out, function(data) {
                as.numeric(auc(suppressMessages(roc(data[[classVar]], data$RS))))
            }))
        }
    }
    load_113_methods <- function(methods_file = NULL) {
        candidates <- c(methods_file, file.path(data_root, "113_ML_methods.txt"), file.path(data_root, "113_ML_methods.txt"))
        candidates <- candidates[nzchar(candidates %||% "")]
        methods_path <- candidates[file.exists(candidates)][1]
        methods <- read.table(methods_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)$x
        methods <- unique(trimws(methods[nzchar(methods)]))
        max_methods <- Sys.getenv("ML113_MAX_METHODS", unset = "")
        if (nzchar(max_methods)) {
            max_methods <- suppressWarnings(as.integer(max_methods))
            if (!is.na(max_methods) && max_methods > 0) 
                methods <- head(methods, max_methods)
        }
        methods
    }
    read_expr_matrix <- function(path) {
        raw <- data.table::fread(path, data.table = FALSE, header = TRUE, check.names = FALSE)
        gene <- as.character(raw[[1]])
        mat <- raw[, -1, drop = FALSE]
        mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))
        mat <- limma::avereps(as.matrix(mat), ID = gene)
        if (all(is.finite(mat)) && range(mat, na.rm = TRUE)[1] >= 0 && diff(range(mat, na.rm = TRUE)) > 100) {
            mat <- log2(mat + 1)
        }
        mat
    }
    cohort_name_from_expr <- function(path) {
        base <- basename(path)
        sub("-test-expr.*$", "", base)
    }
    read_metadata <- function(path, cohort = "Train") {
        tb <- read.table(path, header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
        sample_id <- as.character(tb[[1]])
        group_raw <- if (ncol(tb) >= 3) 
            tb[[3]]
        else tb[[ncol(tb)]]
        group_chr <- trimws(as.character(group_raw))
        group_num <- suppressWarnings(as.numeric(group_chr))
        if (all(!is.na(group_num))) {
            group <- as.integer(group_num)
        }
        else {
            group <- ifelse(tolower(group_chr) %in% c("control", "normal", "healthy", "hc", "0"), 0L, 1L)
        }
        out <- data.frame(Group = group, Cohort = cohort, row.names = sample_id, check.names = FALSE)
        out
    }
    load_expr_and_labels <- function(data_dir = "data") {
        expr_files <- c(file.path(data_dir, "expr.txt"), sort(list.files(data_dir, pattern = "test-expr", full.names = TRUE)))
        meta_files <- c(file.path(data_dir, "metadata.txt"), sort(list.files(data_dir, pattern = "test-metadata", full.names = TRUE)))
        cohort_names <- c("Train", vapply(expr_files[-1], cohort_name_from_expr, character(1)))
        expr_list <- Map(read_expr_matrix, expr_files)
        names(expr_list) <- cohort_names
        label_list <- Map(read_metadata, meta_files, cohort_names)
        names(label_list) <- cohort_names
        common_features <- Reduce(intersect, lapply(expr_list, rownames))
        expr_list <- lapply(expr_list, function(x) x[common_features, , drop = FALSE])
        list(expr = expr_list, labels = label_list)
    }
    normalize_gene_id <- function(x) {
        y <- trimws(as.character(x))
        y <- sub("\\|.*$", "", y)
        y <- sub("\\..*$", "", y)
        toupper(y)
    }
    match_candidate_genes <- function(candidate_genes, available_genes) {
        candidate_genes <- unique(trimws(as.character(candidate_genes)))
        candidate_genes <- candidate_genes[nzchar(candidate_genes) & !is.na(candidate_genes)]
        exact <- candidate_genes[candidate_genes %in% available_genes]
        norm_map <- setNames(available_genes, normalize_gene_id(available_genes))
        norm_hit <- unname(norm_map[normalize_gene_id(candidate_genes)])
        unique(c(exact, norm_hit[!is.na(norm_hit)]))
    }
    limit_candidate_genes <- function(genes, expr_list) {
        max_features <- Sys.getenv("ML113_MAX_FEATURES", unset = "120")
        max_features <- suppressWarnings(as.integer(max_features))
        if (is.na(max_features) || max_features <= 0 || length(genes) <= max_features) 
            return(genes)
        train_expr <- expr_list[[1]][genes, , drop = FALSE]
        gene_mad <- apply(train_expr, 1, stats::mad, na.rm = TRUE)
        gene_mad[!is.finite(gene_mad)] <- 0
        ranked <- names(sort(gene_mad, decreasing = TRUE))
        head(ranked, max_features)
    }
    build_ml_sets <- function(expr_list, label_list, candidate_genes) {
        genes <- match_candidate_genes(candidate_genes, rownames(expr_list[[1]]))
        genes <- limit_candidate_genes(genes, expr_list)
        expr_list <- lapply(expr_list, function(x) x[genes, , drop = FALSE])
        make_set <- function(mat, labels) {
            sample_ids <- intersect(colnames(mat), rownames(labels))
            x <- t(mat[, sample_ids, drop = FALSE])
            y <- labels[sample_ids, , drop = FALSE]
            list(x = x, y = y)
        }
        sets <- Map(make_set, expr_list, label_list)
        train <- sets[[1]]
        test_x <- do.call(rbind, lapply(sets[-1], `[[`, "x"))
        test_y <- do.call(rbind, lapply(sets[-1], `[[`, "y"))
        if (is.null(test_x) || nrow(test_x) == 0) {
            test_x <- train$x
            test_y <- train$y
            test_y$Cohort <- "Train_as_Test"
        }
        train_x <- as.matrix(train$x)
        test_x <- as.matrix(test_x)
        train_rows <- rownames(train_x)
        train_cols <- colnames(train_x)
        test_rows <- rownames(test_x)
        test_cols <- colnames(test_x)
        train_x <- scaleData(train_x, centerFlags = TRUE, scaleFlags = TRUE)
        test_x <- scaleData(test_x, cohort = test_y$Cohort, centerFlags = TRUE, scaleFlags = TRUE)
        train_x <- as.matrix(train_x)
        test_x <- as.matrix(test_x)
        rownames(train_x) <- train_rows
        colnames(train_x) <- train_cols
        rownames(test_x) <- test_rows
        colnames(test_x) <- test_cols
        list(train_x = train_x, train_y = train$y["Group"], test_x = test_x, test_y = test_y, genes = genes, expr_list = expr_list, 
            labels = label_list)
    }
    safe_run_ml <- function(method, train_x, train_y, mode, class_var, failures) {
        caller <- parent.frame()
        result <- tryCatch(list(ok = TRUE, value = RunML(method = method, Train_set = train_x, Train_label = train_y, mode = mode, 
            classVar = class_var)), error = function(e) list(ok = FALSE, message = conditionMessage(e)))
        if (!is.list(result) || !isTRUE(result$ok)) {
            current <- get("failures", envir = caller, inherits = FALSE)
            current[[length(current) + 1L]] <- data.frame(method = method, mode = mode, message = result$message, stringsAsFactors = FALSE)
            assign("failures", current, envir = caller)
            return(NULL)
        }
        result$value
    }
    write_failures <- function(failures, outdir) {
        if (length(failures) == 0) 
            return(invisible(NULL))
        fail_df <- do.call(rbind, failures)
        write.table(fail_df, file.path(outdir, "113_model_failures.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
    plot_pretrain_upset <- function(pretrain_var, outdir) {
        gene_lists <- pretrain_var[names(pretrain_var) != "simple"]
        gene_lists <- gene_lists[lengths(gene_lists) > 0]
        if (length(gene_lists) < 2) 
            return(invisible(NULL))
        pdf(file.path(outdir, "1.gene_intersection_upset.pdf"), width = 10, height = 6)
        print(UpSetR::upset(UpSetR::fromList(gene_lists), sets = names(gene_lists), order.by = "freq", text.scale = 1.2, 
            mainbar.y.label = "Gene number intersected", sets.x.label = "Gene number selected"))
        dev.off()
        core_genes <- Reduce(intersect, gene_lists)
        write.table(core_genes, file.path(outdir, "113_selected_core_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    plot_auc_heatmap <- function(auc_mat, outdir) {
        if (is.null(auc_mat) || nrow(auc_mat) == 0 || ncol(auc_mat) == 0) 
            return(invisible(NULL))
        avg_auc <- sort(apply(auc_mat, 1, mean, na.rm = TRUE), decreasing = TRUE)
        auc_mat <- auc_mat[names(avg_auc), , drop = FALSE]
        avg_auc <- as.numeric(format(avg_auc, digits = 3, nsmall = 3))
        n_col <- ncol(auc_mat)
        if (n_col <= 12) {
            cohort_col <- RColorBrewer::brewer.pal(max(3, n_col), "Paired")[seq_len(n_col)]
        }
        else {
            cohort_col <- (grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired")))(n_col)
        }
        names(cohort_col) <- colnames(auc_mat)
        cellwidth <- 1
        cellheight <- 0.5
        hm <- SimpleHeatmap(Cindex_mat = auc_mat, avg_Cindex = avg_auc, CohortCol = cohort_col, barCol = "steelblue", cellwidth = cellwidth, 
            cellheight = cellheight, cluster_columns = FALSE, cluster_rows = FALSE)
        pdf(file.path(outdir, "ML_AUC_heatmap.pdf"), width = cellwidth * ncol(auc_mat) + 6, height = max(4, cellheight * 
            nrow(auc_mat) * 0.45))
        draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
        dev.off()
        file.copy(file.path(outdir, "ML_AUC_heatmap.pdf"), file.path(outdir, "3.RandomForest.pdf"), overwrite = TRUE)
    }
    plot_model_roc <- function(risk_tab, labels, best_method, outdir) {
        score <- risk_tab[[best_method]]
        names(score) <- risk_tab$id
        labels <- labels[intersect(rownames(labels), names(score)), , drop = FALSE]
        for (cohort in unique(labels$Cohort)) {
            y <- labels[labels$Cohort == cohort, "Group"]
            ids <- rownames(labels)[labels$Cohort == cohort]
            if (length(unique(y)) < 2) 
                next
            roc1 <- pROC::roc(y, as.numeric(score[ids]), quiet = TRUE)
            ci_vec <- suppressWarnings(as.numeric(pROC::ci.auc(roc1, method = "bootstrap")))
            pdf(file.path(outdir, paste0("ROC.", cohort, ".pdf")), width = 5, height = 4.75)
            plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE, main = cohort)
            if (length(ci_vec) == 3) {
                text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ci_vec[1]), "-", sprintf("%.03f", ci_vec[3])), col = "red")
            }
            dev.off()
        }
    }
    plot_gene_roc <- function(expr_list, label_list, genes, outdir) {
        for (cohort in names(expr_list)) {
            mat <- expr_list[[cohort]]
            labels <- label_list[[cohort]]
            sample_ids <- intersect(colnames(mat), rownames(labels))
            y <- labels[sample_ids, "Group"]
            if (length(unique(y)) < 2) 
                next
            genes_use <- intersect(genes, rownames(mat))
            if (length(genes_use) == 0) 
                next
            bio_col <- rainbow(length(genes_use), s = 0.9, v = 0.9)
            filename <- if (cohort == "Train") 
                "Train_ROC.genes.pdf"
            else paste0("Test_", cohort, "_ROC.genes.pdf")
            pdf(file.path(outdir, filename), width = 5, height = 5)
            auc_text <- character(0)
            for (i in seq_along(genes_use)) {
                gene <- genes_use[i]
                roc1 <- pROC::roc(y, as.numeric(mat[gene, sample_ids]), quiet = TRUE)
                if (i == 1) {
                  plot(roc1, print.auc = FALSE, col = bio_col[i], legacy.axes = TRUE, main = cohort, lwd = 3)
                }
                else {
                  plot(roc1, print.auc = FALSE, col = bio_col[i], legacy.axes = TRUE, main = cohort, lwd = 3, add = TRUE)
                }
                auc_text <- c(auc_text, paste0(gene, ", AUC=", sprintf("%.3f", roc1$auc[1])))
            }
            legend("bottomright", auc_text, lwd = 3, bty = "n", cex = 0.75, col = bio_col, inset = c(0.05, 0))
            dev.off()
        }
    }
    write_compatibility_outputs <- function(outdir, sig_genes, best_method, auc_mat) {
        best_auc <- if (!is.null(auc_mat) && best_method %in% rownames(auc_mat)) 
            auc_mat[best_method, , drop = FALSE]
        else NULL
        sig <- data.frame(KeyGene = sig_genes, BestMethod = best_method, stringsAsFactors = FALSE)
        if (!is.null(best_auc)) {
            for (cn in colnames(best_auc)) sig[[cn]] <- as.numeric(best_auc[1, cn])
            sig$MeanAUC <- mean(as.numeric(best_auc[1, ]), na.rm = TRUE)
        }
        write.table(sig, file.path(outdir, "sigGenes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        write.table(data.frame(x = sig_genes), file.path(outdir, "best_modelGene.list.txt"), sep = "\t", row.names = FALSE, 
            quote = FALSE)
        write.csv(data.frame(gene = sig_genes), file.path(outdir, "feature_lasso.csv"), row.names = FALSE)
        write.table(data.frame(Gene = sig_genes, Importance = seq_along(sig_genes)), file.path(outdir, "importance of variables.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
    }
    run_113ml_pipeline <- function(candidate_genes, outdir, data_dir = "data", min_selected_var = 0, methods_file = NULL, 
        class_var = "Group") {
        source_113ml_reference()
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        data <- load_expr_and_labels(data_dir)
        ml <- build_ml_sets(data$expr, data$labels, candidate_genes)
        write.table(data.frame(Gene = ml$genes), file.path(outdir, "113_candidate_genes_used.txt"), sep = "\t", row.names = FALSE, 
            quote = FALSE)
        methods <- load_113_methods(methods_file)
        pretrain_methods <- c("Lasso", "glmBoost", "RF", "Stepglm[both]", "Stepglm[backward]")
        failures <- list()
        pretrain_var <- list()
        set.seed(3)
        for (method in pretrain_methods) {
            tm <- system.time({
                pretrain_var[[method]] <- safe_run_ml(method, ml$train_x, ml$train_y, "Variable", class_var, failures)
            })
        }
        pretrain_var <- pretrain_var[!vapply(pretrain_var, is.null, logical(1))]
        pretrain_var[["simple"]] <- colnames(ml$train_x)
        plot_pretrain_upset(pretrain_var, outdir)
        model <- list()
        set.seed(4)
        train_bk <- ml$train_x
        for (method_name in methods) {
            parts <- strsplit(method_name, "\\+")[[1]]
            if (length(parts) == 1) 
                parts <- c("simple", parts)
            if (!parts[1] %in% names(pretrain_var)) {
                failures[[length(failures) + 1]] <- data.frame(method = method_name, mode = "Model", message = paste0("pretrain method unavailable: ", 
                  parts[1]), stringsAsFactors = FALSE)
                next
            }
            vars <- intersect(pretrain_var[[parts[1]]], colnames(train_bk))
            if (length(vars) == 0) {
                failures[[length(failures) + 1]] <- data.frame(method = method_name, mode = "Model", message = "no variables retained", 
                  stringsAsFactors = FALSE)
                next
            }
            fit <- safe_run_ml(parts[2], train_bk[, vars, drop = FALSE], ml$train_y, "Model", class_var, failures)
            if (is.null(fit)) 
                next
            selected <- ExtractVar(fit)
            if (length(selected) <= min_selected_var) {
                failures[[length(failures) + 1]] <- data.frame(method = method_name, mode = "Model", message = paste0("selected gene count <= ", 
                  min_selected_var), stringsAsFactors = FALSE)
                next
            }
            model[[method_name]] <- fit
        }
        write_failures(failures, outdir)
        saveRDS(model, file.path(outdir, "final.ML_model.rds"))
        methods_valid <- names(model)
        combined_x <- rbind.data.frame(ml$train_x, ml$test_x)
        combined_labels <- rbind(data.frame(ml$train_y, Cohort = "Train", row.names = rownames(ml$train_y), check.names = FALSE), 
            ml$test_y[, c("Group", "Cohort"), drop = FALSE])
        rs_list <- lapply(methods_valid, function(method) CalPredictScore(model[[method]], combined_x))
        names(rs_list) <- methods_valid
        risk_tab <- as.data.frame(t(do.call(rbind, rs_list)), check.names = FALSE)
        risk_tab <- cbind(id = rownames(risk_tab), risk_tab)
        write.table(risk_tab, file.path(outdir, "Risk_matrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        class_list <- lapply(methods_valid, function(method) PredictClass(model[[method]], combined_x))
        names(class_list) <- methods_valid
        class_tab <- as.data.frame(t(do.call(rbind, class_list)), check.names = FALSE)
        class_tab <- cbind(id = rownames(class_tab), class_tab)
        write.table(class_tab, file.path(outdir, "Classification_matrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        fea_list <- lapply(model, ExtractVar)
        fea_df <- do.call(rbind, lapply(names(fea_list), function(method) {
            data.frame(features = fea_list[[method]], algorithm = method, stringsAsFactors = FALSE)
        }))
        write.table(fea_df, file.path(outdir, "113_model_slected_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        auc_list <- lapply(methods_valid, function(method) {
            RunEval(fit = model[[method]], Test_set = ml$test_x, Test_label = ml$test_y, Train_set = ml$train_x, Train_label = ml$train_y, 
                Train_name = "Train", cohortVar = "Cohort", classVar = class_var)
        })
        names(auc_list) <- methods_valid
        auc_mat <- do.call(rbind, auc_list)
        auc_tab <- cbind(Method = rownames(auc_mat), auc_mat)
        write.table(auc_tab, file.path(outdir, "113_model.AUC.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        avg_auc <- sort(apply(auc_mat, 1, mean, na.rm = TRUE), decreasing = TRUE)
        best_method <- names(avg_auc)[1]
        sig_genes <- unique(intersect(fea_list[[best_method]], colnames(ml$train_x)))
        if (length(sig_genes) == 0) 
            sig_genes <- unique(intersect(ml$genes, colnames(ml$train_x)))
        write_compatibility_outputs(outdir, sig_genes, best_method, auc_mat)
        plot_auc_heatmap(auc_mat, outdir)
        plot_model_roc(risk_tab, combined_labels, best_method, outdir)
        plot_gene_roc(ml$expr_list, ml$labels, sig_genes, outdir)
        invisible(list(best_method = best_method, key_genes = sig_genes, auc = auc_mat))
    }
    library(clusterProfiler)
    library(ggrepel)
    library(cowplot)
    library(ggraph)
    library(limma)
    gather_graph_node <- function(df, index = NULL, value = tail(colnames(df), 1), root = NULL) {
        require(dplyr)
        if (!length(index) < 2) {
            list <- lapply(seq_along(index), function(i) {
                dots <- index[1:i]
                df %>% group_by(.dots = dots) %>% summarise(node.size = sum(.data[[value]]), node.level = index[[i]], node.count = n()) %>% 
                  mutate(node.short_name = as.character(.data[[dots[[length(dots)]]]]), node.branch = as.character(.data[[dots[[1]]]])) %>% 
                  tidyr::unite(node.name, dots, sep = "/")
            })
            data <- do.call("rbind", list) %>% as_tibble()
            data$node.level <- factor(data$node.level, levels = index)
            if (is.null(root)) {
                return(data)
            }
            else {
                root_data <- data.frame(node.name = root, node.size = sum(df[[value]]), node.level = root, node.count = 1, 
                  node.short_name = root, node.branch = root, stringsAsFactors = F)
                data <- rbind(root_data, data)
                data$node.level <- factor(data$node.level, levels = c(root, index))
                return(data)
            }
        }
    }
    gather_graph_edge <- function(df, index = NULL, root = NULL) {
        require(dplyr)
        if (!length(index) < 2) 
            if (length(index) == 2) {
                data <- df %>% mutate(from = .data[[index[[1]]]]) %>% tidyr::unite(to, index, sep = "/") %>% dplyr::select(from, 
                  to) %>% mutate_at(c("from", "to"), as.character)
            }
            else {
                list <- lapply(seq(2, length(index)), function(i) {
                  dots <- index[1:i]
                  df %>% tidyr::unite(from, dots[-length(dots)], sep = "/", remove = F) %>% tidyr::unite(to, dots, sep = "/") %>% 
                    dplyr::select(from, to) %>% mutate_at(c("from", "to"), as.character)
                })
                data <- do.call("rbind", list)
            }
        data <- as_tibble(data)
        if (is.null(root)) {
            return(data)
        }
        else {
            root_data <- df %>% group_by(.dots = index[[1]]) %>% summarise(count = n()) %>% mutate(from = root, to = as.character(.data[[index[[1]]]])) %>% 
                dplyr::select(from, to)
            rbind(root_data, data)
        }
    }
    GSEA_module <- function(mat, outpath = Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo"), species = NA, KeyGene = NA, 
        isKEGG = TRUE) {
        dir <- list.files(outpath, pattern = "GSEA", full.names = TRUE, recursive = FALSE)
        if (length(dir) == 0) {
            dir <- paste0(outpath, "/GSEA/")
            dir.create(dir)
        }
        else {
            dir <- dir[1]
        }
        if (!is.na(species)) {
            mouse <- ifelse(species == "mouse", T, F)
        }
        else {
            mouse <- F
        }
        for (GeneIndex in 1:length(KeyGene)) {
            gene <- KeyGene[GeneIndex]
            low <- mat[gene, ] <= median(mat[gene, ])
            high <- mat[gene, ] > median(mat[gene, ])
            lowRT <- mat[, low]
            highRT <- mat[, high]
            conNum <- ncol(lowRT)
            treatNum <- ncol(highRT)
            mat <- cbind(lowRT, highRT)
            Type <- c(rep("con", conNum), rep("treat", treatNum))
            design <- model.matrix(~0 + factor(Type))
            colnames(design) <- c("con", "treat")
            fit <- lmFit(mat, design)
            cont.matrix <- makeContrasts(treat - con, levels = design)
            fit2 <- contrasts.fit(fit, cont.matrix)
            fit2 <- eBayes(fit2)
            allDiff <- topTable(fit2, adjust = "fdr", number = 2e+05)
            gsym.fc <- data.frame(SYMBOL = rownames(allDiff), logFC = allDiff$logFC)
            if (!dir.exists(paste0(dir, "/", GeneIndex, ".", gene))) {
                dir.create(paste0(dir, "/", GeneIndex, ".", gene))
            }
            selectedGeneID <- gene
            mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9", "#7A142C", "#5D90BA", 
                "#431A3D", "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D", "#7CC767")
            if (isKEGG) {
                if (mouse) {
                  gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
                  gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                  gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                  id.fc <- gsym.fc.id.sorted$logFC
                  names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                  kk <- gseKEGG(id.fc, organism = "mmu")
                  kk <- append_kegg_category(kk)
                  kk.gsym <- setReadable(kk, "org.Mm.eg.db", "ENTREZID")
                  sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
                }
                else {
                  gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                  gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                  gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                  id.fc <- gsym.fc.id.sorted$logFC
                  names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                  kegg.local <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
                  term2gene <- unique(data.frame(term = kegg.local$gs_name, gene = as.character(kegg.local$entrez_gene), 
                    stringsAsFactors = FALSE))
                  kk <- GSEA(id.fc, TERM2GENE = term2gene, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, verbose = FALSE)
                  kk@result$category <- "KEGG"
                  kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
                  sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
                }
                write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/GSEA_output.csv"), quote = F, row.names = F)
            }
            else {
                if (mouse) {
                  gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
                  gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                  gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                  id.fc <- gsym.fc.id.sorted$logFC
                  names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                  kk <- gseGO(id.fc, OrgDb = "org.Mm.eg.db", pAdjustMethod = "none")
                  kk.gsym <- setReadable(kk, "org.Mm.eg.db", "ENTREZID")
                  sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
                }
                else {
                  gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
                  gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
                  gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
                  id.fc <- gsym.fc.id.sorted$logFC
                  names(id.fc) <- gsym.fc.id.sorted$ENTREZID
                  kk <- gseGO(id.fc, OrgDb = "org.Hs.eg.db", pAdjustMethod = "none")
                  kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
                  sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
                }
                write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/GSEA_output.csv"), quote = F, row.names = F)
                if (nrow(sortkk) < 3) {
                  Sys.sleep(3)
                  next
                }
            }
            if (isKEGG) {
                pathway_keep <- !is.na(sortkk$category) & sortkk$category != "Human Diseases" & !is.na(sortkk$ID) & nzchar(sortkk$ID)
                geneSetID <- head(sortkk$ID[pathway_keep], 3)
                description.grep <- sortkk$Description[match(geneSetID, sortkk$ID)]
            }
            else {
                geneSetID <- head(sortkk$ID[!is.na(sortkk$ID) & nzchar(sortkk$ID)], 3)
                description.grep <- sortkk$Description[match(geneSetID, sortkk$ID)]
            }
            geneSetID <- geneSetID[!is.na(geneSetID) & nzchar(geneSetID)]
            description.grep <- description.grep[!is.na(description.grep) & nzchar(description.grep)]
            if (length(geneSetID) == 0 || length(description.grep) == 0) {
                write.table(data.frame(gene = gene, reason = "no_plot_eligible_pathway", stringsAsFactors = FALSE), paste0(dir, 
                  "/", GeneIndex, ".", gene, "/", "GSEA_plot_skip.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
                next
            }
            x <- kk
            geneList <- position <- NULL
            gsinfo_list <- lapply(geneSetID, function(pathway_id) {
                tryCatch(enrichplot:::gsInfo(pathway_id, object = x), error = function(e) {
                  write.table(data.frame(gene = gene, pathway_id = pathway_id, reason = conditionMessage(e), stringsAsFactors = FALSE), 
                    paste0(dir, "/", GeneIndex, ".", gene, "/", "GSEA_plot_skip.tsv"), quote = FALSE, row.names = FALSE, 
                    sep = "\t", append = file.exists(paste0(dir, "/", GeneIndex, ".", gene, "/", "GSEA_plot_skip.tsv")), 
                    col.names = !file.exists(paste0(dir, "/", GeneIndex, ".", gene, "/", "GSEA_plot_skip.tsv")))
                  NULL
                })
            })
            gsinfo_list <- gsinfo_list[!vapply(gsinfo_list, is.null, logical(1))]
            if (length(gsinfo_list) == 0) {
                next
            }
            gsdata <- do.call(rbind, gsinfo_list)
            description.grep <- unique(gsdata$Description)
            gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL, length(gsinfo_list))
            p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + geom_line(aes_(y = ~runningScore, color = ~Description), 
                size = 1) + scale_color_manual(values = mycol, labels = function(x) stringr::str_wrap(gsub("_", " ", sub("^KEGG_", 
                "", x)), width = 32)) + geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + ylab("Enrichment\n Score") + 
                theme_bw() + theme(panel.grid = element_blank()) + theme(legend.position = "top", legend.title = element_blank(), 
                legend.text = element_text(size = 8), legend.background = element_rect(fill = "transparent")) + guides(color = guide_legend(nrow = 2, 
                byrow = TRUE)) + theme(axis.text.y = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(), axis.line.x = element_blank(), plot.margin = margin(t = 0.2, r = 0.2, b = 0, 
                  l = 0.2, unit = "cm"))
            rel_heights <- c(1.5, 0.5, 1.5)
            p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + 
                xlab(NULL) + ylab(NULL) + scale_color_manual(values = mycol) + theme_bw() + theme(panel.grid = element_blank()) + 
                theme(legend.position = "none", plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                  axis.text = element_blank(), axis.line.x = element_blank()) + scale_y_continuous(expand = c(0, 0))
            df2 <- p.res$data
            df2$y <- p.res$data$geneList[df2$x]
            df2$gsym <- p.res$data$gsym[df2$x]
            selectgenes <- data.frame(gsym = selectedGeneID)
            selectgenes <- merge(selectgenes, df2, by = "gsym")
            selectgenes <- selectgenes[selectgenes$position == 1, ]
            p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + geom_segment(data = df2, 
                aes_(x = ~x, xend = ~x, y = ~y, yend = 0), color = "grey") + geom_bar(position = "dodge", stat = "identity") + 
                scale_fill_manual(values = mycol, guide = FALSE) + scale_color_manual(values = mycol, guide = FALSE) + geom_hline(yintercept = 0, 
                lty = 2, lwd = 0.2) + ylab("Ranked list\n metric") + xlab("Rank in ordered dataset") + theme_bw() + theme(axis.text.y = element_text(size = 12, 
                face = "bold"), panel.grid = element_blank()) + geom_text_repel(data = selectgenes, show.legend = FALSE, 
                direction = "x", ylim = c(2, NA), angle = 90, size = 2.5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, 
                  "lines")) + theme(plot.margin = margin(t = -0.1, r = 0.4, b = 0.2, l = 0.2, unit = "cm"))
            plotlist <- list(p.res, p2, p.pos)
            plotlistNum <- length(plotlist)
            plotlist[[plotlistNum]] <- plotlist[[plotlistNum]] + theme(axis.line.x = element_line(), axis.ticks.x = element_line(), 
                axis.text.x = element_text(size = 12, face = "bold"))
            p <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
            ggsave(paste0(dir, "/", GeneIndex, ".", gene, "/", "/1.GSEA_multi_pathways.pdf"), width = 10, height = 7, plot = p)
            sortkk <- kk.gsym[kk.gsym@result$Description %in% description.grep, ]
            go <- data.frame(Category = "KEGG", ID = sortkk$ID, Term = sortkk$Description, Genes = gsub("/", ", ", sortkk$core_enrichment), 
                adj_pval = sortkk$p.adjust)
            genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)
            df <- GOplot::circle_dat(go, genelist)[, c(3, 5, 6)]
            very_easy <- data.frame(term = unique(description.grep), stringsAsFactors = FALSE)
            very_easy <- very_easy[!is.na(very_easy$term) & nzchar(very_easy$term), , drop = FALSE]
            if (nrow(very_easy) > 0) {
                write.table(very_easy, paste0(dir, "/", GeneIndex, ".", gene, "/", "very_easy_input.txt"), quote = F, row.names = F, 
                  sep = "\t")
            }
            nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root = "all")
            edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
            nodes <- nodes %>% mutate_at(c("node.level", "node.branch"), as.character)
            graph <- tidygraph::tbl_graph(nodes, edges)
            width <- 600/(nrow(nodes) + 100)
            gc1 <- ggraph(graph, layout = "dendrogram", circular = TRUE) + geom_edge_diagonal(aes(color = node1.node.branch, 
                filter = node1.node.level != "all"), alpha = 0.5, edge_width = 2.5) + scale_edge_color_manual(values = c("#61C3ED", 
                "red", "purple", "darkgreen")) + geom_node_point(aes(size = node.size, filter = node.level != "all"), color = "#61C3ED") + 
                scale_size(range = c(width, 10)) + theme(legend.position = "none") + geom_node_text(aes(x = 1.05 * x, y = 1.05 * 
                y, label = node.short_name, angle = -((-node_angle(x, y) + 90)%%180) + 90, filter = leaf), color = "black", 
                size = 3, hjust = "outward") + geom_node_text(aes(label = node.short_name, filter = !leaf & (node.level != 
                "all")), color = "black", fontface = "bold", size = 6, family = "sans") + theme(panel.background = element_rect(fill = NA)) + 
                coord_cartesian(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3))
            ggsave(paste0(dir, "/", GeneIndex, ".", gene, "/2.Ccgraph.pdf"), width = 14, height = 14, plot = gc1)
        }
    }
    CoreAlg <- function(X, y) {
        svn_itor <- 3
        res <- function(i) {
            if (i == 1) {
                nus <- 0.25
            }
            if (i == 2) {
                nus <- 0.5
            }
            if (i == 3) {
                nus <- 0.75
            }
            model <- svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = F)
            model
        }
        if (Sys.info()["sysname"] == "Windows") {
            out <- mclapply(1:svn_itor, res, mc.cores = 1)
        }
        else {
            out <- mclapply(1:svn_itor, res, mc.cores = svn_itor)
        }
        nusvm <- rep(0, svn_itor)
        corrv <- rep(0, svn_itor)
        t <- 1
        while (t <= svn_itor) {
            weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
            weights[which(weights < 0)] <- 0
            w <- weights/sum(weights)
            u <- sweep(X, MARGIN = 2, w, "*")
            k <- apply(u, 1, sum)
            nusvm[t] <- sqrt((mean((k - y)^2)))
            corrv[t] <- cor(k, y)
            t <- t + 1
        }
        rmses <- nusvm
        mn <- which.min(rmses)
        model <- out[[mn]]
        q <- t(model$coefs) %*% model$SV
        q[which(q < 0)] <- 0
        w <- (q/sum(q))
        mix_rmse <- rmses[mn]
        mix_r <- corrv[mn]
        newList <- list(w = w, mix_rmse = mix_rmse, mix_r = mix_r)
    }
    doPerm <- function(perm, X, Y) {
        itor <- 1
        Ylist <- as.list(data.matrix(Y))
        dist <- matrix()
        while (itor <= perm) {
            yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
            yr <- (yr - mean(yr))/sd(yr)
            result <- CoreAlg(X, yr)
            mix_r <- result$mix_r
            if (itor == 1) {
                dist <- mix_r
            }
            else {
                dist <- rbind(dist, mix_r)
            }
            itor <- itor + 1
        }
        newList <- list(dist = dist)
    }
    CIBERSORT <- function(dir, sig_matrix, mixture_file, perm = 0, QN = TRUE) {
        library(e1071)
        library(parallel)
        library(preprocessCore)
        X <- read.table(sig_matrix, header = T, sep = "\t", row.names = 1, check.names = F)
        Y <- mixture_file
        X <- data.matrix(X)
        Y <- data.matrix(Y)
        X <- X[order(rownames(X)), ]
        Y <- Y[order(rownames(Y)), ]
        P <- perm
        if (max(Y) < 50) {
            Y <- 2^Y
        }
        if (QN == TRUE) {
            tmpc <- colnames(Y)
            tmpr <- rownames(Y)
            Y <- normalize.quantiles(Y)
            colnames(Y) <- tmpc
            rownames(Y) <- tmpr
        }
        if (substr(Sys.Date(), 6, 7) > 5) {
            next
        }
        Xgns <- row.names(X)
        Ygns <- row.names(Y)
        YintX <- Ygns %in% Xgns
        Y <- Y[YintX, ]
        XintY <- Xgns %in% row.names(Y)
        X <- X[XintY, ]
        X <- (X - mean(X))/sd(as.vector(X))
        if (P > 0) {
            nulldist <- sort(doPerm(P, X, Y)$dist)
        }
        header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
        output <- matrix()
        itor <- 1
        mixtures <- dim(Y)[2]
        pval <- 9999
        while (itor <= mixtures) {
            y <- Y[, itor]
            y <- (y - mean(y))/sd(y)
            result <- CoreAlg(X, y)
            if (substr(Sys.Date(), 1, 4) > 3021) {
                next
            }
            w <- result$w
            mix_r <- result$mix_r
            mix_rmse <- result$mix_rmse
            if (P > 0) {
                pval <- 1 - (which.min(abs(nulldist - mix_r))/length(nulldist))
            }
            out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
            if (itor == 1) {
                output <- out
            }
            else {
                output <- rbind(output, out)
            }
            itor <- itor + 1
        }
        write.table(rbind(header, output), file = paste0(dir, "/CIBERSORT-Results.txt"), sep = "\t", row.names = F, col.names = F, 
            quote = F)
        obj <- rbind(header, output)
        obj <- obj[, -1]
        obj <- obj[-1, ]
        obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
        rownames(obj) <- colnames(Y)
        colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
        obj
    }
    save_pdf <- function(path, plot, width, height) {
        ggsave(path, plot = plot, width = width, height = height, device = cairo_pdf, limitsize = FALSE)
    }
    obj <- readRDS(file.path(input_dir, "seurat.annotated.rds"))
    md <- as.data.table(obj@meta.data)
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
    ribo <- grep("^RP[SL]", rownames(counts), value = TRUE)
    denom <- Matrix::colSums(counts)
    md[, `:=`(percent.ribo, if (length(ribo)) 
        100 * Matrix::colSums(counts[ribo, , drop = FALSE])/pmax(denom, 1)
    else 0)]
    md[, `:=`(percent.mt, percent_mt)]
    qc_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
    qsum <- rbindlist(lapply(qc_vars, function(v) {
        md[, .(q05 = quantile(get(v), 0.05, na.rm = TRUE), q25 = quantile(get(v), 0.25, na.rm = TRUE), med = median(get(v), 
            na.rm = TRUE), q75 = quantile(get(v), 0.75, na.rm = TRUE), q95 = quantile(get(v), 0.95, na.rm = TRUE)), by = dataset_id][, 
            `:=`(metric, v)]
    }))
    p_qc <- ggplot(qsum, aes(dataset_id, ymin = q05, ymax = q95, lower = q25, upper = q75, middle = med, fill = dataset_id)) + 
        geom_boxplot(stat = "identity", width = 0.7, linewidth = 0.25) + facet_wrap(~metric, scales = "free_y", nrow = 1) + 
        theme_classic() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), legend.position = "none", strip.background = element_blank()) + 
        labs(x = NULL, y = NULL)
    save_pdf(file.path(out, "1.QCFilter", "1.VlnPlot.pdf"), p_qc, 20, 4)
    p_qc1 <- ggplot(md, aes(nCount_RNA, nFeature_RNA)) + geom_hex(bins = 120) + scale_fill_viridis_c(trans = "log10") + theme_classic()
    p_qc2 <- ggplot(md, aes(nCount_RNA, percent.mt)) + geom_hex(bins = 120) + scale_fill_viridis_c(trans = "log10") + theme_classic()
    p_qc3 <- ggplot(md, aes(nCount_RNA, percent.ribo)) + geom_hex(bins = 120) + scale_fill_viridis_c(trans = "log10") + theme_classic()
    save_pdf(file.path(out, "1.QCFilter", "2.FeatureScatter.pdf"), p_qc1 + p_qc2 + p_qc3, 21, 7)
    sdv <- Stdev(obj, reduction = "pca")
    p_elbow <- ggplot(data.frame(PC = seq_along(sdv), SD = sdv), aes(PC, SD)) + geom_point(size = 2) + geom_line() + theme_classic() + 
        labs(x = "PC", y = "Standard deviation")
    save_pdf(file.path(out, "2.Cluster", "2.ElbowPlot.pdf"), p_elbow, 7, 7)
    p_pca <- DimPlot(obj, reduction = "pca", group.by = "dataset_id", raster = TRUE, pt.size = 1)
    save_pdf(file.path(out, "2.Cluster", "3.PcaPlot.pdf"), p_pca, 10, 8)
    p_harmony <- DimPlot(obj, reduction = "harmony", group.by = "dataset_id", raster = TRUE, pt.size = 1)
    save_pdf(file.path(out, "2.Cluster", "4.HarmonyPlot.pdf"), p_harmony, 10, 8)
    cluster_values <- as.character(obj$seurat_clusters)
    cluster_levels <- as.character(sort(unique(as.integer(cluster_values))))
    obj$preannotation_cluster <- factor(cluster_values, levels = cluster_levels)
    cluster_colors <- setNames((scales::hue_pal())(length(cluster_levels)), cluster_levels)
    set.seed(20260720)
    p_cluster <- DimPlot(obj, reduction = "umap", group.by = "preannotation_cluster", cols = unname(cluster_colors), label = TRUE, 
        repel = TRUE, raster = TRUE, pt.size = 1) + NoLegend()
    save_pdf(file.path(out, "2.Cluster", "5.Cluster.pdf"), p_cluster, 10, 8)
    celltype_levels <- levels(factor(obj$cellType_1))
    obj$cellType_1 <- factor(obj$cellType_1, levels = celltype_levels)
    p_annotate <- DimPlot(obj, reduction = "umap", group.by = "cellType_1", label = TRUE, repel = TRUE, raster = TRUE, pt.size = 1) + 
        NoLegend()
    save_pdf(file.path(out, "3.CellAnnotate", "1.Umap.Cluster.pdf"), p_annotate, 10, 8)
    pdc_all <- readRDS(file.path(input_dir, "pDC_fullcell_no_sampling.rds"))
    pdc_emb <- Embeddings(pdc_all, "umap")
    pdc_meta <- data.frame(UMAP_1 = pdc_emb[, 1], UMAP_2 = pdc_emb[, 2], phenotype = as.character(pdc_all$cellType_2), nFeature_RNA = pdc_all$nFeature_RNA)
    cyto_all <- CytoTRACE::CytoTRACE(as.matrix(GetAssayData(pdc_all, assay = "RNA", slot = "counts")), enableFast = FALSE)$CytoTRACE
    pdc_meta$CytoTRACE <- as.numeric(cyto_all)
    p_cyto_umap <- ggplot(pdc_meta, aes(UMAP_1, UMAP_2, color = CytoTRACE)) + geom_point(size = 0.35, alpha = 0.75) + scale_color_viridis_c() + 
        theme_classic()
    save_pdf(file.path(out, "6.CytoTRACE", "1.CytoTRACE_plot.pdf"), p_cyto_umap, 7, 6)
    p_cyto_box <- ggplot(pdc_meta, aes(phenotype, CytoTRACE, fill = phenotype)) + geom_violin(scale = "width", alpha = 0.65) + 
        geom_boxplot(width = 0.14, outlier.shape = NA) + theme_classic() + labs(x = NULL)
    save_pdf(file.path(out, "6.CytoTRACE", "2.CytoTRACE_by_Phenotype_Boxplot.pdf"), p_cyto_box, 7, 5)
    rm(counts, md, obj, pdc_all, pdc_emb, pdc_meta, cyto_all)
    invisible(gc())
    local({
        obj <- readRDS(file.path(input_dir, "seurat.annotated.rds"))
        stopifnot(all(levels(obj$group) == c("Control", "Disease")))
        if (!file.exists(file.path(out, "2.Cluster", "1.VariableFeaturePlot.pdf"))) {
            obj <- NormalizeData(obj, verbose = FALSE)
            obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
            top10 <- head(VariableFeatures(obj), 10)
            p1 <- VariableFeaturePlot(obj)
            p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
            ggsave(file.path(out, "2.Cluster", "1.VariableFeaturePlot.pdf"), plot = p1 + p2, height = 7, width = 14)
        }
        marker_file <- file.path(stage, "legacy_markers_cellType_1.csv")
        if (file.exists(marker_file)) {
            markers <- read.csv(marker_file, check.names = FALSE)
        }
        else {
            markers <- make_marker_table_sparse(obj, celltype = "cellType_1")
            write.csv(markers, marker_file, row.names = FALSE)
        }
        marker_features <- intersect(unique(head(markers$gene, 80)), rownames(obj))
        if (!file.exists(file.path(out, "3.CellAnnotate", "3.Dotplot.pdf"))) {
            pdot <- DotPlot(obj, features = marker_features, cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 
                10)) + theme(panel.border = element_rect(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, 
                vjust = 0.5), legend.position = "top", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, 
                "cm"))
            ggsave(file.path(out, "3.CellAnnotate", "3.Dotplot.pdf"), plot = pdot, width = 3 + length(marker_features) * 
                0.3, height = 3 + length(unique(obj$cellType_1)) * 0.5)
        }
        tmp_sc <- file.path(stage, "legacy_singlecell_module_out")
        dir.create(tmp_sc, recursive = TRUE, showWarnings = FALSE)
        existing_plotdata <- file.exists(file.path(tmp_sc, "ClusterGvis", "Heatmap.plotdata.rds"))
        heatmap_ok <- if (existing_plotdata) 
            FALSE
        else tryCatch({
            ClusterGVis_module(object = obj, OutPath = tmp_sc, scRNA.markers = markers, celltype = "cellType_1", species = NA, 
                OnlyPlot = FALSE)
            TRUE
        }, error = function(e) {
            FALSE
        })
        if (!heatmap_ok) {
            suppressPackageStartupMessages(library(ClusterGVis))
            cg <- file.path(tmp_sc, "ClusterGvis")
            st.data <- readRDS(file.path(cg, "Heatmap.plotdata.rds"))
            topm <- read.csv(file.path(cg, "TopMarkers.csv"), check.names = FALSE)
            cluster_order <- sort(unique(as.integer(as.character(st.data$wide.res$cluster))))
            pdf(file.path(cg, "1.Heatmap.pdf"), height = 15, width = 17, onefile = FALSE)
            visCluster(object = st.data, plot.type = "both", column_title_rot = 45, markGenes = unique(topm$gene), markGenes.side = "left", 
                line.side = "left", cluster.order = cluster_order, add.bar = TRUE)
            dev.off()
        }
        file.copy(file.path(tmp_sc, "ClusterGvis", "1.Heatmap.pdf"), file.path(out, "4.ClusterGvis", "1.Heatmap.pdf"), overwrite = TRUE)
        sig <- fread("data/signature.txt", check.names = FALSE)
        sig_list <- list(Peroxisome = sig[[1]])
        sig_pdfs <- file.path(out, "5.Signature", c("1.Dotplot of Peroxisome.pdf", "2.Boxplot of Peroxisome.pdf", "3.Volcano diagram of Peroxisome.pdf"))
        sig_done <- all(c("AUCell", "UCell", "singscore", "ssGSEA", "AddModuleScore", "Scoring", "cellType_2") %in% colnames(obj@meta.data)) && 
            all(file.exists(sig_pdfs))
        if (!sig_done) {
            obj <- signature_module(object = obj, outpath = out, signature = sig_list, group_by = "group", celltype = "cellType_1", 
                OnlyPlot = FALSE)
            old_ct <- as.character(obj$cellType_1)
            score_median <- median(obj$Scoring[old_ct == "pDC"], na.rm = TRUE)
            obj$cellType_2 <- ifelse(old_ct == "pDC" & obj$Scoring > score_median, "Hsco_pDC", ifelse(old_ct == "pDC" & obj$Scoring <= 
                score_median, "Lsco_pDC", old_ct))
            obj$cellType_2 <- factor(obj$cellType_2)
            signature_volcano(object = obj, outpath = out, ShowCell = c("Hsco_pDC", "Lsco_pDC"), group_by = "cellType_2", 
                celltype = "cellType_2", idents = c("Hsco_pDC", "Lsco_pDC"), OnlyPlot = FALSE)
        }
        pdc <- subset(obj, subset = cellType_2 %in% c("Hsco_pDC", "Lsco_pDC"))
        monocle_files <- c("1.Trajectory.Pseudotime.pdf", "2.Trajectory.State.pdf", "3.Trajectory.Cluster.pdf")
        if (!all(file.exists(file.path(tmp_sc, "Monocle2", monocle_files)))) {
            monocle2_module(object = pdc, OutPath = tmp_sc, group_by = "cellType_2", ShowGene = NA, OnlyPlot = FALSE)
        }
        for (f in monocle_files) {
            from <- file.path(tmp_sc, "Monocle2", f)
            file.copy(from, file.path(out, "7.Monocle2", f), overwrite = TRUE)
        }
        cellchat_rds <- file.path(tmp_sc, "CellChat", "cellchat.rds")
        cellchat_ok <- if (file.exists(cellchat_rds)) 
            FALSE
        else tryCatch({
            cellchat_module(object = obj, OutPath = tmp_sc, group_by = "all", celltype = "cellType_2", ShowCell = "Hsco_pDC", 
                species = NA, OnlyPlot = FALSE, AllPlot = FALSE)
            TRUE
        }, error = function(e) {
            FALSE
        })
        if (!cellchat_ok) {
            suppressPackageStartupMessages(library(CellChat))
            suppressPackageStartupMessages(library(patchwork))
            ccdir <- file.path(tmp_sc, "CellChat")
            cellchat <- readRDS(file.path(ccdir, "cellchat.rds"))
            groupSize <- as.numeric(table(cellchat@idents))
            count_mat <- cellchat@net$count
            weight_mat <- cellchat@net$weight
            diag(count_mat) <- 0
            diag(weight_mat) <- 0
            pdf(file.path(ccdir, "1.Net_number_strength.pdf"), width = 10, height = 5)
            par(mfrow = c(1, 2))
            netVisual_circle(count_mat, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
            netVisual_circle(weight_mat, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
            dev.off()
            df.net <- read.table(file.path(ccdir, "net_lr.txt"), sep = "\t", check.names = FALSE, header = TRUE)
            bardata <- as.data.frame(table(c(df.net$source, df.net$target)))
            colnames(bardata) <- c("Cell_Type", "all_sum")
            bardata <- bardata[order(bardata$all_sum, decreasing = TRUE), ]
            bardata$Cell_Type <- factor(bardata$Cell_Type, levels = bardata$Cell_Type)
            sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
            sample_color2 <- ifelse(bardata$all_sum > summary(bardata$all_sum)[2], ifelse(bardata$all_sum > summary(bardata$all_sum)[3], 
                ifelse(bardata$all_sum > summary(bardata$all_sum)[5], sample_color[1], sample_color[2]), sample_color[3]), 
                sample_color[4])
            sample_color1 <- scales::alpha(sample_color2, 0.9)
            pB2 <- ggplot(bardata, aes(Cell_Type, all_sum, fill = Cell_Type, colour = Cell_Type)) + geom_bar(stat = "identity", 
                width = 0.8) + scale_fill_manual(values = sample_color1) + scale_colour_manual(values = sample_color2) + 
                theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_text(size = 12, colour = "black"), 
                axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black")) + labs(x = "Cell Type", 
                y = "Counts")
            ggsave(file.path(ccdir, "2.Interaction Count.pdf"), plot = pB2, width = 10, height = 8)
            ident1 <- which(levels(cellchat@idents) %in% "Hsco_pDC")
            ident2 <- which(!levels(cellchat@idents) %in% "Hsco_pDC")
            gg1 <- netVisual_bubble(cellchat, sources.use = ident1, targets.use = ident2, angle.x = 45, remove.isolate = FALSE, 
                font.size = 12, font.size.title = 15, return.data = TRUE, title.name = "Signaling from Hsco_pDC")
            gg2 <- netVisual_bubble(cellchat, sources.use = ident2, targets.use = ident1, angle.x = 45, remove.isolate = FALSE, 
                font.size = 12, font.size.title = 15, return.data = TRUE, title.name = "Signaling to Hsco_pDC")
            ggsave(file.path(ccdir, "3.Bubble.pdf"), plot = gg1$gg.obj + gg2$gg.obj, width = 5 + length(levels(cellchat@idents)), 
                height = 8 + max(nrow(gg1$communication), nrow(gg2$communication)) * 0.05)
        }
        for (f in c("1.Net_number_strength.pdf", "2.Interaction Count.pdf", "3.Bubble.pdf")) {
            from <- file.path(tmp_sc, "CellChat", f)
            file.copy(from, file.path(out, "8.CellChat", f), overwrite = TRUE)
        }
    })
    local({
        ml_out <- file.path(out, "9.Machine learning")
        res_ml <- run_113ml_pipeline(candidate_genes = keygenes, outdir = ml_out, data_dir = "data")
        fwrite(data.table(KeyGene = keygenes), file.path(ml_out, "sigGenes.txt"), sep = "\t")
        ml_data <- load_expr_and_labels("data")
        plot_gene_roc(ml_data$expr, ml_data$labels, keygenes, ml_out)
        pb_pdc_table <- fread(file.path(input_dir, "patient_pseudobulk_pDC.tsv"), check.names = FALSE)
        pb_pdc <- as.matrix(pb_pdc_table[, -1, with = FALSE])
        rownames(pb_pdc) <- pb_pdc_table[[1]]
        patient_meta <- fread(file.path(input_dir, "patient_metadata.tsv"))
        patient_meta <- patient_meta[group %in% c("Healthy", "Sepsis")]
        main_samples <- intersect(colnames(pb_pdc), patient_meta$patient_uid)
        main_group <- factor(ifelse(patient_meta[match(main_samples, patient_uid), group] == "Healthy", "Control", "Disease"), 
            levels = c("Control", "Disease"))
        x_cv <- t(pb_pdc[keygenes, main_samples, drop = FALSE])
        y_cv <- as.integer(main_group == "Disease")
        df_cv <- data.frame(Group = main_group, x_cv, check.names = FALSE)
        set.seed(10007)
        fold <- integer(nrow(df_cv))
        for (lev in levels(main_group)) {
            idx <- which(main_group == lev)
            fold[idx] <- sample(rep(1:5, length.out = length(idx)))
        }
        pred_cv <- matrix(NA_real_, nrow(df_cv), 3, dimnames = list(NULL, c("Logistic", "LASSO", "Random forest")))
        for (k in 1:5) {
            train <- fold != k
            test <- !train
            fit_glm_cv <- glm(Group ~ ., df_cv[train, ], family = binomial())
            pred_cv[test, "Logistic"] <- predict(fit_glm_cv, df_cv[test, ], type = "response")
            fit_lasso_cv <- cv.glmnet(as.matrix(x_cv[train, , drop = FALSE]), y_cv[train], family = "binomial", alpha = 1, 
                nfolds = min(5, sum(train)))
            pred_cv[test, "LASSO"] <- as.numeric(predict(fit_lasso_cv, as.matrix(x_cv[test, , drop = FALSE]), s = "lambda.min", 
                type = "response"))
            fit_rf_cv <- randomForest(x = x_cv[train, , drop = FALSE], y = main_group[train], ntree = 800)
            pred_cv[test, "Random forest"] <- predict(fit_rf_cv, x_cv[test, , drop = FALSE], type = "prob")[, "Disease"]
        }
        roc_cv <- lapply(seq_len(ncol(pred_cv)), function(i) {
            pROC::roc(main_group, pred_cv[, i], levels = c("Control", "Disease"), direction = "<", quiet = TRUE)
        })
        pdf(file.path(ml_out, "4.Test_5foldCV_ROC.genes.pdf"), width = 5, height = 4.75)
        plot(roc_cv[[1]], col = "#E64B35", lwd = 2.5, main = "Out-of-fold ROC")
        plot(roc_cv[[2]], add = TRUE, col = "#4DBBD5", lwd = 2.5)
        plot(roc_cv[[3]], add = TRUE, col = "#00A087", lwd = 2.5)
        abline(0, 1, lty = 2, col = "grey60")
        legend("bottomright", legend = sprintf("%s AUC %.3f", colnames(pred_cv), vapply(roc_cv, function(r) as.numeric(pROC::auc(r)), 
            numeric(1))), col = c("#E64B35", "#4DBBD5", "#00A087"), lwd = 2.5, bty = "n", cex = 0.8)
        dev.off()
        unlink(file.path(ml_out, c("1.gene_intersection_upset.pdf", "3.RandomForest.pdf")), force = TRUE)
        file.rename(file.path(ml_out, "ML_AUC_heatmap.pdf"), file.path(ml_out, "1.ML_AUC_heatmap.pdf"))
        file.rename(file.path(ml_out, "ROC.Train.pdf"), file.path(ml_out, "2.ROC.Train.pdf"))
        file.rename(file.path(ml_out, "Train_ROC.genes.pdf"), file.path(ml_out, "3.Train_ROC.genes.pdf"))
        file.rename(file.path(ml_out, "Test_GSE290679_ROC.genes.pdf"), file.path(ml_out, "5.Test_GSE290679_ROC.genes.pdf"))
        file.copy("data/expr_all_main.txt", "data/expr.txt", overwrite = TRUE)
        file.copy("data/metadata_all_main.txt", "data/metadata.txt", overwrite = TRUE)
        local({
            base_out <- Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo")
            mat.path <- "./data/expr.txt"
            keygene.path <- list.files(path = base_out, recursive = TRUE, pattern = "sigGenes\\.txt$", full.names = TRUE)
            mat <- data.table::fread(mat.path, data.table = FALSE, header = TRUE)
            mat <- limma::avereps(mat[, -1], ID = mat[, 1])
            if (range(mat)[1] >= 0 & range(mat)[2] - range(mat)[1] > 100) {
                mat <- log2(mat + 1)
            }
            KeyGene <- read.table(keygene.path[1], header = TRUE)
            GSEA_module(mat = mat, outpath = base_out, species = NA, KeyGene = KeyGene$KeyGene, isKEGG = TRUE)
        })
        local({
            library(GSVA)
            library(limma)
            library(ggplot2)
            base_out <- Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo")
            keygene.path <- list.files(path = base_out, recursive = TRUE, pattern = "sigGenes\\.txt$", full.names = TRUE)
            mat.path <- "data/expr.txt"
            outpath <- base_out
            keygene <- read.table(keygene.path[1], header = TRUE)
            KeyGene <- keygene$KeyGene
            mat <- data.table::fread(mat.path, data.table = FALSE, header = TRUE)
            mat <- avereps(mat[, -1], ID = mat[, 1])
            if (range(mat)[1] >= 0 & range(mat)[2] - range(mat)[1] > 100) {
                mat <- log2(mat + 1)
            }
            load("./data/hallmark.gs.rdata")
            gsva_es <- gsva(as.matrix(mat), gs)
            cutoff <- 1
            summary_rows <- list()
            x <- list.files(outpath, pattern = "GSVA", full.names = TRUE)
            if (length(x) == 0) {
                outpath <- paste0(outpath, "/GSVA/")
                dir.create(outpath)
            }
            else {
                outpath <- x[1]
            }
            for (i in KeyGene) {
                subexpr <- as.numeric(mat[i, ])
                names(subexpr) <- colnames(mat)
                lsam <- names(subexpr[subexpr < median(subexpr)])
                hsam <- names(subexpr[subexpr >= median(subexpr)])
                group_list <- data.frame(sample = c(lsam, hsam), group = c(rep("Lexp", length(lsam)), rep("Hexp", length(hsam))))
                design <- model.matrix(~0 + factor(group_list$group))
                colnames(design) <- levels(factor(group_list$group))
                rownames(design) <- colnames(gsva_es)
                contrast.matrix <- makeContrasts(Hexp - Lexp, levels = design)
                fit <- lmFit(gsva_es[, group_list$sample], design)
                fit2 <- contrasts.fit(fit, contrast.matrix)
                fit2 <- eBayes(fit2)
                x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
                pathway <- stringr::str_replace(row.names(x), "HALLMARK_", "")
                df <- data.frame(ID = pathway, score = x$t, adj.P.Val = x$adj.P.Val, stringsAsFactors = FALSE)
                summary_df <- df[order(df$adj.P.Val, -abs(df$score)), ]
                top_pathways <- head(summary_df$ID, 3)
                top_pathways <- c(top_pathways, rep(NA_character_, max(0, 3 - length(top_pathways))))
                summary_rows[[length(summary_rows) + 1]] <- data.frame(KeyGene = i, FigIndex = which(KeyGene == i), TopPathway1 = top_pathways[1], 
                  TopPathway2 = top_pathways[2], TopPathway3 = top_pathways[3], stringsAsFactors = FALSE)
                df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf), labels = c(1, 2, 3))
                sortdf <- df[order(df$score), ]
                sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
                p <- ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = c("palegreen3", 
                  "snow3", "dodgerblue4"), guide = FALSE) + geom_hline(yintercept = c(-cutoff, cutoff), color = "white", 
                  linetype = 2, size = 0.3) + geom_text(data = subset(df, score < 0), aes(x = ID, y = 0, label = paste0(" ", 
                  ID), color = group), size = 3, hjust = 0) + geom_text(data = subset(df, score > 0), aes(x = ID, y = -0.05, 
                  label = ID, color = group), size = 3, hjust = 1) + scale_colour_manual(values = c("black", "snow3", "black"), 
                  guide = FALSE) + xlab("") + ylab(paste0("t value of GSVA score\n Hexp vs Lexp group of ", i)) + theme_bw() + 
                  theme(panel.grid = element_blank(), panel.border = element_rect(size = 0.6), axis.line.y = element_blank(), 
                    axis.ticks.y = element_blank(), axis.text.y = element_blank())
                ggsave(paste0(outpath, "/", which(KeyGene == i), ".GSVA plot of ", i, ".pdf"), plot = p, width = 8, height = 7)
            }
            if (length(summary_rows) > 0) {
                gsva_summary <- do.call(rbind, summary_rows)
                write.table(gsva_summary, file = paste0(outpath, "/GSVA_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
            }
        })
        local({
            library(dplyr)
            library(ggplot2)
            library(corrplot)
            library(tidyHeatmap)
            library(tibble)
            library(grid)
            base_out <- Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo")
            outpath <- base_out
            mat.path <- "./data/expr.txt"
            keygene.path <- list.files(path = base_out, recursive = TRUE, pattern = "sigGenes\\.txt$", full.names = TRUE)
            mat <- data.table::fread(mat.path, data.table = FALSE, header = TRUE)
            mat <- limma::avereps(mat[, -1], ID = mat[, 1])
            gene <- read.table(keygene.path[1], header = TRUE)
            x <- list.files(outpath, pattern = "Immune infiltration", full.names = TRUE)
            if (length(x) == 0) {
                outpath <- paste0(outpath, "/Immune infiltration/")
                dir.create(outpath)
            }
            else {
                outpath <- x[1]
            }
            if (file.exists(paste0(outpath, "/CIBERSORT-Results.txt"))) {
                results <- read.table(paste0(outpath, "/CIBERSORT-Results.txt"), header = TRUE, sep = "\t", check.names = FALSE, 
                  row.names = 1)
            }
            else {
                results <- CIBERSORT(dir = outpath, "./data/cibersort.ref.txt", mat, perm = 100, QN = TRUE)
            }
            results <- data.frame(results, check.names = FALSE)
            data <- dplyr::select(results, -c("P-value", "Correlation", "RMSE"))
            data <- t(data)
            GseGroup <- read.table("./data/metadata.txt", header = FALSE, sep = "\t", check.names = FALSE)
            colnames(GseGroup) <- c("Acc", "Cluster", "Tissue")
            ciber <- data %>% t() %>% data.frame(check.names = FALSE)
            ciber$Acc <- rownames(ciber)
            if (sort(unique(GseGroup$Tissue))[1] == "Control") {
                datGroup <- GseGroup[order(GseGroup$Tissue), ]
            }
            else {
                datGroup <- GseGroup[order(GseGroup$Tissue, decreasing = TRUE), ]
            }
            datGroup$Acc <- factor(datGroup$Acc, levels = datGroup$Acc)
            data <- dplyr::inner_join(datGroup, ciber, by = "Acc")
            data_p <- reshape2::melt(data, id.vars = colnames(datGroup))
            data_p$Acc <- factor(data_p$Acc, levels = datGroup$Acc)
            p <- as.tibble(data_p) %>% group_by(Tissue) %>% tidyHeatmap::heatmap(.row = variable, .column = Acc, .value = value, 
                scale = "row", cluster_rows = T, palette_value = circlize::colorRamp2(seq(-2, 2, length.out = 11), RColorBrewer::brewer.pal(11, 
                  "RdBu")), palette_grouping = list(c("#66C2A5", "#FC8D62")), show_column_names = F, row_names_gp = gpar(fontsize = 10), 
                column_names_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 0), row_title_gp = gpar(fontsize = 0))
            pdf(paste0(outpath, "/1.Immune infiltration.pdf"), width = 10, height = 6)
            print(p)
            dev.off()
            data <- dplyr::select(results, -c("P-value", "Correlation", "RMSE"))
            data <- data[, colSums(data) != 0]
            pdf(paste0(outpath, "/2.CorHeatmap.pdf"), height = 13, width = 13)
            corrplot(corr = cor(data), method = "color", order = "hclust", tl.col = "black", addCoef.col = "black", number.cex = 0.8, 
                col = colorRampPalette(c("blue", "white", "red"))(50))
            dev.off()
            data <- dplyr::select(results, -c("P-value", "Correlation", "RMSE"))
            Type <- read.table("./data/metadata.txt", sep = "\t", check.names = FALSE, row.names = 1, header = FALSE)
            Type <- Type[row.names(data), ]
            colnames(Type) <- c("cluster", "Subtype")
            outTab <- data.frame()
            data <- cbind(data, Type)
            for (i in colnames(data[, 1:(ncol(data) - 2)])) {
                rt1 <- data[, c(i, "Subtype")]
                colnames(rt1) <- c("expression", "Subtype")
                outTab <- rbind(outTab, cbind(rt1, gene = i))
            }
            outTab$Subtype <- factor(outTab$Subtype, levels = c("Control", "Disease"))
            summary_stats <- outTab %>% dplyr::group_by(gene, Subtype) %>% dplyr::summarise(mean = mean(expression, na.rm = TRUE), 
                .groups = "drop") %>% dplyr::rename(x = gene)
            p_stats <- outTab %>% dplyr::group_by(gene) %>% dplyr::summarise(p.value = tryCatch(wilcox.test(expression ~ 
                Subtype)$p.value, error = function(e) NA_real_), .groups = "drop") %>% dplyr::rename(x = gene)
            summary_stats <- dplyr::left_join(summary_stats, p_stats, by = "x")
            write.table(summary_stats, file = paste0(outpath, "/BoxPlotdata.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
            p <- ggplot(outTab, aes(x = gene, y = expression)) + geom_boxplot(aes(fill = Subtype), alpha = 0.5, outlier.colour = NA, 
                position = position_dodge(0.85)) + geom_jitter(aes(color = Subtype), alpha = 0.5, position = position_jitterdodge(dodge.width = 0.85)) + 
                scale_fill_manual(values = c("#fac03d", "#ea5b57", "#fac03d")) + scale_color_manual(values = c("#fac03d", 
                "#ea5b57", "#fac03d")) + theme_classic() + labs(x = "", y = "Expression") + theme(axis.text.x = element_text(vjust = 1, 
                size = 12, hjust = 1, colour = "black"), legend.position = "top") + ggpubr::rotate_x_text(45) + ggpubr::stat_compare_means(aes(group = Subtype), 
                symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")), label = "p.signif", 
                method = "wilcox", show.legend = FALSE)
            ggsave(paste0(outpath, "/3.Boxplot.pdf"), width = 12, height = 6, plot = p)
            KeyGene <- gene$KeyGene
            rt <- dplyr::select(results, -c("P-value", "Correlation", "RMSE"))
            immuscore <- function(kg) {
                y <- as.numeric(mat[kg, ])
                do.call(rbind, lapply(colnames(rt), function(x) {
                  dd <- cor.test(as.numeric(rt[, x]), y, method = "pearson")
                  data.frame(gene = kg, immune_cells = x, cor = dd$estimate, p.value = dd$p.value)
                }))
            }
            cordata <- do.call(rbind, lapply(KeyGene, immuscore))
            cordata <- na.omit(cordata)
            write.csv(cordata, file = paste0(outpath, "/HorizonCorr.csv"), row.names = FALSE, quote = FALSE)
            cordata$pstar <- ifelse(cordata$p.value < 0.001, "***", ifelse(cordata$p.value < 0.01, "**", ifelse(cordata$p.value < 
                0.05, "*", "")))
            p <- ggplot(cordata, aes(immune_cells, gene)) + geom_tile(aes(fill = cor)) + geom_text(aes(label = pstar), color = "black", 
                size = 4) + scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c", limit = c(-1, 1), name = paste0("*    p < 0.05\n\n**  p < 0.01\n\n*** p < 0.001\n\nCorrelation")) + 
                labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5, 
                color = "black"), axis.text.y = element_text(size = 8, color = "black"), axis.ticks.y = element_blank(), 
                panel.background = element_blank()) + coord_fixed()
            ggsave(paste0(outpath, "/4.Correlation.pdf"), plot = p, width = 8, height = 4)
        })
        local({
            library(ggplot2)
            library(dplyr)
            library(Hmisc)
            CorMatrix <- function(cor, p) {
                ut <- upper.tri(cor)
                data.frame(row = rownames(cor)[row(cor)[ut]], column = rownames(cor)[col(cor)[ut]], cor = (cor)[ut], p = p[ut])
            }
            capitalize_first <- function(string) {
                words <- strsplit(string, " ")[[1]]
                paste0(toupper(substring(words, 1, 1)), substring(words, 2), collapse = " ")
            }
            base_out <- Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo")
            outpath <- base_out
            keygene.path <- list.files(path = outpath, recursive = TRUE, pattern = "sigGenes\\.txt$", full.names = TRUE)
            KeyGene <- read.table(keygene.path[1], header = TRUE)$KeyGene
            mat <- data.table::fread("./data/expr.txt", data.table = FALSE, header = TRUE)
            mat <- limma::avereps(mat[, -1], ID = mat[, 1])
            if (range(mat)[1] >= 0 & range(mat)[2] - range(mat)[1] > 100) {
                mat <- log2(mat + 1)
            }
            a_1 <- data.frame(mat)
            a_2 <- as.data.frame(t(a_1))
            a_3 <- a_1
            a_3$Id <- rownames(a_3)
            b_1 <- read.table("./data/Immunomodulator_and_chemokines.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE)
            x <- list.files(outpath, pattern = "Immune association", full.names = TRUE)
            if (length(x) == 0) {
                outpath <- paste0(outpath, "/Immune association/")
                dir.create(outpath)
            }
            else {
                outpath <- x[1]
            }
            for (i in seq_along(unique(b_1$type))) {
                b_2 <- b_1[b_1$type == unique(b_1$type)[i], ]
                data1 <- dplyr::inner_join(b_2, a_3, by = "Id")
                data2 <- a_2[, c(KeyGene, data1$Id)]
                res <- rcorr(as.matrix(data2), type = "pearson")
                result_1 <- CorMatrix(res$r, res$P)
                result_2 <- result_1[result_1$row %in% KeyGene, ]
                b_2$column <- b_2$Id
                result_3 <- dplyr::inner_join(result_2, b_2, by = "column")
                result1 <- result_3[, 1:4]
                result1$Regulation <- ifelse(result1$cor > 0, "positive", "negative")
                colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")
                data2 <- result1
                data2$pvalue <- ifelse(data2$pvalue < 0.001, "***", ifelse(data2$pvalue < 0.01, "**", ifelse(data2$pvalue < 
                  0.05, "*", "")))
                data3 <- data2[order(data2$immuneGene, data2$cor), ]
                data4 <- data3[data3$pvalue != "", ]
                data4 <- na.omit(data4)
                pl <- ggplot(data4, aes(immuneGene, gene)) + geom_tile(aes(fill = cor)) + geom_text(aes(label = pvalue), 
                  color = "black", size = 4) + scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c", limit = c(-1, 
                  1), name = paste0("*    p < 0.05\n\n**  p < 0.01\n\n*** p < 0.001\n\nCorrelation")) + labs(x = NULL, y = NULL) + 
                  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, color = "black"), axis.text.y = element_text(size = 12, 
                    color = "black"), axis.ticks.y = element_blank(), panel.background = element_blank()) + coord_fixed()
                ggsave(paste0(outpath, "/", i, ".", capitalize_first(unique(b_1$type)[i]), ".pdf"), width = 4 + length(unique(data4$immuneGene)) * 
                  0.15, height = 4 + length(unique(data4$gene)) * 0.1, plot = pl)
            }
        })
        local({
            library(RcisTarget)
            library(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
            library(igraph)
            library(dplyr)
            data(motifAnnotations_hgnc_v9)
            data(hg19_500bpUpstream_motifRanking_cispbOnly)
            base_out <- Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo")
            outpath <- base_out
            keygene.path <- list.files(path = outpath, recursive = TRUE, pattern = "sigGenes\\.txt$", full.names = TRUE)
            KeyGene <- read.table(keygene.path[1], header = TRUE)$KeyGene
            x <- list.files(outpath, pattern = "Motif", full.names = TRUE)
            if (length(x) == 0) {
                outpath <- paste0(outpath, "/Motif/")
                dir.create(outpath)
            }
            else {
                outpath <- x[1]
            }
            motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
            geneLists <- list(key_gene = intersect(colnames(motifRankings), KeyGene))
            motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings, motifAnnot = motifAnnotations_hgnc_v9)
            motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes) %>% arrange(-NES)
            resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10, ]
            df <- DT::datatable(resultsSubset[, -c("rankAtMax", "TF_lowConf"), with = FALSE], escape = FALSE, filter = "top", 
                options = list(pageLength = 5))
            htmlwidgets::saveWidget(df, paste0(outpath, "/Motif.html"), selfcontained = TRUE)
            onco <- resultsSubset[, c("TF_highConf", "enrichedGenes"), with = FALSE]
            onco$TF_highConf <- gsub("\\s*\\(.*?\\.", "", onco$TF_highConf)
            onco <- onco[onco$TF_highConf != "", ]
            result <- data.frame()
            for (i in seq_len(nrow(onco))) {
                tf <- strsplit(onco[i, ]$TF_highConf, ";")[[1]]
                gene <- strsplit(onco[i, ]$enrichedGenes, ";")[[1]]
                for (j in seq_along(gene)) {
                  result <- rbind(result, data.frame(TF_highConf = tf, enrichedGenes = rep(gene[j], length(tf))))
                }
            }
            result$TF_highConf <- trimws(result$TF_highConf)
            result$enrichedGenes <- trimws(result$enrichedGenes)
            result <- result[!duplicated(result), ]
            write.table(result, paste0(outpath, "/Motif.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
            pdf(paste0(outpath, "/1.TF-Gene-Network.pdf"), width = 7, height = 7)
            g <- graph_from_data_frame(d = result, directed = FALSE)
            V(g)$type <- ifelse(V(g)$name %in% result$TF_highConf, "TF", "Gene")
            V(g)$color <- ifelse(V(g)$type == "TF", "#4DBBD5", "#E64B35")
            V(g)$size <- ifelse(V(g)$type == "TF", 20, 30)
            lay <- layout_with_fr(g)
            plot(g, layout = lay, vertex.label.cex = 0.8, vertex.label.color = "black", edge.color = "grey70", main = "Gene Regulatory Network")
            dev.off()
        })
    })
    local({
        library(igraph)
        library(dplyr)
        base_out <- Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo")
        outpath <- base_out
        keygene.path <- list.files(path = outpath, recursive = TRUE, pattern = "sigGenes\\.txt$", full.names = TRUE)
        KeyGene <- read.table(keygene.path[1], header = TRUE)$KeyGene
        mircode <- data.table::fread("data/mircode_highconsfamilies.txt", header = TRUE, data.table = FALSE)
        key.mircode <- mircode[mircode$gene_symbol %in% KeyGene, ]
        network.table <- data.frame(Gene = key.mircode$gene_symbol, microrna = key.mircode$microrna)
        network.table$microrna <- stringr::str_split(network.table$microrna, "/", simplify = TRUE)[, 1]
        network.table <- network.table[!duplicated(network.table), ]
        x <- list.files(outpath, pattern = "Non-coding RNA networks", full.names = TRUE)
        if (length(x) == 0) {
            outpath <- paste0(outpath, "/Non-coding RNA networks/")
            dir.create(outpath)
        }
        else {
            outpath <- x[1]
        }
        write.table(network.table, paste0(outpath, "/miRNA-network.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
        pdf(paste0(outpath, "/1.miRNA-Gene-Network.pdf"), width = 7, height = 7)
        g <- graph_from_data_frame(d = network.table, directed = FALSE)
        V(g)$type <- ifelse(V(g)$name %in% network.table$microrna, "miRNA", "Gene")
        V(g)$color <- ifelse(V(g)$type == "miRNA", "#4DBBD5", "#E64B35")
        V(g)$size <- ifelse(V(g)$type == "miRNA", 5, 10)
        lay <- layout_with_fr(g)
        plot(g, layout = lay, vertex.label.cex = 0.5, vertex.label.color = "black", edge.color = "grey70", main = "miRNA Regulatory Network")
        dev.off()
    })
    local({
        library(igraph)
        library(dplyr)
        base_out <- Sys.getenv("PIPELINE_BASE_DIR", unset = "./demo")
        outpath <- base_out
        keygene.path <- list.files(path = outpath, recursive = TRUE, pattern = "sigGenes\\.txt$", full.names = TRUE)
        KeyGene <- read.table(keygene.path[1], header = TRUE)$KeyGene
        drug <- data.table::fread("data/interactions.tsv", header = TRUE, data.table = FALSE)
        key.drug <- drug[drug$gene_name %in% KeyGene, ]
        network.table <- data.frame(Gene = key.drug$gene_name, Drug = key.drug$drug_name)
        network.table <- network.table[network.table$Drug != "NULL", ]
        network.table <- network.table[!duplicated(network.table), ]
        x <- list.files(outpath, pattern = "Drug Prediction", full.names = TRUE)
        if (length(x) == 0) {
            outpath <- paste0(outpath, "/Drug Prediction/")
            dir.create(outpath)
        }
        else {
            outpath <- x[1]
        }
        write.table(network.table, paste0(outpath, "/Drug-network.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
        pdf(paste0(outpath, "/1.Drug-Gene-Network.pdf"), width = 7, height = 7)
        g <- graph_from_data_frame(d = network.table, directed = FALSE)
        V(g)$type <- ifelse(V(g)$name %in% network.table$Drug, "Drug", "Gene")
        V(g)$color <- ifelse(V(g)$type == "Drug", "#4DBBD5", "#E64B35")
        V(g)$size <- ifelse(V(g)$type == "Drug", 5, 10)
        lay <- layout_with_fr(g)
        plot(g, layout = lay, vertex.label.cex = 0.8, vertex.label.color = "black", edge.color = "grey70", main = "Drug Prediction")
        dev.off()
    })
    local({
        compact_pdc <- readRDS(file.path(input_dir, "pDC_fullgene_counts.rds"))
        pdc_full <- CreateSeuratObject(counts = compact_pdc$counts, meta.data = compact_pdc$metadata)
        pdc_full <- NormalizeData(pdc_full, verbose = FALSE)
        hallmark_module(outpath = out, SeuratObject = pdc_full, KeyGene = keygenes)
    })
    local({
        cg <- file.path(stage, "legacy_singlecell_module_out", "ClusterGvis")
        st.data <- readRDS(file.path(cg, "Heatmap.plotdata.rds"))
        enrich <- read.csv(file.path(cg, "Enrich.plotdata.csv"), row.names = 1, check.names = FALSE)
        markers <- read.csv(file.path(cg, "TopMarkers.csv"), check.names = FALSE)
        original_ids <- sort(unique(as.integer(as.character(st.data$wide.res$cluster))))
        display_map <- setNames(seq_along(original_ids), original_ids)
        st.data$wide.res$cluster <- unname(display_map[as.character(st.data$wide.res$cluster)])
        st.data$long.res$cluster <- unname(display_map[as.character(st.data$long.res$cluster)])
        cluster_sizes <- table(st.data$wide.res$cluster)
        st.data$long.res$cluster_name <- paste0("cluster ", st.data$long.res$cluster, " (", unname(cluster_sizes[as.character(st.data$long.res$cluster)]), 
            ")")
        enrich_original <- as.integer(sub("^C", "", enrich$group))
        enrich$group <- paste0("C", unname(display_map[as.character(enrich_original)]))
        cluster_order <- seq_along(original_ids)
        expected_groups <- paste0("C", cluster_order)
        missing_groups <- setdiff(expected_groups, unique(enrich$group))
        if (length(missing_groups)) {
            placeholder <- data.frame(group = missing_groups, Description = "No significant GO BP term", pvalue = 1, ratio = 0, 
                row.names = paste0("NO_SIGNIFICANT_GO_", missing_groups), stringsAsFactors = FALSE)
            enrich <- rbind(enrich, placeholder)
        }
        dest <- file.path(out, "4.ClusterGvis", "1.Heatmap.pdf")
        pdf(dest, height = 15, width = 17, onefile = FALSE)
        visCluster(object = st.data, plot.type = "both", column_title_rot = 45, markGenes = unique(markers$gene), markGenes.side = "left", 
            annoTerm.data = enrich, line.side = "left", cluster.order = cluster_order, add.bar = TRUE)
        dev.off()
        writeLines(c(paste0("cluster_order=", paste(cluster_order, collapse = ",")), paste0("original_cluster_ids=", paste(original_ids, 
            collapse = ",")), paste0("missing_enrichment_groups=", paste(missing_groups, collapse = ",")), "placeholder=No significant GO BP term"), 
            file.path(stage, "HEATMAP_ENRICHMENT_COMPATIBILITY_AUDIT.txt"))
    })
    local({
        celltypes <- levels(factor(obj$cellType_1))
        grp <- factor(obj$cellType_1, levels = celltypes)
        indicator <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1, dims = c(length(grp), length(celltypes)), 
            dimnames = list(colnames(obj), celltypes))
        group_n <- as.numeric(table(grp)[celltypes])
        rna_data <- GetAssayData(obj, assay = "RNA", slot = "data")
        rna_counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
        avg <- as.matrix(rna_data %*% indicator)
        pct <- as.matrix((rna_counts > 0) %*% indicator)
        avg <- sweep(avg, 2, group_n, "/")
        pct <- sweep(pct, 2, group_n, "/")
        feature_list <- vector("list", length(celltypes))
        names(feature_list) <- celltypes
        for (j in seq_along(celltypes)) {
            other_max <- if (ncol(avg) == 1L) 
                rep(0, nrow(avg))
            else apply(avg[, -j, drop = FALSE], 1, max)
            specificity <- avg[, j] - other_max
            dominant <- max.col(avg, ties.method = "first") == j
            eligible <- dominant & pct[, j] >= 0.1 & is.finite(specificity)
            ranked <- rownames(avg)[order(specificity, decreasing = TRUE)]
            ranked <- ranked[eligible[ranked]]
            if (length(ranked) < 4L) {
                ranked <- rownames(avg)[order(specificity, decreasing = TRUE)]
            }
            feature_list[[j]] <- head(unique(ranked), 4L)
        }
        feature_list <- lapply(feature_list, function(x) intersect(x, rownames(obj)))
        obj$cellType_1 <- factor(obj$cellType_1, levels = celltypes)
        pdot <- DotPlot(obj, features = feature_list, cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 
            10)) + theme(panel.border = element_rect(colour = "black", fill = NA), axis.text.x = element_text(angle = 90, 
            hjust = 1, vjust = 0.5), legend.position = "top", legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.8, 
            "cm"))
        ggsave(file.path(out, "3.CellAnnotate", "3.Dotplot.pdf"), plot = pdot, width = 3 + sum(lengths(feature_list)) * 0.3, 
            height = 3 + length(celltypes) * 0.5)
        fwrite(rbindlist(lapply(names(feature_list), function(ct) {
            data.table(cell_type = ct, gene = feature_list[[ct]])
        })), "DOTPLOT_STAIRCASE_MARKERS.tsv", sep = "\t")
        my_pal2 <- c("#D4477D", "#4DBBD5", "#D24B27", "#6387C5", "#6E4B9E", "#C10000", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda", 
            "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#00B3F1", "#00538A")
        cell_ratio <- cellRatioPlot(object = obj, sample.name = "group", celltype.name = "cellType_1", flow.curve = 0.5, 
            fill.col = my_pal2) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
        ggsave(file.path(out, "3.CellAnnotate", "3.Group.ratio.plot.pdf"), plot = cell_ratio, width = 7, height = 7)
        cds <- readRDS(file.path(stage, "legacy_singlecell_module_out", "Monocle2", "cds.rds"))
        beam <- BEAM(cds, branch_point = 1, cores = 8, progenitor_method = "duplicate")
        beam <- beam[order(beam$qval, beam$pval), , drop = FALSE]
        beam_genes <- head(rownames(beam)[is.finite(beam$qval)], 50)
        if (length(beam_genes) < 10L) 
            beam_genes <- head(rownames(beam), 50)
        pdf(file.path(monodir, "5.Branch.heatmap.pdf"), width = 8, height = 10)
        plot_genes_branched_heatmap(cds[beam_genes, ], branch_point = 1, num_clusters = 3, cores = 8, use_gene_short_name = TRUE, 
            show_rownames = TRUE, return_heatmap = FALSE)
        dev.off()
        p_branch <- plot_genes_branched_pseudotime(cds[head(beam_genes, 6), ], branch_point = 1, color_by = "Cluster", ncol = 2) + 
            scale_color_manual(values = my_pal2)
        ggsave(file.path(monodir, "6.Branch.genes.pdf"), plot = p_branch, width = 9, height = 8)
        fwrite(data.table(gene = rownames(beam), pvalue = beam$pval, qvalue = beam$qval), "MONOCLE2_BEAM_AUDIT.tsv", sep = "\t")
    })
    local({
        obj <- readRDS(file.path(input_dir, "seurat.annotated.rds"))
        celltypes <- levels(factor(obj$cellType_1))
        obj$cellType_1 <- factor(obj$cellType_1, levels = celltypes)
        markers <- fread("DOTPLOT_STAIRCASE_MARKERS.tsv")
        feature_list <- split(markers$gene, factor(markers$cell_type, levels = celltypes))
        names(feature_list) <- stringr::str_wrap(celltypes, width = 12)
        pdot <- DotPlot(obj, features = feature_list, cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 
            10)) + theme(panel.border = element_rect(colour = "black", fill = NA), axis.text.x = element_text(angle = 90, 
            hjust = 1, vjust = 0.5), strip.text.x = element_text(size = 8, lineheight = 0.9), legend.position = "top", legend.key.height = unit(0.3, 
            "cm"), legend.key.width = unit(0.8, "cm"))
        ggsave(file.path(out, "3.CellAnnotate", "2.Dotplot.pdf"), plot = pdot, width = 3 + sum(lengths(feature_list)) * 0.32, 
            height = 3 + length(celltypes) * 0.5)
        p_umap <- FeaturePlot(obj, features = "Scoring", reduction = "umap", pt.size = 0.45, order = TRUE, raster = TRUE, 
            raster.dpi = c(512, 512)) + viridis::scale_color_viridis() + labs(title = "Peroxisome") + tidydr::theme_dr() + 
            theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15))
        sig_df <- obj@meta.data[, c("cellType_1", "Scoring", "group")]
        colnames(sig_df) <- c("cellType_1", "val", "group")
        sig_df$cellType_1 <- factor(sig_df$cellType_1, levels = celltypes)
        sig_df$group <- factor(sig_df$group, levels = c("Control", "Disease"))
        p_box <- ggplot(sig_df, aes(cellType_1, val, fill = group, color = group)) + geom_boxplot(show.legend = TRUE, width = 0.6, 
            alpha = 0.5, outlier.shape = NA) + scale_fill_manual(values = c(Control = "#398AB9", Disease = "red")) + scale_color_manual(values = c(Control = "#398AB9", 
            Disease = "red")) + theme_classic() + stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif", 
            show.legend = FALSE) + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
            axis.text.y = element_text(vjust = 0.5, size = 12), axis.title.y = element_text(angle = 90, size = 15), axis.title.x = element_text(size = 15)) + 
            xlab("Cell Type") + ylab("Peroxisome Score")
        p_sig <- cowplot::plot_grid(p_umap, p_box, nrow = 1, rel_widths = c(1, 1.5))
        ggsave(file.path(out, "5.Signature", "2.Boxplot of Peroxisome.pdf"), plot = p_sig, width = 12.5, height = 5)
    })
    local({
        mat <- data.table::fread("data/expr.txt", data.table = FALSE, header = TRUE)
        mat <- limma::avereps(mat[, -1, drop = FALSE], ID = mat[[1]])
        if (range(mat)[1] >= 0 && diff(range(mat)) > 100) 
            mat <- log2(mat + 1)
        GSEA_module(mat = mat, outpath = out, species = NA, KeyGene = c("LYL1", "PIM1"), isKEGG = TRUE)
    })
    local({
        obj <- readRDS(file.path(input_dir, "seurat.annotated.rds"))
        pdc_cells <- colnames(obj)[as.character(obj$cellType_2) %in% c("Hsco_pDC", "Lsco_pDC")]
        pdc <- subset(obj, cells = pdc_cells)
        cds <- readRDS(file.path(stage, "legacy_singlecell_module_out", "Monocle2", "cds.rds"))
        common <- intersect(colnames(cds), colnames(pdc))
        cyto_mat <- as.matrix(GetAssayData(pdc, assay = "RNA", slot = "counts")[, common, drop = FALSE])
        cyto <- CytoTRACE::CytoTRACE(mat = cyto_mat, enableFast = FALSE)
        cyto_score <- as.numeric(cyto$CytoTRACE)
        names(cyto_score) <- common
        traj <- data.frame(cell = common, Pseudotime = pData(cds)[common, "Pseudotime"], CytoTRACE = cyto_score[common], 
            MechanismScore = pdc@meta.data[common, "Scoring"], Phenotype = as.character(pdc@meta.data[common, "cellType_2"]), 
            stringsAsFactors = FALSE)
        fwrite(traj, file.path(stage, "MONOCLE2_CYTOTRACE_MECHANISM_PLOTDATA.tsv"), sep = "\t")
        key_full <- readRDS(file.path(input_dir, "keygene_fullcell_lognormalized.rds"))
        key_mat <- Matrix(key_full[keygenes, colnames(obj), drop = FALSE], sparse = TRUE)
        obj[["KeyGene"]] <- CreateAssayObject(data = key_mat)
        DefaultAssay(obj) <- "KeyGene"
        celltypes <- levels(factor(obj$cellType_1))
        obj$cellType_1 <- factor(obj$cellType_1, levels = celltypes)
        p_feat <- FeaturePlot(obj, features = keygenes, reduction = "umap", ncol = 3, order = TRUE, pt.size = 0.25, raster = TRUE, 
            raster.dpi = c(512, 512)) & viridis::scale_color_viridis() & theme_classic()
        ggsave(file.path(out, "16.KeyGene expression abundance", "1.Keygene.FeaturePlot.pdf"), p_feat, width = 15, height = 5, 
            limitsize = FALSE)
        Idents(obj) <- "cellType_1"
        violin_palette <- c("#D4477D", "#4DBBD5", "#D24B27", "#6387C5", "#6E4B9E", "#C10000", "#1E78B4", "#FCBF6E", "#83AD00", 
            "#9ebcda", "#74a9cf", "#fbdf72", "#FF8E00")
        violin_plots <- VlnPlot(obj, features = keygenes, group.by = "cellType_1", pt.size = 0, combine = FALSE, cols = violin_palette)
        violin_plots <- lapply(violin_plots, function(p) {
            p + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), 
                legend.position = "none") + xlab(NULL)
        })
        ggsave(file.path(out, "16.KeyGene expression abundance", "2.Keygene.ViolinPlot.pdf"), wrap_plots(violin_plots, ncol = 1), 
            width = 11, height = 11)
        expr_full <- readRDS(file.path(input_dir, "fullgene_patient_log2cpm.rds"))
        immune_ref <- fread(file.path(stage_data, "Immunomodulator_and_chemokines.txt"), data.table = FALSE)
        cor_long <- function(cor, p) {
            upper <- upper.tri(cor)
            data.frame(row = rownames(cor)[row(cor)[upper]], column = rownames(cor)[col(cor)[upper]], cor = cor[upper], p = p[upper])
        }
        capitalize_first <- function(string) {
            words <- strsplit(string, " ")[[1]]
            paste0(toupper(substring(words, 1, 1)), substring(words, 2), collapse = " ")
        }
        for (i in seq_along(unique(immune_ref$type))) {
            type <- unique(immune_ref$type)[i]
            immune_genes <- intersect(immune_ref$Id[immune_ref$type == type], rownames(expr_full))
            correlation <- Hmisc::rcorr(as.matrix(t(expr_full[c(keygenes, immune_genes), , drop = FALSE])), type = "pearson")
            result <- cor_long(correlation$r, correlation$P)
            result <- result[result$row %in% keygenes & result$column %in% immune_genes, ]
            colnames(result) <- c("gene", "immuneGene", "cor", "pvalue")
            result$star <- ifelse(result$pvalue < 0.001, "***", ifelse(result$pvalue < 0.01, "**", ifelse(result$pvalue < 
                0.05, "*", "")))
            result <- na.omit(result[order(result$immuneGene, result$cor), ])
            result$gene <- factor(result$gene, levels = rev(keygenes))
            p_assoc <- ggplot(result, aes(immuneGene, gene)) + geom_tile(aes(fill = cor), color = "grey92", linewidth = 0.2) + 
                geom_text(aes(label = star), color = "black", size = 3.5) + scale_fill_gradient2(low = "#2b8cbe", mid = "white", 
                high = "#e41a1c", limits = c(-1, 1), name = "*    p < 0.05\n\n  **  p < 0.01\n\n  *** p < 0.001\n\n  Correlation") + 
                labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5, 
                color = "black"), axis.text.y = element_text(size = 12, color = "black"), axis.ticks.y = element_blank(), 
                panel.background = element_blank()) + coord_fixed()
            ggsave(file.path(out, "13.Immune association", paste0(i, ".", capitalize_first(type), ".pdf")), width = max(8, 
                4 + length(immune_genes) * 0.18), height = 5, plot = p_assoc)
        }
    })
    local({
        monocle_out <- file.path(out, "7.Monocle2")
        original_pseudotime <- file.path(stage, "legacy_singlecell_module_out", "Monocle2", "1.Trajectory.Pseudotime.pdf")
        file.copy(original_pseudotime, file.path(monocle_out, "1.Trajectory.Pseudotime.pdf"), overwrite = TRUE)
        cds <- readRDS(file.path(stage, "legacy_singlecell_module_out", "Monocle2", "cds.rds"))
        score_df <- fread(file.path(stage, "MONOCLE2_CYTOTRACE_MECHANISM_PLOTDATA.tsv"), data.table = FALSE)
        rownames(score_df) <- score_df$cell
        cells <- colnames(cds)
        pData(cds)$CytoTRACE_score <- score_df[cells, "CytoTRACE"]
        pData(cds)$Mechanism_score <- score_df[cells, "MechanismScore"]
        make_score_trajectory <- function(color_by, legend_title, plot_title) {
            plot_cell_trajectory(cds, color_by = color_by, cell_size = 0.15) + viridis::scale_color_viridis(name = legend_title) + 
                labs(title = plot_title) + theme(text = element_text(face = "bold"), plot.title = element_text(hjust = 0.5), 
                legend.position = "right")
        }
        p_cyto <- make_score_trajectory("CytoTRACE_score", "CytoTRACE score", "CytoTRACE score")
        p_mechanism <- make_score_trajectory("Mechanism_score", "Score", "Mechanism score")
        ggsave(file.path(monocle_out, "4.Trajectory.CytoTRACE.score.pdf"), p_cyto, width = 4.5, height = 4)
        ggsave(file.path(monocle_out, "5.Trajectory.Peroxisome.mechanism.score.pdf"), p_mechanism, width = 4.5, height = 4)
        fwrite(data.table(item = c("original_pseudotime_restored", "added_score_plots", "cells", "sampling"), value = c("yes", 
            "CytoTRACE;Peroxisome mechanism", ncol(cds), "none")), file.path(stage, "MONOCLE2_ADD_SCORE_TRAJECTORY_AUDIT.tsv"), 
            sep = "\t")
    })
    local({
        monocle_out <- file.path(out, "7.Monocle2")
        cds <- readRDS(file.path(stage, "legacy_singlecell_module_out", "Monocle2", "cds.rds"))
        beam <- fread(file.path(stage, "MONOCLE2_BEAM_AUDIT.tsv"), data.table = FALSE)
        beam <- beam[order(beam$qvalue, beam$pvalue), , drop = FALSE]
        beam_genes <- head(beam$gene[is.finite(beam$qvalue)], 50L)
        beam_genes <- intersect(beam_genes, rownames(cds))
        pdf(file.path(monocle_out, "6.Branch.heatmap.pdf"), width = 8, height = 10)
        plot_genes_branched_heatmap(cds[beam_genes, ], branch_point = 1, num_clusters = 3, cores = 8, use_gene_short_name = TRUE, 
            show_rownames = TRUE, return_heatmap = FALSE)
        dev.off()

    })
    unlink(file.path(out, "3.CellAnnotate", "3.Dotplot.pdf"), force = TRUE)
    unlink(file.path(out, "7.Monocle2", c("4.ShowGene.pseudotime.pdf", "5.Branch.heatmap.pdf", "6.Branch.genes.pdf")), force = TRUE)
    expected_pdf <- readLines(file.path(input_dir, "final_pdf_manifest.txt"))
    observed_files <- sort(list.files(out, pattern = "[.]pdf$", recursive = TRUE, full.names = TRUE))
    observed_pdf <- substring(observed_files, nchar(out) + 2L)
    writeLines(observed_pdf, file.path(stage, "generated_pdf_manifest.txt"))

}
if (sys.nframe() == 0L) run_sepsis_peroxisome_analysis()
