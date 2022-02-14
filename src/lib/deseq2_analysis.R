library("tximport")
library("readr")
library("tximportData")
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("DESeq2"))
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(EnhancedVolcano))
library(ggfittext)
suppressPackageStartupMessages(library(scales))
source("lib/edger_analysis.R")
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Rsamtools))

get_tx2gene = function (gene_annotations) {
  tx2gene_file = paste0(gene_annotations, ".tx2gene.tsv")
  if(!file.exists(tx2gene_file)) {
    tx_db = makeTxDbFromGFF(file = gene_annotations)
    k = keys(tx_db, keytype = "TXNAME")
    tx2gene = select(tx_db, k, "GENEID", "TXNAME")
    write.table(tx2gene, tx2gene_file, row.names = F, quote = F, sep="\t")
  }

  tx2gene = read_delim(paste0(gene_annotations, ".tx2gene.tsv"), delim="\t")
  return(tx2gene)
}


get_samples = function (counts_folder, groups, samples_to_remove = NULL, sample_name_map = NULL, groups_regex = F,
                        allow_skip = F) {
  count_folders = file.path(counts_folder, dir(counts_folder))
  count_folders = count_folders[file.info(count_folders)$isdir]
  stopifnot(length(count_folders) > 0)
  sample_tags = basename(count_folders)

  if(!is.null(sample_name_map)) {
    new_sample_tags = sample_name_map[sample_tags]
    stopifnot(all(!is.na(new_sample_tags)))
    stopifnot(length(sample_tags) == length(sample_name_map))
    specified_order = 1:length(sample_name_map)
    names(specified_order) = names(sample_name_map)
    count_folders[specified_order[sample_tags]] = count_folders
    sample_tags[specified_order[sample_tags]] = new_sample_tags
  }

  if(!is.null(samples_to_remove))
    for(sr in samples_to_remove) {
      remove_idx = grepl(sr, sample_tags)
      stopifnot(any(remove_idx))
      cat("Removing", sum(remove_idx), "samples matching", sr, ":", sample_tags[remove_idx], "\n")
      count_folders = count_folders[!remove_idx]
      sample_tags = sample_tags[!remove_idx]
    }

  quant_files = file.path(count_folders, "quant.sf")
  stopifnot(all(file.exists(quant_files)))
  folder_groups = rep(NA, length(count_folders))
  folder_names = rep(NA, length(count_folders))
  for (g in groups) {
    idx = grepl(g, sample_tags, fixed = !groups_regex)
    stopifnot(any(idx))
    stopifnot(all(is.na(folder_groups[idx])))
    folder_groups[idx] = g
    folder_names[idx] = paste0(g, "_", 1:sum(idx))
  }

  if(!allow_skip) {
    stopifnot(all(!is.na(folder_groups)))
  } else {
    cat(sum(is.na(folder_groups)), "samples removed because they don't match provided groups\n")
    found_idx = !is.na(folder_groups)
    folder_groups = folder_groups[found_idx]
    quant_files = quant_files[found_idx]
    sample_tags = sample_tags[found_idx]
    folder_names = folder_names[found_idx]
  }

  folder_groups = factor(folder_groups, levels = groups)
  return(data.frame(file = quant_files, group = folder_groups, sample_tag = sample_tags,
                    row.names = folder_names, stringsAsFactors = F))
}


import_salmon_into_dds = function (all_samples, tx2gene, min_reads=10, design = ~group) {
  txi = tximport(all_samples$file, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T)
  dds = DESeqDataSetFromTximport(txi,
                                 colData = all_samples,
                                 design = design)

  # keep rows with at least min_reads reads total
  dds = dds[rowSums(counts(dds)) >= min_reads, ]
  dds = DESeq(dds)
  return(dds)
}


vsd_pca_plot = function (vsd, label = F, make_ellipse = T) {
  colnames(vsd) = vsd$sample_tag
  pcaData = plotPCA(vsd, intgroup="group", returnData=TRUE)
  percentVar = round(100 * attr(pcaData, "percentVar"))
  p = ggplot(pcaData, aes(PC1, PC2, color = group, label = name, group = group, fill = group))

  p = p +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()

  if(make_ellipse) {
    ellipse.prob = 0.95
    theta <- seq(-2 * pi, 2 * pi, length = 50)
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(pcaData, "group", function(x) {
        if (nrow(x) <= 2) {
            return(NULL)
        }
        sigma <- var(cbind(x$PC1, x$PC2))
        mu <- c(mean(x$PC1), mean(x$PC2))
        ed <- sqrt(qchisq(ellipse.prob, df = 2))
        data.frame(sweep(circle %*% chol(sigma) * ed, 2,
            mu, FUN = "+"), group = x$group[1])
    })

    names(ell)[1:2] <- c("PC1", "PC2")
    ell$name = ""
    p <- p +
      geom_polygon(data = ell, alpha = 0.3, color = "black")
  }

  if(label)
    p = p + geom_text(vjust = "inward", hjust = "inward")
  else
    p = p + geom_point(size = 3, shape = 21, color = "black")

  return(p)
}


global_picture_deseq = function (dds, file_prefix, contrasts, group_colors, make_ellipse = T) {
  vsd = vst(dds, blind=FALSE)
  sampleDists = dist(t(assay(vsd)))
  sampleDistMatrix = as.matrix(sampleDists)
  rownames(sampleDistMatrix) = vsd$sample_tag
  colnames(sampleDistMatrix) = NULL
  colors = colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
  ph = pheatmap(sampleDistMatrix,
                clustering_distance_rows=sampleDists,
                clustering_distance_cols=sampleDists,
                filename = paste0(file_prefix, "sample_distance.png"),
                main = "Distance between samples",
                col=colors)

  nobg = theme(panel.background = element_blank(), legend.key = element_blank(),
      axis.line = element_line(colour = "black"), text = element_text(size = 10, family = "Arial"))

  p = vsd_pca_plot(vsd, make_ellipse = make_ellipse, label = F) + nobg
  if(!is.null(group_colors))
    p = p + scale_fill_manual(name = "group", values = group_colors) +
      scale_color_manual(name = "group", values = group_colors)

  ggsave(plot = p, paste0(file_prefix, "global_picture.png"), width = 5, height = 5)

  p = vsd_pca_plot(vsd, make_ellipse = make_ellipse, label = T) + nobg
  if(!is.null(group_colors))
    p = p + scale_fill_manual(name = "group", values = group_colors) +
      scale_color_manual(name = "group", values = group_colors)

  ggsave(plot = p, paste0(file_prefix, "global_picture_with_labels.png"), width = 5, height = 5)  

  if(length(contrasts) > 1)
    for (contrast in contrasts) {
      if (length(contrast) == 2 || (length(contrast) == 3 && "name" %in% names(contrast))) {
        condition_group = contrast[1]
        control_group = contrast[2]
        p = vsd_pca_plot(vsd[, vsd$group %in% contrast], label = T, make_ellipse = make_ellipse) + nobg
        if(!is.null(group_colors))
          p = p + scale_fill_manual(name = "group", values = group_colors) +
            scale_color_manual(name = "group", values = group_colors)

        ggsave(plot = p, paste0(file_prefix, condition_group, "_vs_", control_group, "_distribution.png"),
              width=5, height=5)
      }
    }
}

compute_de = function (dds, contrast) {
  #resLFC = lfcShrink(dds, coef=paste0("group_", condition_group, "_vs_", control_group), type="apeglm")
  #resLFC = lfcShrink(dds, contrast = c("group", condition_group, control_group), type = "normal")
  if(length(contrast) == 3 && "name" %in% names(contrast))
    contrast = contrast[1:2]

  if(length(contrast) != 3)
    contrast = c("group", contrast)


  resLFC = lfcShrink(dds, contrast = contrast, type = "ashr")
}


pairwise_venn = function (set1, set2, name1, name2) {
  diagram = venn.diagram(list(set1, set2),
                         category.names = c(name1, name2), filename = NULL,
                         fill=c("blue", "red"), lwd=0,
                         cat.dist=c(0.07, 0.07))

  return(diagram)
}


cleanup_name = function (str) {
  return(gsub("[/ ]", "_", str))
}


write_de_results = function(contrast_results, sample_sizes, sample_unit, alpha, min_abs_logfc,
                                parameters_info,
                                pic_prefix, result_xlsx,
                                min_p_val = NULL, vector_pics = F) {
  unlink(result_xlsx)
  main_parameters_info = list("Main parameters" = " ", "FDR cutoff" = alpha, "Fold change cutoff" = 2 ^ min_abs_logfc)
  main_parameters_info[[paste0("# ", sample_unit, "s")]] = sum(sample_sizes)
  parameters_info = c(main_parameters_info, parameters_info)

  parameters_df = data.frame("X1" = names(parameters_info),
                             "X2" = as.character(parameters_info),
                             stringsAsFactors = F)

  parameters_df["X3"] = ""
  parameters_df["X4"] = ""
  parameters_df["X5"] = ""
  parameters_df["X6"] = ""
  parameters_df["X7"] = ""
  parameters_df[1, 5] = "DE genes summary"
  parameters_df[2, 5] = "Contrast"
  parameters_df[2, 6] = "# DE genes UP"
  parameters_df[2, 7] = "# DE genes DOWN"

  parameters_df[10, 1] = paste0("# ", sample_unit, "s")
  row = 11
  for (g in names(sample_sizes)) {
    parameters_df[row, 1] = g
    parameters_df[row, 2] = sample_sizes[[g]]
    row = row + 1
  }

  de_sheets = list()
  row = 3
  for (contrast_res in contrast_results) {
    contrast = contrast_res$contrast
    if (length(contrast) == 2 || (length(contrast) == 3 && "name" %in% names(contrast))) {
      group_name = "group"
      condition_group = contrast[1]
      control_group = contrast[2]
      if ("name" %in% names(contrast)) {
        contrast_pics_prefix = paste0(pic_prefix, cleanup_name(contrast["name"]))
        title = contrast["name"]
      } else {
        contrast_pics_prefix = paste0(pic_prefix, cleanup_name(condition_group), "_vs_", cleanup_name(control_group))
        title = paste(condition_group, "vs.", control_group)
      }
    } else {
      group_name = contrast[1]
      condition_group = contrast[2]
      control_group = contrast[3]
      contrast_pics_prefix = paste0(pic_prefix, cleanup_name(contrast[1]), "_", cleanup_name(condition_group), "_vs_",
                                    cleanup_name(control_group))

      title = paste0(contrast[1], " ", condition_group, " vs. ", control_group)
    }

    message("Control is: ", control_group, ", condition is: ", condition_group)

    resLFC = contrast_res$res
    selection = (resLFC$padj < alpha) & (abs(resLFC$log2FoldChange) > min_abs_logfc)

    selection[is.na(selection)] = F
    names(selection) = rownames(resLFC)
    cat(sum(selection, na.rm = T), "DE genes in", title, "\n")

    ylab = bquote(~-Log[10]~italic(P[adj]))
    selection_up = selection & (resLFC$log2FoldChange > 0)
    selection_down = selection & (resLFC$log2FoldChange < 0)
    num_up = sum(selection_up, na.rm = T)
    num_down = sum(selection_down, na.rm = T)
    lim = max(abs(resLFC$log2FoldChange[selection]))
    xlim = c(-lim, lim)
    if(!is.null(min_p_val)) {
        resLFC$padj = pmax(min_p_val, resLFC$padj)
    }

    EnhancedVolcano(resLFC, x="log2FoldChange", y="padj", lab=resLFC$gene_name,
        ylab=ylab,
        xlab = bquote(~Log[2]~italic("Fold change")),
        title=title, FCcutoff=min_abs_logfc, pCutoff=alpha, gridlines.major=F, gridlines.minor=F,
                    legendPosition = 'none',
        subtitle=NULL,
        caption=paste0(num_up, ' up, ', num_down, " down"),
        labSize=3,
        drawConnectors=F,
        xlim=xlim)

    ggsave(filename = paste0(contrast_pics_prefix, "_volcano.png"), width = 10, height = 10)

    EnhancedVolcano(resLFC, x="log2FoldChange", y="padj", lab=resLFC$gene_name,
        ylab=ylab,
        xlab = bquote(~Log[2]~italic("Fold change")),
        title="", FCcutoff=min_abs_logfc, pCutoff=alpha, gridlines.major=F, gridlines.minor=F,
                    legendPosition = 'none',
        subtitle=NULL,
        caption="",
        labSize=3,
        drawConnectors=F,
        xlim=xlim)

    ggsave(filename = paste0(contrast_pics_prefix, "_volcano.pdf"), width = 10, height = 10)

    resLFC_sig = resLFC[!is.na(resLFC$padj),]
    resLFC_sig = resLFC_sig[(resLFC_sig$padj < alpha) & (abs(resLFC_sig$log2FoldChange) >= min_abs_logfc),]
    keyvals <- ifelse(resLFC_sig$log2FoldChange < 0, 'royalblue', 'red')
    names(keyvals)[keyvals == 'royalblue'] = 'down'
    names(keyvals)[keyvals == 'red'] = 'up'
    EnhancedVolcano(resLFC_sig, x="log2FoldChange", y="padj", lab=NA,
                    ylab = ylab,
                    xlab = bquote(~Log[2]~italic("Fold change")),
                    title = title, FCcutoff = min_abs_logfc, pCutoff = alpha, gridlines.major = F, gridlines.minor = F,
                    legendPosition = 'none',
                    subtitle=NULL,
                    caption=paste0(num_up, ' up, ', num_down, " down"),
                    labSize=3,
                    drawConnectors=F,
                    xlim=xlim,
                    cutoffLineType = "blank",
                    colCustom = keyvals)

    ggsave(filename=paste0(contrast_pics_prefix, "_volcano_small.png"), width=5, height=5)

    keyvals <- ifelse(resLFC_sig$log2FoldChange < 0, 'darkgray', 'black')
    names(keyvals)[keyvals == 'darkgray'] = 'down'
    names(keyvals)[keyvals == 'black'] = 'up'
    EnhancedVolcano(resLFC_sig, x="log2FoldChange", y="padj", lab=NA,
                    ylab = ylab,
                    xlab = bquote(~Log[2]~italic("Fold change")),
                    title = title, FCcutoff = min_abs_logfc, pCutoff = alpha, gridlines.major = F, gridlines.minor = F,
                    legendPosition = 'none',
                    subtitle=NULL,
                    caption = paste0(num_up, ' up, ', num_down, " down"),
                    labSize = 3,
                    drawConnectors=F,
                    xlim=xlim,
                    cutoffLineType = "blank",
                    colCustom = keyvals, colAlpha = 1)

    ggsave(filename=paste0(contrast_pics_prefix, "_volcano_small_bw.png"), width=5, height=5)

    long_title = paste(title, "DE genes")

    selected_res = resLFC[selection, c("gene_name", "log2FoldChange", "padj")]
    selected_res = selected_res[order(selected_res$padj),]
    de_sheets[[long_title]] = selected_res
    parameters_df[row, 5] = title
    parameters_df[row, 6] = num_up
    parameters_df[row, 7] = num_down
    row = row + 1

    resLFC$genes = rownames(resLFC)
  }

  write.xlsx2(parameters_df, file = result_xlsx, sheetName = "Summary", append = T, row.names = F, col.names = F)
  for(long_title in names(de_sheets)) {
    sheet_title = gsub(":", "|", long_title)
    write.xlsx2(de_sheets[[long_title]], file = result_xlsx, sheetName = cleanup_name(sheet_title), append = T,
                row.names = F)

  }
}


de_add_go_enrichment = function(contrast_results, organism, alpha, min_abs_logfc, go_enrichment_cutoff, pic_prefix,
                                result_xlsx,
                                gene_ids_file = NULL, go_top_n = NULL) {
  if (!is.null(gene_ids_file))
    gene_id_map_arg = paste("--gene_name_map", gene_ids_file)
  else
    gene_id_map_arg = ""

  for (contrast_res in contrast_results) {
    contrast = contrast_res$contrast
    if (length(contrast) == 2 || (length(contrast) == 3 && "name" %in% names(contrast))) {
      group_name = "group"
      condition_group = contrast[1]
      control_group = contrast[2]
      if ("name" %in% names(contrast)) {
        contrast_pics_prefix = paste0(pic_prefix, cleanup_name(contrast["name"]))
        title = contrast["name"]
      } else {
        contrast_pics_prefix = paste0(pic_prefix, cleanup_name(condition_group), "_vs_", cleanup_name(control_group))
        title = paste(condition_group, "vs.", control_group)
      }
    } else {
      group_name = contrast[1]
      condition_group = contrast[2]
      control_group = contrast[3]
      contrast_pics_prefix = paste0(pic_prefix, cleanup_name(contrast[1]), "_", cleanup_name(condition_group), "_vs_",
                                    cleanup_name(control_group))

      title = paste0(contrast[1], " ", condition_group, " vs. ", control_group)
    }

    resLFC = contrast_res$res

    selection = (resLFC$padj < alpha) & (abs(resLFC$log2FoldChange) > min_abs_logfc)
    selection[is.na(selection)] = F
    names(selection) = rownames(resLFC)

    contrast_pics_prefix = paste0(pic_prefix, gsub(" ", "_", cleanup_name(title)), "_GO")
    if(file.exists(contrast_pics_prefix))
      unlink(contrast_pics_prefix, recursive = T)

    if(is.null(go_top_n))
      top_n_option = ""
    else
      top_n_option = paste0(" --top_n ", go_top_n, " ")

    contrast_ws_name = cleanup_name(substr(title, 1, 18))
    command_line =
    output = system(paste0("python -m lib.pathway_analysis ", organism, " '", contrast_ws_name, " down|GO' ",
                           result_xlsx,
                           " '", contrast_pics_prefix, "' --enrichment_cutoff ", go_enrichment_cutoff, top_n_option,
                           " ", gene_id_map_arg),
                    input=resLFC[selection & (resLFC$log2FoldChange < 0), "genes"], intern = T)

    output = system(paste0("python -m lib.pathway_analysis ", organism, " '", contrast_ws_name, " up|GO' ",
                           result_xlsx,
                           " '", contrast_pics_prefix, "' --enrichment_cutoff ", go_enrichment_cutoff, top_n_option,
                           " ", gene_id_map_arg),
                    input=resLFC[selection & (resLFC$log2FoldChange > 0), "genes"], intern = T)
  }
}


write_results_deseq = function (dds, organism, contrasts, alpha, min_abs_logfc,
                                go_enrichment_cutoff,
                                pic_prefix, result_xlsx,
                                gene_ids_file = NULL, min_fpkm = 0, min_p_val = NULL, go_top_n = NULL,
                                group_colors = NULL, make_ellipse = T, compare_contrasts = F, heatmap_margins = c(8, 8),
                                vector_pics = F) {
  stopifnot(is.null(gene_ids_file) || file.exists(gene_ids_file))
  global_picture_deseq(dds, pics_prefix, contrasts, group_colors, make_ellipse = make_ellipse)

  fpkms = fpkm(dds)
  colnames(fpkms) = dds@colData$sample_tag
  genes_expressed_high_enough = rownames(fpkms)[rowSums(fpkms > min_fpkm) > 0]
  message(length(genes_expressed_high_enough), "/", nrow(fpkms), " genes with at least one sample having FPKM > ",
    min_fpkm)

  contrast_results = list()
  de_sheets = list()
  row = 3
  for (contrast in contrasts) {
    if (length(contrast) == 2 || (length(contrast) == 3 && "name" %in% names(contrast))) {
      group_name = "group"
      condition_group = contrast[1]
      control_group = contrast[2]
      if ("name" %in% names(contrast)) {
        contrast_pics_prefix = paste0(pic_prefix, cleanup_name(contrast["name"]))
        title = contrast["name"]
      } else {
        contrast_pics_prefix = paste0(pic_prefix, cleanup_name(condition_group), "_vs_", cleanup_name(control_group))
        title = paste(condition_group, "vs.", control_group)
      }
    } else {
      group_name = contrast[1]
      condition_group = contrast[2]
      control_group = contrast[3]
      contrast_pics_prefix = paste0(pic_prefix, cleanup_name(contrast[1]), "_", cleanup_name(condition_group), "_vs_",
                                    cleanup_name(control_group))

      title = paste0(contrast[1], " ", condition_group, " vs. ", control_group)
    }

    resLFC = compute_de(dds, contrast)
    if (!is.null(gene_ids_file)) {
      gene_ids = read.table(gene_ids_file, sep="\t", header=T, stringsAsFactors = F)
      rownames(gene_ids) = gene_ids$id
      resLFC$gene_name = gene_ids[rownames(resLFC), "name"]
    }
    else {
      resLFC$gene_name = rownames(resLFC)
      gene_ids = NULL
    }

    resLFC = resLFC[rownames(resLFC) %in% genes_expressed_high_enough, ]

    selection = (resLFC$padj < alpha) & (abs(resLFC$log2FoldChange) > min_abs_logfc)

    selection[is.na(selection)] = F
    names(selection) = rownames(resLFC)
    cat(sum(selection, na.rm = T), "DE genes in", title, "\n")

    selection_up = selection & (resLFC$log2FoldChange > 0)
    selection_down = selection & (resLFC$log2FoldChange < 0)
    if(!is.null(min_p_val)) {
        resLFC$padj = pmax(min_p_val, resLFC$padj)
    }

    nobg = theme(panel.background = element_blank(), legend.key = element_blank(),
        axis.line = element_line(colour = "black"), text = element_text(size = 20, family = "Arial"))

    control_fpkm = rowMeans(fpkms[rownames(resLFC), dds[[group_name]] == control_group])
    condition_fpkm = rowMeans(fpkms[rownames(resLFC), dds[[group_name]] == condition_group])

    plot_data = data.frame(
        control_fpkm = control_fpkm,
        condition_fpkm = condition_fpkm,
        name = resLFC$gene_name,
        id = rownames(resLFC)
    )

    max_limit = max(c(plot_data$control_fpkm, plot_data$condition_fpkm))
    plot = ggplot(plot_data[plot_data$id %in% genes_expressed_high_enough, ],
                  aes(x = control_fpkm + 1, y = condition_fpkm + 1)) +
        geom_point(color = "gray") +
        geom_point(data = plot_data[selection_up,], color = "red") +
        geom_point(data = plot_data[selection_down,], color = "royalblue") +
        labs(x = paste0("FPKM(", control_group, ")"), y = paste0("FPKM(", condition_group, ")")) +
        nobg +
        scale_x_log10(limits = c(1, max_limit), breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_y_log10(limits = c(1, max_limit), breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))

    ggsave(plot = plot, filename = paste0(contrast_pics_prefix, "_scatter.png"), width=5, height=5)

    plot = ggplot(plot_data[plot_data$id %in% genes_expressed_high_enough, ],
                  aes(x = control_fpkm + 1, y = condition_fpkm + 1)) +
        geom_point(color = "lightgray") +
        geom_point(data = plot_data[selection_up,], color = "black", shape = 21, fill = "black") +
        geom_point(data = plot_data[selection_down,], color = "black", shape = 21, fill = "darkgray") +
        labs(x = paste0("FPKM(", control_group, ")"), y = paste0("FPKM(", condition_group, ")")) +
        nobg +
        scale_x_log10(limits = c(1, max_limit), breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_y_log10(limits = c(1, max_limit), breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))

    ggsave(plot = plot, filename = paste0(contrast_pics_prefix, "_scatter_bw.png"), width=5, height=5)
    #print(plot)
                      
    long_title = paste(title, "DE genes")
    if(sum(selection) > 1) {
      if(length(contrasts) > 1) {
          par(cex.main=0.5)
          selected_fpkms = fpkms[selection[rownames(fpkms)],]
          rownames(selected_fpkms) = gene_ids[rownames(selected_fpkms), "name"]
          plot_rpkms_heatmap(selected_fpkms, long_title,
                             paste0(contrast_pics_prefix, "_de_genes_heatmap_all_samples.png"), margins = heatmap_margins,
                             cluster_cols = F)
      }

      par(cex.main=0.5)
      selected_fpkms = fpkms[selection[rownames(fpkms)], dds[[group_name]] %in% contrast]
      rownames(selected_fpkms) = gene_ids[rownames(selected_fpkms), "name"]
      plot_rpkms_heatmap(selected_fpkms, long_title,
                         paste0(contrast_pics_prefix, "_de_genes_heatmap.png"), margins = heatmap_margins, cexCol = 1,
                         cluster_cols = F)
    }

    resLFC$genes = rownames(resLFC)
    contrast_res = list(contrast = contrast, res = resLFC, title = title)
    contrast_results[[length(contrast_results) + 1]] = contrast_res
  }

  sample_sizes = NULL
  sn = unique(dds$group)
  for (g in sn)
    sample_sizes = c(sample_sizes, sum(dds$group == g))

  names(sample_sizes) = sn

  write_de_results(contrast_results, sample_sizes, "sample", alpha, min_abs_logfc,
                   list("GO enrichment cutoff" = go_enrichment_cutoff), pic_prefix, result_xlsx, min_p_val = min_p_val,
                   vector_pics = vector_pics)

  fpkms_result = make_rpkms_table(fpkms, gene_ids)
  for (contrast_res in contrast_results) {
    fpkms_result = append_de_result_to_rpkms(fpkms_result, contrast_res$res, contrast_res$title,
                                             c("log2FoldChange", "padj"))
  }

  fpkms_result = fpkms_result[order(fpkms_result[, ncol(fpkms_result)]),]

  write.xlsx2(fpkms_result, file = result_xlsx, sheetName = "FPKM", append = T, row.names = F)
  de_add_go_enrichment(contrast_results, organism, alpha, min_abs_logfc, go_enrichment_cutoff, pic_prefix, result_xlsx,
                       gene_ids_file, go_top_n)

  if(compare_contrasts) {
    stopifnot(length(contrast_results) == 2)
    res1 = contrast_results[[1]]$res
    res2 = contrast_results[[2]]$res
    name1 = contrast_results[[1]]$title
    name1_clean = gsub("\\.", "",  cleanup_name(name1))
    name2 = contrast_results[[2]]$title
    name2_clean = gsub("\\.", "",  cleanup_name(name2))
    selection1 = (res1$padj < alpha) & (abs(res1$log2FoldChange) > min_abs_logfc)
    selection1[is.na(selection1)] = F
    names(selection1) = rownames(res1)

    selection2 = (res2$padj < alpha) & (abs(res2$log2FoldChange) > min_abs_logfc)
    selection2[is.na(selection2)] = F
    names(selection2) = rownames(res2)

    diag_up = pairwise_venn(rownames(res1[selection1 & (res1$log2FoldChange > 0),]),
                            rownames(res2[selection2 & (res2$log2FoldChange > 0),]),
                            paste0(name1, " up"), paste0(name2, " up"))

    diag_down = pairwise_venn(rownames(res1[selection1 & (res1$log2FoldChange < 0),]),
                              rownames(res2[selection2 & (res2$log2FoldChange < 0),]),
                              paste0(name1, " down"), paste0(name2, " down"))

    png(paste0(pic_prefix, name1_clean, "_compared_to_", name2_clean, "_venn.png"), width=2000, height=1500, res=250)
    grid.arrange(grobs=list(grobTree(diag_down), grobTree(diag_up)), respect=T,
                 top=paste0("Intersection of DE genes between\n", name1, " and\n", name2))

    dev.off()

    both_up = intersect(rownames(res1[selection1 & (res1$log2FoldChange > 0),]), rownames(res2[selection2 & (res2$log2FoldChange > 0),]))
    comparison_table_up = cbind(res1[both_up, c("genes", "gene_name", "log2FoldChange")], res2[both_up, "log2FoldChange", drop=F])
    if(nrow(comparison_table_up) > 0)
      comparison_table_up$direction = "up"

    both_down = intersect(rownames(res1[selection1 & (res1$log2FoldChange < 0),]), rownames(res2[selection2 & (res2$log2FoldChange < 0),]))
    comparison_table_down = cbind(res1[both_down, c("genes", "gene_name", "log2FoldChange")], res2[both_down, "log2FoldChange", drop=F])
    if(nrow(comparison_table_down) > 0)
      comparison_table_down$direction = "down"

    comparison_table = rbind(comparison_table_up, comparison_table_down)
    colnames(comparison_table)[3:4] = c(paste(name1, "log2FC"), paste(name2, "log2FC"))
    write.xlsx(comparison_table, file = result_xlsx, sheetName = "Contrasts intersection", append = T, row.names = F)
  }
}


calc_gene_lengths = function(GTFfile, genome) {
  #Load the annotation and reduce it
  GTF <- import.gff(GTFfile, format="gtf", genome=genome, feature.type="exon")
  grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
  reducedGTF <- unlist(grl, use.names=T)
  elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

  #Add the GC numbers
  elementMetadata(reducedGTF)$widths <- width(reducedGTF)

  #Create a list of the ensembl_id/GC/length
  calc_GC_length <- function(x) {
      width = sum(elementMetadata(x)$widths)
      return(width)
  }
  output <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length)
  return(output)
}