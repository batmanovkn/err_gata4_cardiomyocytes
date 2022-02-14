suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(gsubfn))
suppressPackageStartupMessages(library("RColorBrewer"))


group_analysis = function(dg, groups) {
  design_matrix = model.matrix(~0 + groups)
  dg = estimateDisp(dg, design=design_matrix)
  return(list(fit=glmFit(dg, design_matrix), design=design_matrix))
}


pairwise_analysis = function(group_res, contrast1, contrast2) {
  contrasts = makeContrasts(contrasts=paste0("groups", contrast2, " - groups", contrast1), levels=group_res$design)
  lrt = glmLRT(group_res$fit, contrast=contrasts)
  res = topTags(lrt, n=100000, p.value=1)$table
  res$genes = as.character(res$genes)
  return(res)
}


load_counts = function(counts_data_file, additional_columns=NULL, start_column = 7, gene_id_column = "Geneid") {
    counts_data = read.csv(counts_data_file, sep="\t", comment.char="#", check.names=F, stringsAsFactors=F)
    if (is.null(additional_columns)) {
        end = ncol(counts_data)
    } else {
        end = ncol(counts_data) - additional_columns
    }

    counts = as.matrix(counts_data[, start_column:end])
    sample_names = colnames(counts)
    sample_names = gsub(".*/", "", sample_names)
    sample_names = gsub("\\.bam", "", sample_names)
    colnames(counts) = sample_names
    rownames(counts) = as.character(counts_data[, gene_id_column])
    if(is.null(additional_columns))
        return(list(counts=counts, gene_lengths=counts_data$Length))

    return(list(counts=counts, gene_lengths=counts_data$Length, additional=counts_data[, (end + 1):ncol(counts_data)]))
}


global_picture = function(counts, gene_lengths, sample_groups, min_cpm=1, min_columns_detected=2, plot_names=F,
                          legend_place="bottomright",
                          main_title="Sample distribution based on gene expression", samples_to_plot=NULL,
                          group_colors=NULL) {
    cat("Counts data:", nrow(counts), "genes,", ncol(counts), "samples\n")
    dg_list = DGEList(counts=counts, genes=rownames(counts))
    
    counts_per_million = cpm(dg_list)
    count_check = counts_per_million > min_cpm
    keep = which(rowSums(count_check) >= min_columns_detected)
    dg_list = dg_list[keep,]
    gene_lengths = gene_lengths[keep]
    cat(nrow(dg_list$counts), "genes after filtering out low counts\n")
    
    dg_list = calcNormFactors(dg_list, method="TMM")    
    if (is.null(group_colors)) {
        all_groups = unique(sample_groups)
        group_colors = brewer.pal(n=length(all_groups), name="Set1")[1:length(all_groups)]
        names(group_colors) = all_groups
    }
    
    if (is.null(samples_to_plot)) {
        dg_list_to_plot = dg_list
    } else {
        dg_list_to_plot = dg_list[, samples_to_plot]
        sample_groups = sample_groups[samples_to_plot]
    }
    
    sample_colors = rep("black", length(sample_groups))
    for (i in 1:length(group_colors))
        sample_colors[sample_groups == names(group_colors)[i]] = group_colors[i]
    
    if(plot_names) {
        plotMDS(dg_list_to_plot, col=sample_colors)
    } else {
        plotMDS(dg_list_to_plot, col=sample_colors, pch=16)
        legend(legend_place, legend=unique(sample_groups), col=group_colors[unique(sample_groups)], pch=16)
    }
    
    
    title(main_title, sub="distances between points ≈ typical logFC between samples")
    rpkms = rpkm(dg_list, gene_lengths)
    
    return(list(dg_list=dg_list, rpkms=rpkms))
}


simple_analysis = function(dg_list, sample_groups, pairs, gene_ids=NULL) {
    group_res = group_analysis(dg_list, sample_groups)
    results = list()
    for (name in names(pairs)) {
        gg = pairs[[name]]
        g1 = gg[1]
        g2 = gg[2]
        cat("Analyzing", name, ":", g1, "vs.", g2, "\n")
        stopifnot(g1 %in% sample_groups)
        stopifnot(g2 %in% sample_groups)
        pa = pairwise_analysis(group_res, g1, g2)
        if (!is.null(gene_ids))
            pa$gene_name = gene_ids[as.character(pa$genes), "name"]
        else
            pa$gene_name = as.character(pa$genes)
        
        results[[name]] = pa
    }
    
    return(results)
}


plot_rpkms_heatmap_old = function(rpkms_to_plot, title, file_name, cluster_cols=T, colorbar=F) {
    pheatmap::pheatmap(rpkms_to_plot, scale="row", show_rownames=nrow(rpkms_to_plot) < 30, treeheight_row=0, colorRampPalette(c("navy", "white", "firebrick3"))(50),
                    angle_col=90, legend=colorbar, main=title, cluster_cols=cluster_cols,
                    filename=file_name)
}


plot_rpkms_heatmap = function(rpkms_to_plot, title, file_name, cluster_cols=T, colorbar=F, main_font_scale=1,
                              cluster_rows=T, ...) {
    if(nrow(rpkms_to_plot) < 30) {
        rownames = rownames(rpkms_to_plot)
    } else {
        rownames = NA
    }
    
    if(cluster_cols) {
        dendrogram = "column"
    }
    else {
        dendrogram = "none"
    }

    if(!is.null(file_name))
        png(file_name, width=2000, height=2000, res=300)

    par(cex.main=main_font_scale)
    gplots::heatmap.2(t(scale(t(rpkms_to_plot))), scale="row", labRow=rownames, dendrogram=dendrogram,
                      col=colorRampPalette(c("navy", "white", "firebrick3"))(50),
                      main=title, Rowv=cluster_rows, Colv=cluster_cols, trace="none", density.info="none", keysize=1,
                      key=colorbar,
                      ...)
    
    if(!is.null(file_name))
        dev.off()
}


make_rpkms_table = function (rpkms, gene_ids=NULL) {
    rpkms_result = as.data.frame(rpkms)
    if (is.null(gene_ids)) {
        rpkms_result$gene_id = rownames(rpkms_result)
        rpkms_result = rpkms_result[, c(ncol(rpkms_result), 1:(ncol(rpkms_result) - 1))]
    } else {
        rpkms_result$gene_name = gene_ids[rownames(rpkms_result), "name"]
        rpkms_result$gene_id = rownames(rpkms_result)
        rpkms_result = rpkms_result[, c(ncol(rpkms_result) - 1, ncol(rpkms_result), 1:(ncol(rpkms_result) - 2))]
    }

    rpkms_sorted = t(apply(rpkms, 1, sort))
    rpkms_result$outlier = F
    second_nonzero = rpkms_sorted[, ncol(rpkms_sorted) - 1] != 0
    rpkms_result$outlier[second_nonzero] =
        rpkms_sorted[second_nonzero, ncol(rpkms_sorted)] / rpkms_sorted[second_nonzero, ncol(rpkms_sorted) - 1] > 3

    return(rpkms_result)
}


append_de_result_to_rpkms = function (rpkms_result, res, contrast, columns_to_append, id_column="genes") {
    for(col in columns_to_append)
        rpkms_result[as.character(res[, id_column]), paste(contrast, col)] = res[, col]

    return(rpkms_result)
}


write_results = function(results, rpkms, min_abs_logfc, alpha, suffix="", gene_ids=NULL, data_prefix="data/",
                         add_raw_pvalue=F, add_title=F, samples_to_plot=NULL, show_colorbar=F, volcano_genes=NULL,
                         p_value_column="FDR",
                         draw_connectors=T, heatmap_font_scale=1, heatmap_margins=c(8, 8), gene_features=NULL,
                         cluster_cols = T,
                         pics_prefix = "pics/", vector_pics = F) {
    result_tables_file = file.path(data_prefix, paste0("result_tables", suffix, ".xlsx"))
    for (contrast in names(results)) {
        res = results[[contrast]]
        write.table(res, file=file.path(data_prefix, paste0("de_genes_", contrast, suffix, ".tsv")),
                    sep="\t", quote=F, row.names=F)
        
        selection = (res[, p_value_column] < alpha) & (abs(res$logFC) > min_abs_logfc)
        rpkms_de = rpkms[as.character(res$genes[selection]), ]
        if(!is.null(gene_ids))
            rownames(rpkms_de) = gene_ids[rownames(rpkms_de), "name"]
        
        cat(contrast, ":", sum(selection), "DE genes\n")
        if (add_title) {
            title = paste0(contrast, " DE genes")
        } else {
            title = ""
        }
        
        if (is.null(volcano_genes)) {
            transcript_lab_size = 3
        } else {
            transcript_lab_size = 5            
        }

        if(p_value_column == "FDR") {
            ylab = bquote(~-Log[10]~italic(P[adj]))
        } else {
            ylab = bquote(~-Log[10]~italic(P))
        }

        plot = EnhancedVolcano(res, x="logFC", y=p_value_column, lab=res$gene_name,
            ylab=ylab,
            xlab = bquote(~Log[2]~italic("Fold change")),
            title=title, FCcutoff=min_abs_logfc, pCutoff=alpha, gridlines.major=F, gridlines.minor=F, legendVisible=F,
            subtitle=NULL,
            caption=paste0(sum(selection & res$logFC > 0), ' up, ', sum(selection & res$logFC < 0), " down"),
            pLabellingCutoff=alpha * 0.1, selectLab=volcano_genes, labSize=transcript_lab_size,
            drawConnectors=draw_connectors)

        if (vector_pics) {
            ggsave(plot, filename=file.path(pics_prefix, paste0("volcano_", contrast, suffix, ".pdf")), width=10,
                height=10)
        } else {
            ggsave(plot, filename=file.path(pics_prefix, paste0("volcano_", contrast, suffix, ".png")), width=10,
                height=10)
        }

        if(sum(selection) == 0) {
            unlink(file.path(pics_prefix, paste0("de_heatmap_", contrast, suffix, ".png")))
        }
        else {
            if (add_title) {
                title = paste0("RPKM of DE genes for ", contrast, " (", nrow(rpkms_de), ")")
            } else {
                title = ""
            }
            
            if(!is.null(samples_to_plot))
                rpkms_de = rpkms_de[, samples_to_plot]
            
            plot_rpkms_heatmap(rpkms_de, title, file.path(pics_prefix, paste0("de_heatmap_", contrast, suffix, ".png")),
                               colorbar=show_colorbar, main_font_scale = heatmap_font_scale, margins = heatmap_margins,
                               cluster_cols = cluster_cols)
        }

        to_write = res[selection, c("genes", "gene_name", "logFC", p_value_column)]
        if (!is.null(gene_features))
            to_write = merge(to_write, gene_features, by.x = "genes", by.y = "row.names", all.x = T)

        to_write$genes = NULL
        to_write = to_write[order(to_write[, p_value_column]),]

        write.xlsx2(to_write, file=result_tables_file,
            sheetName=contrast, append=T, row.names=F)
    }

    rpkms_result = make_rpkms_table(rpkms, gene_ids)

    columns_to_append = c("logFC", p_value_column)
    if(add_raw_pvalue)
        columns_to_append = c(columns_to_append, "PValue")

    for (contrast in names(results))
        rpkms_result = append_de_result_to_rpkms(rpkms_result, results[[contrast]], contrast, columns_to_append)

    write.xlsx2(rpkms_result, file=result_tables_file, sheetName="RPKM", append=T, row.names=F)
}


update_xls_sheet = function (new_data, file, sheet_name) {
    wb = loadWorkbook(file)
    writeDataTable(wb, sheet = sheet_name, new_data, rowNames = F)
    saveWorkbook(wb, file, overwrite = T)
}


make_gene_plotting_data = function(rpkms, genes_to_plot, sample_groups, gene_ids=NULL, by_names=F) {
    if(by_names) {
        given_gene_names = genes_to_plot
        genes_to_plot = gene_ids[gene_ids$name %in% genes_to_plot, "id"]
    }

    stopifnot(all(!is.na(genes_to_plot)))
    if(!all(genes_to_plot %in% rownames(rpkms))) {
        if(by_names)
            stop("These genes are not in RPKMS: ",
                 paste(gene_ids[genes_to_plot[!(genes_to_plot %in% rownames(rpkms))], "name"], collapse = ", "))

        stop("These genes are not in RPKMS: ",
            paste(genes_to_plot[!(genes_to_plot %in% rownames(rpkms))], collapse = ", "))
    }

    df = as.data.frame(rpkms[genes_to_plot, , drop=F])
    if(!is.null(gene_ids))
        df$gene_name = gene_ids[rownames(df), "name"]
    else
        df$gene_name = rownames(df)
    
    if(by_names)
        df$gene_name = factor(df$gene_name, levels=given_gene_names)
    
    df = melt(df, id.vars="gene_name")
    df$group = NA
    sample_names = colnames(rpkms)
    for (i in 1:length(sample_names))
        df$group[df$variable == sample_names[i]] = sample_groups[i]
    
    return(df)
}


plot_genes = function(gene_plotting_data) {
    plot = ggplot(gene_plotting_data, aes(x=group, y=value)) + geom_boxplot(position="identity") +
        labs(y="RPKM", color="Gene") +
        theme(
            panel.border=element_blank(),
            panel.grid=element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1)
        ) + 
        facet_wrap("gene_name", scales="free")
    
    return(plot)
}
