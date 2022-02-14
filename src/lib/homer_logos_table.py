import os
import argparse
import re

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import logomaker
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = "white"
from scipy.cluster.hierarchy import dendrogram, linkage

from . import misc


def read_motifs_file(motifs_file):
    header = None
    data = []
    result = []
    with open(motifs_file) as f:
        for line in f:
            line = line.replace("\n", "")
            if line == "":
                continue

            if line[0] == ">":
                if header is not None:
                    result.append((header, np.array(data, dtype=float)))

                data = []
                header = line[1:]
            else:
                data.append(line.split("\t"))

    if header is not None:
        result.append((header, np.array(data, dtype=float)))

    return result


def convert_motif_matrix(motif_matrix, reverse_complement):
    if reverse_complement:
        motif_matrix = np.flipud(np.fliplr(motif_matrix))

    ic = 2 + (motif_matrix * np.log2(motif_matrix + 0.0001)).sum(axis=1)
    motif_matrix_weighted = motif_matrix * ic[:, np.newaxis]
    motif_data = pd.DataFrame(motif_matrix_weighted)
    motif_data = motif_data.rename(columns={0: "A", 1: "C", 2: "G", 3: "T"})
    motif_data["pos"] = range(len(motif_data))
    motif_data = motif_data.set_index("pos")
    return motif_data


def read_motif_data(folder, motif_idx, reverse_complement, motif_file=None):
    if motif_file is None:
        motif_file = os.path.join(folder, f"motif{motif_idx}.motif")

    if not os.path.exists(motif_file):
        motif_file = os.path.join(folder, f"known{motif_idx}.motif")

    motifs = read_motifs_file(motif_file)
    assert motifs
    header, motif_matrix = motifs[0]
    info = header.split("\t")[-1].split(",")
    target_percent = float(info[0].split("(")[-1].split(")")[0].replace("%", ""))
    pval = info[-1].replace("P:", "")

    return convert_motif_matrix(motif_matrix, reverse_complement), target_percent, pval


def plot_logo(motif_data, ax_logo):
    logomaker.Logo(motif_data,
                   shade_below=0,
                   fade_below=1,
                   font_name='Arial',
                   show_spines=False,
                   baseline_width=0,
                   ax=ax_logo)

    plt.axis("off")


def make_logos_table(input_folder, output_file, motif_specs, title=None, figsize=None, adjust_top=None,
                     label_fontsize=20):
    if figsize is None:
        figsize = (8, len(motif_specs) + 1)

    motifs_to_plot = []
    for motif_spec in motif_specs:
        idx = re.match(r"^\d+", motif_spec)[0]
        reverse_complement = motif_spec[len(idx)] != "+"
        name = motif_spec[len(idx) + 1:]
        motifs_to_plot.append((idx, reverse_complement, name))

    fig = plt.figure(figsize=figsize)
    if title is not None:
        fig.suptitle(title, fontsize=15, fontfamily="Arial")

    if adjust_top is not None:
        fig.subplots_adjust(top=adjust_top)

    ncol = 4
    gs = GridSpec(len(motifs_to_plot) + 1, ncol, width_ratios=[1, 6, 1, 1])
    ax_header = fig.add_subplot(gs[0])
    ax_header.text(0.5, 0.5, "TF", horizontalalignment='center', verticalalignment='center', fontsize=label_fontsize,
                   fontfamily="Arial")

    plt.axis("off")

    ax_header = fig.add_subplot(gs[1])
    ax_header.text(0.5, 0.5, "Logo", horizontalalignment='center', verticalalignment='center', fontsize=label_fontsize,
                   fontfamily="Arial")

    plt.axis("off")

    ax_header = fig.add_subplot(gs[2])
    ax_header.text(0.5, 0.5, "% peaks", horizontalalignment='center', verticalalignment='center',
                   fontsize=label_fontsize,
                   fontfamily="Arial")

    plt.axis("off")

    ax_header = fig.add_subplot(gs[3])
    ax_header.text(0.5, 0.5, "P", horizontalalignment='center', verticalalignment='center', fontsize=label_fontsize,
                   fontfamily="Arial")

    plt.axis("off")

    for i, (motif_idx, rc, name) in enumerate(motifs_to_plot):
        ax_name = fig.add_subplot(gs[(i + 1) * ncol])
        ax_name.text(0.5, 0.5, name, horizontalalignment='center', verticalalignment='center', fontsize=label_fontsize,
                     fontfamily="Arial")

        plt.axis("off")

        ax_logo = fig.add_subplot(gs[(i + 1) * ncol + 1])
        motif_data, target_percent, pval = read_motif_data(input_folder, motif_idx, rc)
        plot_logo(motif_data, ax_logo)

        ax_info = fig.add_subplot(gs[(i + 1) * ncol + 2])
        ax_info.text(0.5, 0.5, f"{target_percent:0.1f}%", horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=15, fontfamily="Arial")

        plt.axis("off")

        ax_info = fig.add_subplot(gs[(i + 1) * ncol + 3])
        ax_info.text(0.5, 0.5, f"{pval}", horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=15, fontfamily="Arial")

        plt.axis("off")

    if adjust_top is None:
        plt.tight_layout()

    plt.savefig(output_file, dpi=600)


def homer_motif_to_meme(homer_motif_file, meme_file):
    with misc.TempFileName() as fixed_motif:
        with open(homer_motif_file) as fi:
            with open(fixed_motif, "w") as fo:
                for line in fi:
                    if line[0] == ">":
                        line = ">" + line.split("\t")[1] + "\n"
                        
                    fo.write(line)
        
        os.system(f"chen2meme {fixed_motif} >> {meme_file}")    


def plot_known_motifs_clustering(motifs_res_folder, alpha=0.05, min_fg_percent=5, min_enrichment=1.3,
                                 title="Known motif scan"):
    known_results = pd.read_csv(os.path.join(motifs_res_folder, "knownResults.txt"), sep="\t")
    motifs_db = pd.read_csv("/home/kirill/Downloads/homer/motifs/extras/motifTable.txt", sep="\t")
    motifs_db = motifs_db.set_index("Name")
    meme_file = os.path.join(motifs_res_folder, "all_known_motifs.meme")
    if os.path.exists(meme_file):
        os.remove(meme_file)

    motif_names = []
    motif_files = []
    factors = []
    qs = []
    fgs = []
    enrichments = []
    for _, row in known_results.iterrows():
        q = row["q-value (Benjamini)"]
        if q > alpha:
            continue
                    
        fg = float(row["% of Target Sequences with Motif"].replace("%", ""))
        if fg < min_fg_percent:
            continue
            
        bg = float(row["% of Background Sequences with Motif"].replace("%", ""))
        enrichment = fg / bg
        if enrichment < min_enrichment:
            continue
            
        motif_name = row["Motif Name"]
        if motif_name not in motifs_db.index:
            print(f"{motif_name} not found in motifs DB; skipping")
            continue

        motif_names.append(motif_name)
        factor = motifs_db.loc[motif_name, "Factor Name"]
        factors.append(factor)
        qs.append(q)    
        fgs.append(fg)
        enrichments.append(enrichment)
        #print(f"{motif_name}: q={q} fg={fg} enrichment={enrichment} factor: {factor}")
        motif_file = os.path.join("/home/kirill/Downloads/homer/motifs", motifs_db.loc[motif_name, "Filename"])
        motif_files.append(motif_file)
        assert os.path.exists(motif_file)
        homer_motif_to_meme(motif_file, meme_file)
        
    with misc.TempFileName() as dist_file:
        os.system(f"tomtom {meme_file} {meme_file} -text -thresh 1 > {dist_file}")
        dist_data = pd.read_csv(dist_file, sep="\t", comment="#")    
        
    dists_matrix = np.zeros((len(motif_names), len(motif_names)))
    name_idx = dict(zip(motif_names, range(len(motif_names))))
    #print(motif_names)
    for _, row in dist_data.iterrows():
        i1 = name_idx[row["Query_ID"]]
        i2 = name_idx[row["Target_ID"]]
        d = row["E-value"]
        dists_matrix[i1, i2] = d

    fig = plt.figure(figsize=(10, (len(factors) + 2) * 0.7))
    dendro_width = 0.07
    factor_width = 0.15
    p_width = 0.1
    fg_width = 0.1
    enrichment_width = 0.1
    logo_width = 1 - (dendro_width + factor_width + p_width + fg_width + enrichment_width)
    header_height = 1 / (len(factors) + 2)
    row_height = (1 - header_height * 2) / len(factors)
    space = row_height * 0.2

    ax = fig.add_axes([0, 1 - header_height, 1, header_height])
    plt.text(0.5, 0.5, title, verticalalignment="center", horizontalalignment="center")
    plt.axis("off")

    ax = fig.add_axes([dendro_width, 1 - header_height * 2, factor_width, header_height])
    plt.text(0.5, 0.5, "TF", verticalalignment="center", fontweight="bold", horizontalalignment="center")
    plt.axis("off")

    ax = fig.add_axes([dendro_width + factor_width, 1 - header_height * 2, logo_width, header_height])
    plt.text(0.5, 0.5, "Logo", verticalalignment="center", fontweight="bold", horizontalalignment="center")
    plt.axis("off")

    ax = fig.add_axes([dendro_width + factor_width + logo_width, 1 - header_height * 2, p_width, header_height])
    plt.text(0.5, 0.5, "FDR", verticalalignment="center", fontweight="bold", horizontalalignment="center")
    plt.axis("off")

    ax = fig.add_axes([dendro_width + factor_width + logo_width + p_width, 1 - header_height * 2, fg_width, header_height])
    plt.text(0.5, 0.5, "% peaks", verticalalignment="center", fontweight="bold", horizontalalignment="center")
    plt.axis("off")

    ax = fig.add_axes([dendro_width + factor_width + logo_width + p_width + fg_width, 1 - header_height * 2, enrichment_width, header_height])
    plt.text(0.5, 0.5, "Enr.", verticalalignment="center", fontweight="bold", horizontalalignment="center")
    plt.axis("off")


    ax = fig.add_axes([0.0, 0, dendro_width, 1 - header_height * 2])
    lm = linkage(dists_matrix, method='centroid')
    res = dendrogram(lm, orientation="left", labels=factors)
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)

    for i, fi in enumerate(res["leaves"]):
        ax = fig.add_axes([dendro_width, i * row_height + space / 2, factor_width, row_height - space])
        plt.text(0.0, 0.5, factors[fi], verticalalignment="center")
        plt.axis("off")
        
        ax = fig.add_axes([dendro_width + factor_width, i * row_height + space / 2, logo_width, row_height - space])
        motif_file_data = read_motifs_file(motif_files[fi])
        assert len(motif_file_data) == 1
        motif_matrix = convert_motif_matrix(motif_file_data[0][1], False)
        plot_logo(motif_matrix, ax)
        
        ax = fig.add_axes([dendro_width + factor_width + logo_width, i * row_height + space / 2, p_width, row_height - space])
        q = qs[fi]
        if q == 0.0:
            q = "< 1e-4"

        plt.text(0.5, 0.5, q, verticalalignment="center", horizontalalignment="center")
        plt.axis("off")
        
        ax = fig.add_axes([dendro_width + factor_width + logo_width + p_width, i * row_height + space / 2, fg_width, row_height - space])
        plt.text(0.5, 0.5, f"{fgs[fi]:.1f}%", verticalalignment="center", horizontalalignment="center")
        plt.axis("off")    
        
        ax = fig.add_axes([dendro_width + factor_width + logo_width + p_width + fg_width, i * row_height + space / 2, enrichment_width, row_height - space])
        plt.text(0.5, 0.5, f"{enrichments[fi]:.1f}", verticalalignment="center", horizontalalignment="center")
        plt.axis("off")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make a publication-ready figure out of Homer motif search results")
    parser.add_argument("input_folder", type=str, help="a homerResults folder")
    parser.add_argument("output_file", type=str, help="name of the output figure file")
    parser.add_argument("motif_spec", type=str, nargs="+", help="specifications for the motifs to plot: "
                                                                "<Motif idx><Strand><Name>. Strand can be + or -, + "
                                                                "being the original strand. Example: 1-DR4")

    parser.add_argument("--title", required=False, type=str, help="plot title")

    args = parser.parse_args()
    make_logos_table(args.input_folder, args.output_file, args.motif_spec, args.title,
                     figsize=(8, len(args.motif_spec) + 1))
