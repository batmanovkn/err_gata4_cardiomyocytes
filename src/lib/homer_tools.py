import os
import glob
import re
from typing import Union, List
from collections import Counter
import random
import shutil

import pandas as pd
import numpy as np
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import openpyxl
import scipy.stats

from . import misc, parallel_run, bedtools, gprofiler


def make_tag_directory(bam_file, runner=None, genome_file=None, tag_directory=None, single=True):
    if runner is None:
        runner = parallel_run.SimpleRunner()

    if tag_directory is None:
        a_file = bam_file[0] if isinstance(bam_file, list) else bam_file
        tag_directory = misc.replace_extension(a_file, "tags")

    if isinstance(bam_file, list):
        assert all(os.path.exists(bf) for bf in bam_file)
        bam_file = " ".join(f"'{f}'" for f in bam_file)
    else:
        assert os.path.exists(bam_file)
        bam_file = f"'{bam_file}'"

    assert not os.path.exists(tag_directory)
    genome_option = '' if genome_file is None else f"-genome '{genome_file}' -checkGC"
    runner.add(f"makeTagDirectory '{tag_directory}' "
               f"{genome_option} {'-single' * single} "
               f"{bam_file}")

    return tag_directory


def pool_tag_directories(source_dirs: list, output, runner=parallel_run.SimpleRunner()):
    assert len(source_dirs) > 0, "No source directories given"
    assert all(os.path.exists(d) for d in source_dirs)
    inputs = " ".join(f"'{f}'" for f in source_dirs)
    runner.add(f"makeTagDirectory '{output}' -d {inputs} -keepAll -single")


def peak_file_to_bed(peak_file, runner=parallel_run.SimpleRunner(), output_bed=None):
    if output_bed is None:
        output_bed = misc.replace_extension(peak_file, "bed")

    runner.add(f"pos2bed.pl '{peak_file}' -o '{output_bed}'")
    return output_bed


def bed_file_to_peak(bed_file, runner=parallel_run.SimpleRunner(), output_peaks=None):
    if output_peaks is None:
        output_peaks = misc.replace_extension(bed_file, "peaks")

    runner.add(f"bed2pos.pl '{bed_file}' -o '{output_peaks}'")
    return output_peaks


def find_peaks(tags_dir, output_peaks, runner=parallel_run.SimpleRunner(), input_tags_dir=None, min_tags_per_peak=None,
               peak_size=None,
               fold_over_input=None,
               tags_to_define_peak=None, style=None, region=False, min_dist_between_peaks=None, p_value_over_input=None,
               fold_over_local=None, fdr_over_local=None):
    assert os.path.exists(tags_dir)
    size_option = f"-size {peak_size}" * (peak_size is not None)
    input_option = f"-i '{input_tags_dir}'" * (input_tags_dir is not None)
    style_option = f"-style {style}" * (style is not None)
    region_option = "-region" * region
    mindist_option = f"-minDist {min_dist_between_peaks}" * (min_dist_between_peaks is not None)
    p_value_over_input_option = f"-P {p_value_over_input}" * (p_value_over_input is not None)
    fold_over_local_option = f"-L {fold_over_local}" * (fold_over_local is not None)
    fdr_over_local_option = f"-LP {fdr_over_local}" * (fdr_over_local is not None)
    runner.add(f"findPeaks '{tags_dir}' {size_option} "
               f"{input_option} "
               f"{style_option} {region_option} {mindist_option} -o '{output_peaks}' "
               f"{f'-minTagThreshold {min_tags_per_peak}' * (min_tags_per_peak is not None)} "
               f"{f'-F {fold_over_input}' * (fold_over_input is not None)} "
               f"{f'-tagThreshold {tags_to_define_peak}' * (tags_to_define_peak is not None)} {p_value_over_input_option} {fold_over_local_option} "
               f"{fdr_over_local_option}")


def merge_peaks(peak_files_to_merge: List[str], output_peak_file: str, max_distance: Union[int, str],
                runner=parallel_run.SimpleRunner(),
                venn=False, genome_size=None):
    venn_option = f"-venn '{misc.replace_extension(output_peak_file, 'venn')}'" * venn
    matrix_option = f"-matrix {misc.join_quoted_paths(misc.replace_extension(output_peak_file, 'matrix'))}" + \
                    (f" -gsize {genome_size}" * (genome_size is not None))

    runner.add(f"mergePeaks -d {max_distance} {misc.join_quoted_paths(peak_files_to_merge)} {matrix_option} "
               f"{venn_option} > "
               f"'{output_peak_file}'")


def get_reproducible_peaks(replicate_peaks: List[str], output_peak_file: str, min_replicates: int,
                           max_distance: Union[int, str] = "given"):
    merged_peaks_file = output_peak_file + ".tmp"
    merge_peaks(replicate_peaks, merged_peaks_file, max_distance)
    peaks_data = pd.read_csv(merged_peaks_file, sep="\t")
    peaks_data = peaks_data[peaks_data["Total subpeaks"] >= min_replicates]
    peaks_data.to_csv(output_peak_file, sep="\t", index=False)
    os.remove(merged_peaks_file)


def compare_peaks(peak_file_1: str, peak_file_2: str, name1: str, name2: str, offsets=None, max_distance="given",
                  genome_size=None):
    assert os.path.exists(peak_file_1) and os.path.exists(peak_file_2)

    merged_peak_file = peak_file_1 + ".merged.peaks"
    merge_peaks([peak_file_1, peak_file_2], output_peak_file=merged_peak_file, max_distance=max_distance,
                runner=parallel_run.SimpleRunner(),
                venn=True, genome_size=genome_size)

    venn_file = misc.replace_extension(merged_peak_file, "venn")
    venn_data = pd.read_csv(venn_file, sep="\t")
    count_1 = venn_data.Total[venn_data.Name == peak_file_1].values[0]
    count_2 = venn_data.Total[venn_data.Name == peak_file_2].values[0]
    idx_both = venn_data.Name == f"{peak_file_1}|{peak_file_2}"
    if idx_both.sum() > 0:
        count_both = venn_data.Total[idx_both].values[0]
    else:
        count_both = 0

    v = venn2((count_1, count_2, count_both), set_labels=[name1, name2])
    if offsets:
        for label, offset in offsets.items():
            l = v.get_label_by_id(label)
            l.set_transform(l.get_transform() + mpl.transforms.Affine2D().translate(*offset))

    pval_matrix_file = misc.replace_extension(merged_peak_file, "matrix.logPvalue.matrix.txt")
    pval_data = pd.read_csv(pval_matrix_file, sep="\t")
    p_value = np.exp(float(pval_data.iloc[0, 2]))
    if p_value < 1e-100:
        pval_text = "< 1e-100"
    elif p_value >= 1:
        pval_text = "non-sig."
    else:
        pval_text = f"= {p_value:.3g}"

    plt.figtext(0, -0.15, f"co-occurrence p {pval_text}", transform=plt.gca().transAxes)
    return venn_file, (count_1, count_2, count_both), merged_peak_file, p_value


def compare_peaks_3way(peak_file_1: str, peak_file_2: str, peak_file_3: str, name1: str, name2: str, name3: str, max_distance="given"):
    assert os.path.exists(peak_file_1) and os.path.exists(peak_file_2) and os.path.exists(peak_file_3)

    merged_peak_file = peak_file_1 + ".merged.peaks"
    merge_peaks([peak_file_1, peak_file_2, peak_file_3], output_peak_file=merged_peak_file, max_distance=max_distance,
                runner=parallel_run.SimpleRunner(),
                venn=True)

    venn_file = misc.replace_extension(merged_peak_file, "venn")
    venn_data = pd.read_csv(venn_file, sep="\t")
    venn_data = venn_data.fillna(" ")
    venn_sequence = ["X  ", " X ", "XX ", "  X", "X X", " XX", "XXX"]
    counts_sequence = []
    for s in venn_sequence:
        idxs = [venn_data.iloc[:, i] == s[i] for i in range(3)]
        idx = idxs[0] & idxs[1] & idxs[2]
        if sum(idx) == 0:
            counts_sequence.append(0)
        else:
            assert sum(idx) == 1
            counts_sequence.append(venn_data.Total[idx].values[0])

    venn3(counts_sequence, [name1, name2, name3])


def intersect_peaks(peak_file_1: str, peak_file_2: str, intersection_peaks: str):
    runner = parallel_run.SimpleRunner()
    bed1 = peak_file_to_bed(peak_file_1, runner)
    bed2 = peak_file_to_bed(peak_file_2, runner)
    intersection_bed = misc.replace_extension(intersection_peaks, "bed")
    bedtools.intersect(bed1, bed2, intersection_bed, runner)
    bed_file_to_peak(intersection_bed, runner, intersection_peaks)
    os.remove(bed1)
    os.remove(bed2)
    os.remove(intersection_bed)


def intersecting_peak_ids(peaks_to_test: str, peaks_whose_presence_to_check: str, temp_file_suffix: str,
                          distance="given"):
    temp_files = [f"{temp_file_suffix}_{i}" for i in range(3)]
    runner = parallel_run.SimpleRunner()
    misc.copy_file(peaks_to_test, temp_files[0], runner)
    misc.copy_file(peaks_whose_presence_to_check, temp_files[1], runner)
    runner.add(f"mergePeaks -d {distance} {temp_files[0]} {temp_files[1]} > {temp_files[2]}")
    peaks_data = pd.read_csv(temp_files[2], sep="\t", comment=None)
    peaks_data = peaks_data.loc[~pd.isnull(peaks_data[temp_files[0]]) & ~pd.isnull(peaks_data[temp_files[1]]), :]
    peaks_id_in_intersection = set()
    for ids in peaks_data[temp_files[0]]:
        peaks_id_in_intersection.update(str(ids).split(","))

    misc.delete_files(temp_files, runner)
    return peaks_id_in_intersection


def call_peaks_vs_background_with_replicates(tag_dirs, background_tag_dirs, output_peak_file, runner,
                                             input_tag_dirs=None,
                                             style=None, genome=None, additional_options=""):
    style_option = f"-style {style}" * (style is not None)
    genome_option = f"-genome {genome}" * (genome is not None)
    if input_tag_dirs is not None:
        input_option = "-i " + misc.join_quoted_paths(input_tag_dirs)
    else:
        input_option = ""

    runner.add(f"getDifferentialPeaksReplicates.pl {style_option} {genome_option} "
               f"-t {misc.join_quoted_paths(tag_dirs)} "
               f"-b {misc.join_quoted_paths(background_tag_dirs)} {input_option} {additional_options} > "
               f"{output_peak_file}")


def get_differential_peaks(peaks, foreground_tags, background_tags, output_peak_file, runner, min_fold_change=None):
    min_fold_change_option = f"-F {min_fold_change}" * (min_fold_change is not None)
    runner.add(f"getDifferentialPeaks '{peaks}' '{foreground_tags}' '{background_tags}' {min_fold_change_option} > "
               f"'{output_peak_file}'")


def annotate_peaks(input_peaks, annotated_file, genome, runner=parallel_run.SimpleRunner(), motifs_file=None,
                   center_on_motif=None, size=None, short_chrom_names=False, custom_options=""):
    assert os.path.exists(input_peaks)

    if short_chrom_names:
        peaks_lc = misc.replace_extension(input_peaks, "lchrom.peaks")
        to_long_chromnames(input_peaks, peaks_lc)
        input_peaks = peaks_lc

    center_option = f"-center '{center_on_motif}'" * (center_on_motif is not None)
    if motifs_file is None:
        motifs_option = ""
    else:
        motifs_option = f"-m {misc.join_quoted_paths(motifs_file)}"

    size_option = f"-size {size}" * (size is not None)
    runner.add(f"annotatePeaks.pl '{input_peaks}' {genome} {motifs_option} {center_option} {size_option} {custom_options} "
               f"> '{annotated_file}'")


def annotate_tss_with_tags(genome, size, bin_size, is_heatmap, tag_dirs, output_file, gtf_file=None, strand=None):
    os.system(f"annotatePeaks.pl tss {genome} -size {size} -hist {bin_size} {'-ghist' * is_heatmap} "
              f"-d {' '.join(tag_dirs)} {'' if gtf_file is None else f'-gtf {gtf_file}'}"
              f" {f'-strand {strand}' * (strand is not None)} > {output_file}")


def make_tss_peaks(genome, size, output_file, runner=parallel_run.SimpleRunner(), gtf_file=None):
    gtf_option = "" if gtf_file is None else f"-gtf '{gtf_file}'"
    runner.add(f"annotatePeaks.pl tss {genome} -size {size} {gtf_option} > '{output_file}'")


def find_motifs(peaks, genome, sequence_size, output_folder, runner=parallel_run.SimpleRunner(), de_novo=True,
                number_of_motifs=None, lengths=None,
                background_peaks=None, redundant=None):
    lengths = "" if lengths is None else f"-len {','.join(str(l) for l in lengths)}"
    redundant_option = f"-redundant {redundant}" * (redundant is not None)
    if background_peaks is not None:
        bg_option = f'-h -bg {misc.join_quoted_paths(background_peaks)}'
    else:
        bg_option = ""

    runner.add(f"findMotifsGenome.pl {misc.join_quoted_paths(peaks)} {genome} {misc.join_quoted_paths(output_folder)}"
               f" -size {sequence_size} "
               f"{'-nomotif' * (not de_novo)} "
               f"{f'-S {number_of_motifs}' * (number_of_motifs is not None)} "
               f"{lengths} -p {runner.threads} {bg_option} "
               f"{redundant_option}")

    zip_file = output_folder + ".zip"
    runner.add(f"zip -r {misc.join_quoted_paths(zip_file)} {misc.join_quoted_paths(output_folder)}")


def find_motifs_in_promoters(genes_file, genome, output_folder, runner, de_novo=True,
                             number_of_motifs=None,
                             lengths=None, background=None, threads=15):
    lengths = "" if lengths is None else f"-len {','.join(str(l) for l in lengths)}"
    runner.add(f"findMotifs.pl {genes_file} {genome} {output_folder} "
               f"{'-nomotif' * (not de_novo)} "
               f"{f'-S {number_of_motifs}' * (number_of_motifs is not None)} "
               f"{lengths} -p {threads} {f'-bg {background}' * (background is not None)}")


def count_tags(peaks, tags_dirs, genome, sequence_size, output, histogram=False, runner=parallel_run.SimpleRunner(),
               bin_size=None, norm=True):
    assert isinstance(tags_dirs, list) or isinstance(tags_dirs, str)
    hist_options = f"-hist {bin_size} -ghist" * histogram
    tags_list = " ".join(f"'{d}'" for d in tags_dirs) if isinstance(tags_dirs, list) else f"'{tags_dirs}'"
    norm_option = None
    if norm == "rlog":
        norm_option = "-rlog"
    elif norm == "none":
        norm_option = "-raw"
    elif norm == True:
        norm_option = ""

    assert norm_option is not None
    runner.add(
        f"annotatePeaks.pl '{peaks}' {genome} -size {sequence_size} {hist_options} -cpu {runner.threads} "
        f"{norm_option} "
        f"-d {tags_list} > "
        f"'{output}'")


def motif_histogram_in_regions(peaks, motif_files, genome, sequence_size, bin_size, output):
    os.system("annotatePeaks.pl %s %s -m %s -size %d -hist %d > %s" %
              (peaks, genome, " ".join(motif_files), sequence_size, bin_size, output))


def grep_known_motifs(pattern, motifs_result_folder, best_only, max_p=1):
    matching_files = []
    for motif_file in glob.glob(os.path.join(motifs_result_folder, "knownResults", "known*.motif")):
        with open(motif_file) as f:
            line = f.readline()
            name = line.split("\t")[1]
            p = line.split("\t")[-1].split(",")[-1]
            assert p.startswith("P:")
            p = float(p.replace("P:", "").replace("\n", ""))
            if p < max_p and re.match(pattern, name, re.IGNORECASE):
                matching_files.append((motif_file, p, name))

    if best_only:
        if not matching_files:
            return None

        best_match = min(matching_files, key=lambda x: x[1])
        return best_match[0], best_match[1]

    return [(f, n) for f, p, n in matching_files]


def make_bigwig(tags_dir, chrom_sizes, runner=None, style=None, sort=False, remake=False):
    if runner is None:
        runner = parallel_run.SimpleRunner()

    bigwig_file = os.path.join(tags_dir, os.path.basename(tags_dir) + ".ucsc.bigWig")
    if remake and os.path.exists(bigwig_file):
        os.remove(bigwig_file)

    if not os.path.exists(bigwig_file):
        assert os.path.exists(chrom_sizes)
        runner.add(f"makeUCSCfile '{tags_dir}' -o auto -bigWig '{chrom_sizes}' "
                   f"{f'-style {style}' * (style is not None)}")

        if sort:
            runner.add(f"LC_COLLATE=C sort -k1,1 -k2,2n -o '{bigwig_file}' '{bigwig_file}'")
            tmp_file = bigwig_file + ".tmp"
            runner.add(f"bedGraphToBigWig '{bigwig_file}' '{chrom_sizes}' '{tmp_file}'")
            misc.delete_files(bigwig_file, runner)
            misc.rename_file(tmp_file, bigwig_file, runner)

    return bigwig_file


def extract_fasta(peaks, genome_fasta_file, output_fasta, runner=parallel_run.SimpleRunner(), tab_delimited=False):
    runner.add(f"homerTools extract '{peaks}' '{genome_fasta_file}' {'-fa' * (not tab_delimited)} > '{output_fasta}'")


def to_long_chromnames(peaks: str, long_chromnames_file=None):
    if long_chromnames_file is None:
        if peaks.endswith(".peaks"):
            long_chromnames_file = peaks.replace(".peaks", ".lchrom.peaks")
        else:
            long_chromnames_file = peaks + ".lchrom.peaks"

    peaks_data = pd.read_csv(peaks, sep="\t", comment="#", header=None, dtype=str)
    peaks_data.iloc[:, 1] = [f"chr{chrom}" for chrom in peaks_data.iloc[:, 1]]
    peaks_data.to_csv(long_chromnames_file, sep="\t", index=False, header=False)
    return long_chromnames_file


def to_short_chromnames(peaks, short_chromnames_file=None, chrom_column=1):
    if short_chromnames_file is None:
        if peaks.endswith(".peaks"):
            short_chromnames_file = peaks.replace(".peaks", ".schrom.peaks")
        else:
            short_chromnames_file = peaks + ".schrom.peaks"

    peaks_data = pd.read_csv(peaks, sep="\t", comment="#", header=None, dtype=str)
    peaks_data.iloc[:, chrom_column] = [chrom.replace("chr", "") for chrom in peaks_data.iloc[:, chrom_column]]
    peaks_data.to_csv(short_chromnames_file, sep="\t", index=False, header=False)
    return short_chromnames_file


def plot_area_kind_distribution(peaks_kind, kind_name, frequency=False):
    ko_regions = []
    ko_dir = []
    for d, annotated_peak_file in peaks_kind:
        annotated_data = pd.read_csv(annotated_peak_file, sep="\t")
        region_type = [t.split(" ")[0] for t in annotated_data.Annotation if not pd.isna(t)]
        ko_regions.extend(region_type)
        ko_dir.extend([d] * len(region_type))

    df = pd.DataFrame({kind_name: ko_dir, "region": ko_regions})
    counts_df = (df["region"]
               .groupby(df[kind_name])
               .value_counts(normalize=True)
               .rename("prop")
               .reset_index())

    region_order = ["Intergenic", "intron", "promoter-TSS", "exon", "non-coding", "TTS", "5'", "3'"]
    if len(peaks_kind) == 1:
        sns.countplot(x="region", data=df, order=region_order)
        plt.xticks(rotation="vertical")
        plt.xlabel("genome region")
    else:
        if not frequency:
            sns.countplot(x=kind_name, hue="region", data=df,
                          hue_order=region_order)

            plt.ylabel("# peaks")
        else:
            counts_df.prop *= 100
            sns.barplot(x=kind_name, y="prop", hue="region", data=counts_df,
                        hue_order=region_order,
                        order=[d for d, _ in peaks_kind])

            plt.ylabel("% peaks")

    if frequency:
        plt.ylabel("% peaks")
    else:
        plt.ylabel("# peaks")

    counts_df = (df["region"]
               .groupby(df[kind_name])
               .value_counts(normalize=False)
               .rename("prop")
               .reset_index())

    return counts_df


def call_updown_changes(condition1_tag_dirs, condition2_tag_dirs, base_peaks_file_name, runner,
                        additional_peak_calling_options=""):
    peak_up_file = base_peaks_file_name + "_up.peaks"
    call_peaks_vs_background_with_replicates(condition1_tag_dirs, condition2_tag_dirs, peak_up_file, runner,
                                             additional_options=additional_peak_calling_options)
    peak_down_file = base_peaks_file_name + "_down.peaks"
    call_peaks_vs_background_with_replicates(condition2_tag_dirs, condition1_tag_dirs, peak_down_file, runner,
                                             additional_options=additional_peak_calling_options)

    return peak_up_file, peak_down_file


def call_updown_changes_no_replicates(peaks_condition1, peaks_condition2,
                                      condition1_tag_dirs, condition2_tag_dirs, base_peaks_file_name, runner,
                                      min_fold_change=None):
    peak_up_file = base_peaks_file_name + "_down.peaks"
    get_differential_peaks(peaks_condition1, foreground_tags=condition1_tag_dirs, background_tags=condition2_tag_dirs,
                           output_peak_file=peak_up_file, runner=runner, min_fold_change=min_fold_change)

    peak_down_file = base_peaks_file_name + "_up.peaks"
    get_differential_peaks(peaks_condition2, foreground_tags=condition2_tag_dirs, background_tags=condition1_tag_dirs,
                           output_peak_file=peak_down_file, runner=runner, min_fold_change=min_fold_change)

    return peak_up_file, peak_down_file


def peak_density_heatmap(peak_file, tag_dirs: dict, runner=parallel_run.SimpleRunner(), recount=True, log_scale=False,
                         normalize_separately=False, tick_rotation="horizontal", sort_by=None):
    counts_file = peak_file + ".counts"
    if recount or not os.path.exists(counts_file):
        count_tags(peak_file, list(tag_dirs.values()), "none", 2000, counts_file, runner=runner, histogram=True,
                   bin_size=10)

    counts_data = pd.read_csv(counts_file, sep="\t")

    plotting_data = []
    counts = counts_data.iloc[:, 1:].values
    assert counts.shape[1] % len(tag_dirs) == 0
    cols = counts.shape[1] // len(tag_dirs)
    idx = np.arange(counts.shape[0])
    for ds_name, s in zip(tag_dirs.keys(), range(0, counts.shape[1], cols)):
        dataset = counts[:, s:s + cols]
        if log_scale:
            dataset = np.log1p(dataset)

        if normalize_separately:
            dataset /= np.quantile(dataset, 0.95)

        plotting_data.append(dataset)
        if s + cols < counts.shape[1]:
            plotting_data.append(np.zeros((counts.shape[0], 10)))

        if sort_by is not None and ds_name == sort_by:
            idx = np.argsort(-dataset.sum(axis=1))

    plotting_data = np.concatenate(plotting_data, axis=1)[idx, :]
    if not normalize_separately:
        plotting_data /= np.quantile(plotting_data, 0.95)

    palette = sns.cubehelix_palette(start=0.9, rot=0, dark=0.2, light=1, n_colors=100, as_cmap=True, hue=0.8)
    plt.imshow(plotting_data, cmap=palette, aspect="auto", vmin=0, vmax=1)
    plt.xticks([(cols + 10) * i + cols // 2 for i in range(len(tag_dirs))], labels=list(tag_dirs.keys()),
               rotation=tick_rotation)

    plt.yticks([])
    ax = plt.gca()
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis=u'both', which=u'both', length=0)
    return counts.shape[0]


def peak_average_density_plot(peak_file, tag_dirs: dict, runner=parallel_run.SimpleRunner(), size=2000, bin_size=10,
                              recount=True, palette=None, **kwargs):
    counts_file = peak_file + ".counts"
    border_size = 200
    if recount or not os.path.exists(counts_file):
        count_tags(peak_file, list(tag_dirs.values()), "none", size + border_size * 2, counts_file, runner=runner,
                   histogram=True,
                   bin_size=bin_size)

    counts_data = pd.read_csv(counts_file, sep="\t")
    counts = counts_data.iloc[:, 1:].values
    assert counts.shape[1] % len(tag_dirs) == 0
    cols = counts.shape[1] // len(tag_dirs)
    border_cols = border_size // bin_size
    col_positions = counts_data.columns[1 + border_cols:1 + cols - border_cols]
    for ds_name, s in zip(tag_dirs.keys(), range(0, counts.shape[1], cols)):
        dataset = counts[:, s + border_cols:s + cols - border_cols]
        dataset = pd.DataFrame(dataset, columns=col_positions)
        dataset = pd.melt(dataset)
        dataset["variable"] = dataset["variable"].astype(int)
        if palette is not None:
            sns.lineplot(data=dataset, x="variable", y="value", label=ds_name, color=palette[ds_name], **kwargs)
        else:
            sns.lineplot(data=dataset, x="variable", y="value", label=ds_name, **kwargs)

    plt.ylabel("Reads per 10M")
    plt.xlabel("Distance from peak center")
    return counts.shape[0]


def reads_scatterplot(peak_file, tag_dir_x, tag_dir_y, name_x, name_y, size=1000, runner=parallel_run.SimpleRunner(),
                      recount=True, log_scale=False):
    counts_file = peak_file + ".counts_sum"
    if recount or not os.path.exists(counts_file):
        count_tags(peak_file, [tag_dir_x, tag_dir_y], "none", size, counts_file, runner=runner)

    counts_data = pd.read_csv(counts_file, sep="\t")
    x_column = [c for c in counts_data.columns if c.startswith(tag_dir_x + " Tag Count")]
    assert len(x_column) == 1
    x_column = x_column[0]

    y_column = [c for c in counts_data.columns if c.startswith(tag_dir_y + " Tag Count")]
    assert len(y_column) == 1
    y_column = y_column[0]
    p = sns.scatterplot(data=counts_data, x=x_column, y=y_column, s=5, linewidth=0)
    plt.xlabel(name_x)
    plt.ylabel(name_y)
    if log_scale:
        p.set(xscale="log")
        p.set(yscale="log")

    r = np.corrcoef(counts_data[x_column], counts_data[y_column])[0, 1]
    return len(counts_data), r


def annotate_and_export(peak_file, genome, organism, gene_id_map, prefix, peaks_folder, pics_folder, tables_folder=None,
                        tables_file=None,
                        chrom_names_are_short=True):
    if isinstance(gene_id_map, str):
        gene_id_map = pd.read_csv(gene_id_map, sep="\t")
        gene_id_map = dict(zip(gene_id_map.id, gene_id_map.name))

    peaks_annotated = os.path.join(peaks_folder, f"{prefix}.annotated.peaks")
    print(f"Create {peaks_annotated}")
    annotate_peaks(peak_file, peaks_annotated, genome=genome, short_chrom_names=chrom_names_are_short)
    assert os.path.exists(peaks_annotated)

    bed_file = os.path.join(peaks_folder, f"{prefix}.bed")
    print(f"Save {bed_file}")
    peak_file_to_bed(peaks_annotated, output_bed=bed_file)

    if tables_file is None:
        assert tables_folder is not None
        tables_file = os.path.join(tables_folder, f'{prefix}_peaks.xlsx')

    sheet_name = f"{prefix} peaks"
    print(f"Add [{sheet_name}] to {tables_file}")
    writer = pd.ExcelWriter(tables_file, engine='openpyxl', mode="a" if os.path.exists(tables_file) else "w")
    wb = writer.book
    if sheet_name in wb:
        del wb[sheet_name]

    peak_data = pd.read_csv(peaks_annotated, sep="\t")
    peak_data = peak_data.rename(columns={peak_data.columns[0]: "peak_id"})
    peak_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

    pic_file = os.path.join(pics_folder, f"{prefix}_peak_regions.png")
    print(f"Create {pic_file}")
    plt.figure(figsize=(6, 5))
    plot_area_kind_distribution([("", peaks_annotated)], "")
    plt.title(f"{prefix} peak regions")
    plt.subplots_adjust(bottom=0.4, left=0.2)
    plt.savefig(pic_file, dpi=300)

    print(f"Add [{prefix}|GO] to {tables_file}")
    wb = openpyxl.load_workbook(tables_file)
    genes = list(set(peak_data["Nearest Ensembl"][~pd.isna(peak_data["Nearest Ensembl"])]))
    gprofiler.go_analysis_for_genes(genes, prefix, gene_id_map, wb, 0.05, 2, os.path.join(pics_folder, "GO", prefix),
                                    organism)

    wb.save(tables_file)

    return peaks_annotated


def compare_and_export_peaks(peaks1, peaks2, genome, organism, gene_id_map, prefix1, prefix2, prefix_both, peaks_folder,
                             pics_folder, tables_folder,
                             chrom_names_are_short=True, max_distance="given", venn_offsets=None, genome_size=None):
    tables_file = os.path.join(tables_folder, f"{prefix_both}_peaks.xlsx")
    if os.path.exists(tables_file):
        os.remove(tables_file)

    misc.make_sure_folder_exists(pics_folder)
    misc.make_sure_folder_exists(peaks_folder)
    plt.figure(figsize=(6, 6))
    _, _, merged_peaks_file, _ = compare_peaks(peaks1, peaks2, prefix1, prefix2, max_distance=max_distance,
                                               offsets=venn_offsets, genome_size=genome_size)

    plt.savefig(os.path.join(pics_folder, f"{prefix_both}_venn.png"), dpi=300)
    peaks_data = pd.read_csv(merged_peaks_file, sep="\t")
    peaks1_present = ~pd.isnull(peaks_data[peaks1])
    peaks2_present = ~pd.isnull(peaks_data[peaks2])
    peaks1_only = peaks1_present & ~peaks2_present
    peaks2_only = peaks2_present & ~peaks1_present
    both_present = peaks1_present & peaks2_present
    for idx, name in [(peaks1_only, f"{prefix1}_only"), (peaks2_only, f"{prefix2}_only"),
                      (both_present, f"{prefix_both}")]:
        print(f"Exporting {name}")
        selected_data = peaks_data[idx]
        selected_peaks_file = os.path.join(peaks_folder, f"{name}.peaks")
        selected_data.to_csv(selected_peaks_file, sep="\t", index=False, header=False)
        annotate_and_export(selected_peaks_file, genome, organism, gene_id_map, name, peaks_folder, pics_folder,
                            tables_file=tables_file,
                            chrom_names_are_short=chrom_names_are_short)


def pseudoreplicate(tag_dir, out1, out2):
    tmp1 = out1 + ".tmp"
    tmp2 = out2 + ".tmp"
    os.mkdir(tmp1)
    os.mkdir(tmp2)
    with open(os.path.join(tag_dir, "genome.tags.tsv"), "r") as in_f:
        with open(os.path.join(tmp1, "genome.tags.tsv"), "w") as out1_f:
            with open(os.path.join(tmp2, "genome.tags.tsv"), "w") as out2_f:
                for line in in_f:
                    fields = line.split("\t")
                    num_of_reads = round(float(fields[4]) * 2)
                    assert num_of_reads > 0
                    reads_to_1 = random.randint(0, num_of_reads)
                    reads_to_2 = num_of_reads - reads_to_1
                    if reads_to_1 > 0:
                        out1_f.write("\t".join(fields[:4] + [str(reads_to_1 / 2), fields[5]]))
                    if reads_to_2 > 0:
                        out2_f.write("\t".join(fields[:4] + [str(reads_to_2 / 2), fields[5]]))

    pool_tag_directories([tmp1], out1)
    pool_tag_directories([tmp2], out2)
    shutil.rmtree(tmp1)
    shutil.rmtree(tmp2)


def compare_gene_sets(genes1, genes2, name1, name2, background):
    if not isinstance(genes1, set):
        genes1 = set(genes1)

    if not isinstance(genes2, set):
        genes2 = set(genes2)

    if not isinstance(background, set):
        background = set(background)

    genes1 = genes1 & background
    genes2 = genes2 & background

    venn2([genes1, genes2], [name1, name2])
    _, p_value = scipy.stats.fisher_exact([[len(background - genes1 - genes2), len(genes1)],
                                           [len(genes2), len(genes1 & genes2)]])

    if p_value < 1e-100:
        pval_text = "< 1e-100"
    elif p_value >= 1:
        pval_text = "non-sig."
    else:
        pval_text = f"= {p_value:.3g}"

    plt.figtext(0, -0.15, f"association p {pval_text}", transform=plt.gca().transAxes)
