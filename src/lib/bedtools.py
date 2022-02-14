import os
from typing import List
import shutil

import HTSeq
import pandas as pd

from . import misc, parallel_run


def slop(bed_file, size, chrom_sizes, output):
    os.system(f"bedtools slop -b {size} -i {bed_file} -g {chrom_sizes} > {output}")


def intersect(bed_a, bed_b, output, runner, write_orig_a=False, write_orig_b=False, left_outer_join=False):
    runner.add(f"bedtools intersect -a '{bed_a}' -b '{bed_b}' "
               f"{'-wa' * write_orig_a} {'-wb' * write_orig_b} {'-loj' * left_outer_join} > '{output}'")


def intersect_all(beds, output, runner):
    assert len(beds) != 0
    tmp_file = beds[0] + ".tmp_intersection.bed"
    runner.add(f"cp '{beds[0]}' '{output}'")
    for bed in beds[1:]:
        intersect(bed, output, tmp_file, runner)
        runner.add(f"rm '{output}'")
        runner.add(f"mv '{tmp_file}' '{output}'")


def dump_bed(bed_file, chroms, begs, ends, ids, strands, infos):
    with open(bed_file, "w") as f:
        for chrom, beg, end, id, strand, info in zip(chroms, begs, ends, ids, strands, infos):
            f.write(f"{chrom}\t{beg}\t{end}\t{id}\t.\t{strand}\t{info}\n")


def load_bed(bed_file):
    columns = misc.read_tsv(bed_file, skip_comments=True)
    return zip(*columns[:4], *columns[5:])


def compare_peaks(bed1, bed2, intersection_bed):
    intersect(bed1, bed2, intersection_bed, parallel_run.SimpleRunner(), write_orig_a=True, write_orig_b=True)

    peaks1 = pd.read_csv(bed1, sep="\t")
    peaks2 = pd.read_csv(bed2, sep="\t")
    intersection_peaks = pd.read_csv(intersection_bed, sep="\t")

    names1 = frozenset(peaks1.iloc[:, 3])
    names2 = frozenset(peaks2.iloc[:, 3])

    intersection_names1 = frozenset(intersection_peaks.iloc[:, 3])
    intersection_names2 = frozenset(intersection_peaks.iloc[:, 9])

    return names1, names2, intersection_names1, intersection_names2
                    

def chromosome_sizes(genome_fasta_file):
    chrom_sizes_file = genome_fasta_file + ".chrom_sizes"
    if not os.path.exists(chrom_sizes_file):
        print(f"Extracting chromosome sizes from {genome_fasta_file}")
        res = os.system(f"samtools faidx {genome_fasta_file}")
        assert res == 0
        res = os.system(f"cut -f1,2 {genome_fasta_file}.fai > {chrom_sizes_file}")
        assert res == 0
    
    return chrom_sizes_file


def sort_bed(bed, output, runner):
    runner.add(f"sort-bed '{bed}' > '{output}'")


def union(bed_files: List[str], output: str, runner):
    sorted_beds = [b + ".sorted" for b in bed_files]
    for b, s in zip(bed_files, sorted_beds):
        sort_bed(b, s, runner)
        
    runner.add(f"bedops --merge --ec {misc.join_quoted_paths(sorted_beds)} > '{output}'")
    misc.delete_files(sorted_beds, runner)


def read_bed_value(bed, value_column, return_sum=False):
    ga = HTSeq.GenomicArray(chroms="auto", stranded=False)
    s = 0
    with open(bed) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.replace("\n", "").split("\t")
            v = float(fields[value_column])
            ga[HTSeq.GenomicInterval(fields[0], int(fields[1]), int(fields[2]),)] = v
            s += v

    if return_sum:
        return ga, s

    return ga


def add_genomic_values(ga1, ga2):
    ga_sum = HTSeq.GenomicArray(chroms="auto", stranded=False)
    for iv, v in ga1.steps():
        ga_sum[iv] = v

    for iv, v in ga2.steps():
        ga_sum[iv] += v

    return ga_sum


def scale_genomic_values(ga, scale):
    ga_scaled = HTSeq.GenomicArray(chroms="auto", stranded=False)
    for iv, v in ga.steps():
        ga_scaled[iv] = v * scale

    return ga_scaled


def average_in(ga, iv):
    s = 0
    for siv, v in ga[iv].steps():
        s += siv.length * v

    return s / iv.length


def liftover(bed_in, transform, bed_out, runner=parallel_run.SimpleRunner()):
    assert os.path.exists(bed_in)
    chain_file = f"reference/liftOver/{transform}.over.chain.gz"
    assert os.path.exists(chain_file)
    runner.add(f"liftOver -bedPlus=5 {misc.join_quoted_paths(bed_in)} {chain_file} {misc.join_quoted_paths(bed_out)} /dev/null")
