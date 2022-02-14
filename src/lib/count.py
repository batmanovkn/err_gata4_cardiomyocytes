import os
import glob
import json
import argparse

from tqdm.autonotebook import tqdm

from . import misc, parallel_run


def make_feature_counts(bam_files, gene_annotation_file, out_file, runner):
    featurecounts_binary = "/home/batmanov/software/subread-1.6.4-Linux-x86_64/bin/featureCounts"
    if not os.path.exists(featurecounts_binary):
        featurecounts_binary = "featureCounts"

    bam_files = " ".join(f"'{f}'" for f in bam_files)
    runner.add(f"{featurecounts_binary} -T {runner.threads} -t exon -g gene_id -a '{gene_annotation_file}' "
               f"-o '{out_file}' "
               f"{bam_files}")

               
def make_salmon_decoy(genome_fasta, transcriptome_fasta, gene_annotation_gtf, decoy_file, runner):
    decoy_gen_file = "/home/kirill/git_repos/by_others/SalmonTools/scripts/generateDecoyTranscriptome.sh"
    if not os.path.exists(decoy_gen_file):
        decoy_gen_file = "/home/batmanov/software/SalmonTools/scripts/generateDecoyTranscriptome.sh"
    
    assert os.path.exists(decoy_gen_file)
               
    runner.add(f"sh {decoy_gen_file} -g {genome_fasta} "
               f"-t {transcriptome_fasta} -a {gene_annotation_gtf} -j {runner.threads} -o {decoy_file}")
               
               
def make_salmon_index(genome_fasta, transcriptome_fasta, gene_annotation_gtf, runner=parallel_run.SimpleRunner()):
    index_folder = transcriptome_fasta + ".salmon_index"
    if not os.path.exists(index_folder):
        decoy_folder = transcriptome_fasta + ".salmon_decoy"
        if not os.path.exists(decoy_folder) or not os.path.exists(os.path.join(decoy_folder, 'gentrome.fa')):
            make_salmon_decoy(genome_fasta, transcriptome_fasta, gene_annotation_gtf, decoy_folder, runner)

        runner.add(f"salmon index -t {os.path.join(decoy_folder, 'gentrome.fa')} "
                   f"-d {os.path.join(decoy_folder, 'decoys.txt')} "
                   f"-i {index_folder} --keepDuplicates -p {runner.threads}")
               
    return index_folder
                   

def salmon_count(fastq_read1, salmon_index, output_folder, runner=parallel_run.SimpleRunner(), fastq_read2=None):
    if fastq_read2 is None:
        if isinstance(fastq_read1, str):
            assert os.path.exists(fastq_read1)
            input_arg = f"-r '{fastq_read1}'"
        else:
            assert all(os.path.exists(f) for f in fastq_read1)
            input_arg = "-r " + misc.join_quoted_paths(fastq_read1)
    else:
        assert os.path.exists(fastq_read1) and os.path.exists(fastq_read2)
        input_arg = f"-1 '{fastq_read1}' -2 '{fastq_read2}'"
               
    runner.add(f"salmon quant -i '{salmon_index}' -l a {input_arg} --validateMappings -o '{output_folder}' "
               f"-p {runner.threads}")


def make_salmon_counts_summary(all_counts_folder):
    with open(os.path.join(all_counts_folder, "summary.tsv"), "w") as f:
        f.write("sample\t%mapped\treads_mapped\n")
        for subfolder in glob.glob(os.path.join(all_counts_folder, "*")):
            if not os.path.isdir(subfolder):
                continue

            with open(os.path.join(subfolder, "aux_info", "meta_info.json")) as ff:
                info = json.load(ff)

            f.write(f"{os.path.basename(subfolder)}\t{info['percent_mapped']}\t{info['num_mapped']}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads in transcripts using Salmon")
    parser.add_argument("--index", required=True, help="salmon index folder")
    parser.add_argument("--raw_folder", required=True, help="folder with fastq files")
    parser.add_argument("--single", required=False, action="store_true", help="single read libraries")
    parser.add_argument("--output_folder", required=False, help="folder where to put counts")

    args = parser.parse_args()

    output_folder = args.output_folder
    if output_folder is None:
        output_folder = os.path.join(os.path.dirname(os.path.abspath(args.raw_folder)), "counts")

    misc.make_sure_folder_exists(output_folder)
    assert os.path.exists(args.raw_folder)

    fastas = set(glob.glob(os.path.join(args.raw_folder, "*.fastq.gz")))
    assert fastas
    pairs = []
    while fastas:
        a_file = fastas.pop()
        if "_R1_" in os.path.basename(a_file):
            paired = os.path.join(args.raw_folder, os.path.basename(a_file).replace("_R1_", "_R2_"))
            assert paired in fastas
            fastas.remove(paired)
        elif "_R2_" in os.path.basename(a_file):
            paired = os.path.join(args.raw_folder, os.path.basename(a_file).replace("_R2_", "_R1_"))
            assert paired in fastas
            fastas.remove(paired)
            a_file, paired = paired, a_file
        else:
            raise RuntimeError(f"Could not determine a pair for {a_file}")

        pairs.append([a_file, paired])

    misc.make_sure_folder_exists("logs")
    runner = parallel_run.LocalRunner(os.path.join("logs", "salmon_count.log"))
    for file_1, file_2 in tqdm(pairs):
        basename = os.path.basename(file_1)
        basename = basename[:basename.index("_R1_")]
        salmon_count(fastq_read1=file_1, fastq_read2=file_2, salmon_index=args.index,
                     output_folder=os.path.join(output_folder, basename),
                     runner=runner)

        runner.run()

    make_salmon_counts_summary(output_folder)
