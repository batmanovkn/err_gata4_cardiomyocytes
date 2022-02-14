import os
import subprocess

from tqdm.autonotebook import tqdm

from . import misc
from . import parallel_run


def star_executable():
    star_path = "/home/batmanov/software/STAR-2.7.1a/source/STAR"
    if not os.path.exists(star_path):
        star_path = "STAR"

    return star_path


def star_index_file(genome_file, gene_annotations, runner=None):
    index_folder = genome_file + ".STAR_index"
    if not os.path.exists(index_folder):
        if runner is None:
            raise RuntimeError(f"STAR index is missing for {genome_file}")
            
        print(f"Making STAR index for {genome_file}")
        runner.add(f"mkdir {index_folder}")
        runner.add(f"{star_executable()} --runThreadN {runner.threads} --runMode genomeGenerate "
                   f"--genomeDir {index_folder} "
                   f"--genomeFastaFiles {genome_file} --sjdbGTFfile {gene_annotations}"
                   f" --sjdbOverhang 100"
                   #f" --limitGenomeGenerateRAM 35000000000"
                   )

    return index_folder


def bowtie_index_file(genome_file, make_index, runner):
    index_folder = genome_file + ".Bowtie_index"
    prefix = os.path.join(index_folder, "idx")
    if not os.path.exists(index_folder):
        if not make_index:
            raise RuntimeError(f"bowtie2 index is missing for {genome_file}")
            
        print(f"Making Bowtie index for {genome_file}")
        runner.add(f"mkdir '{index_folder}'")
        runner.add(f"bowtie2-build --threads {runner.threads} -f '{genome_file}' '{prefix}'")

    return prefix


def index_bam(bam_file, runner):
    runner.add(f"samtools index '{bam_file}'")


def num_reads_in_bam(bam_file, properly_paired=True):
    paired_flag = "-f 3" * properly_paired
    return int(parallel_run.capture_stdout(f"samtools view -c {paired_flag} {misc.join_quoted_paths(bam_file)}"))


def align_rnaseq_star(fastq_file_read1, fastq_file_read2, out_dir, genome_file=None, gene_annotations=None, runner=None,
                      star_index_folder=None, multimap_max=1, res_file=None):
    assert runner is not None

    if out_dir is None:
        out_dir = os.path.dirname(fastq_file_read1)

    misc.make_sure_folder_exists(out_dir)

    prefix = os.path.join(out_dir, misc.remove_extension(os.path.basename(fastq_file_read1).replace(".gz", "")))
    aligned_bam = prefix + "Aligned.sortedByCoord.out.bam"

    if not os.path.exists(aligned_bam):
        files = "'" + fastq_file_read1 + "'"
        if fastq_file_read2 is not None:
            files += " '" + fastq_file_read2 + "'"
        
        runner.add(f"echo Aligning {files}")
        is_local = isinstance(runner, parallel_run.LocalRunner)
        read_files_option = "--readFilesCommand zcat" * (fastq_file_read1.endswith(".gz"))

        if star_index_folder is None:
            star_index_folder = star_index_file(genome_file, gene_annotations, is_local)
            
        runner.add(f"{star_executable()} --runMode alignReads --runThreadN {runner.threads} --genomeDir "
                   f"{star_index_folder} "
                   f"--readFilesIn {files} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate "
                   f"--outFilterMultimapNmax {multimap_max} --outFileNamePrefix {prefix} {read_files_option}")

        misc.delete_files([prefix + temp_suffix for temp_suffix in ["Log.out", "Log.progress.out", "SJ.out.tab"]],
                          runner)

    if res_file is not None:
        misc.rename_file(aligned_bam, res_file, runner)
        aligned_bam = res_file

    if not os.path.exists(aligned_bam + ".bai"):
        index_bam(aligned_bam, runner)

    return aligned_bam, aligned_bam.replace("Aligned.sortedByCoord.out.bam", "Log.final.out")


def align_rnaseq_experiment(fastq_files, rename_function, genome_file, gene_annotations, local_run, dry_run=False,
                            aligned_folder=None):
    summary = []
    first_file = fastq_files[0]
    if not isinstance(first_file, str):
        first_file = first_file[0]

    fastq_folder = os.path.dirname(first_file)
    if aligned_folder is None:
        aligned_folder = os.path.join(fastq_folder, "aligned")
        
    print("Aligned files folder:", aligned_folder)
    for i, file1 in enumerate(tqdm(fastq_files)):
        if isinstance(file1, str):
            file2 = None
        else:
            file1, file2 = file1
        
        log_file = misc.replace_extension(file1, "alignment_log")
        if local_run:
            runner = parallel_run.LocalRunner(log_file)
        else:
            runner = parallel_run.PMACSRunner(log_file, threads=16, memory_gb=100, name=f"STAR_{i + 1}_{len(fastq_files)}")
            
            runner.add("source /etc/profile.d/modules.sh")
            # runner.add("module add samtools-1.1")
            # runner.add("export PATH=~/software/STAR-2.7.1a/source:$PATH")

        bam_file, bam_log_file = align_rnaseq_star(file1, file2, aligned_folder, genome_file, gene_annotations, runner)
        good_name = rename_function(bam_file)
        if good_name is None:
            print(f"Skipping {file1} as requested by renaming function")
            continue
        
        summary.append((os.path.basename(good_name), os.path.basename(bam_file)))

        if os.path.exists(good_name):
            print(f"Skipping {file1} because destination {good_name} exists")
            continue
        
        misc.rename_file(bam_file, good_name, runner)
        misc.rename_file(bam_file + ".bai", good_name + ".bai", runner)
        good_log_name = good_name + ".log"
        misc.rename_file(bam_log_file, good_log_name, runner)
        runner.run(dry=dry_run)

    with open(os.path.join(aligned_folder, "summary.tsv"), "w") as f:
        f.write("\t".join(["file name", "original sample name"]) + "\n")
        for file, orig_name in sorted(summary):
            f.write(file + "\t" + orig_name + "\n")
            
            
def make_summary(aligned_folder):
    grab_fields = ["Number of input reads", "Average input read length", "Uniquely mapped reads %",
                   "Mismatch rate per base, %", "% of reads mapped to too many loci", "% of reads unmapped: too short",
                   "% of reads unmapped: other"]
    
    summary = []
    summary_file = os.path.join(aligned_folder, "summary.tsv")
    with open(summary_file) as f:
        f.readline()
        for line in f:
            file, orig_name = line.replace("\n", "").split("\t")
            log_file = os.path.join(aligned_folder, file + ".log")
            
            summary_for_file = {}
            with open(log_file) as ff:
                for line in ff:
                    for field in grab_fields:
                        if field in line:
                            summary_for_file[field] = line.split("|")[1].strip()
                            break

            summary.append((file, orig_name, summary_for_file))
            
    with open(summary_file, "w") as f:
        f.write("\t".join(["file name", "original sample name"] + grab_fields) + "\n")
        for file, orig_name, s in sorted(summary):
            f.write(file + "\t" + orig_name + "\t" + "\t".join(s[field] for field in grab_fields) + "\n")


def align_chipseq(fastq, out_dir, genome_file, runner, max_mem_gb=25, unique_only=True):
    file1 = fastq
    if isinstance(file1, list):
        file1 = file1[0]

    if out_dir is None:
        out_dir = os.path.dirname(file1)

    misc.make_sure_folder_exists(out_dir)

    prefix = os.path.join(out_dir, misc.remove_extension(os.path.basename(file1)))
    unique1_file = prefix + ".unique1.bam"
    if not os.path.exists(unique1_file):
        bam_file = prefix + ".bam"

        if isinstance(fastq, str):
            input_spec = f"-U '{fastq}'"
            assert os.path.exists(fastq)
        else:
            input_spec = f"-1 '{fastq[0]}' -2 '{fastq[1]}'"
            assert all(os.path.exists(fq) for fq in fastq)

        if not os.path.exists(bam_file):
            runner.add(f"bowtie2 -p {runner.threads} -N 1 -t -x {bowtie_index_file(genome_file, False, runner)} "
                       f"{input_spec} | samtools view -bSF4 - > '{bam_file}'")

        sorted_bam = prefix + ".sorted.bam"
        mem = max_mem_gb // runner.threads
        if mem == 0:
            mem = 1

        runner.add(f"samtools sort -@ {max_mem_gb // mem} -o '{sorted_bam}' -O bam -T '{sorted_bam}' "
                   f"-m {mem}G '{bam_file}'")

        if unique_only:
            dedup_file = prefix + ".unique.bam"
            runner.add(f"samtools rmdup -s '{sorted_bam}' '{dedup_file}'")
        else:
            dedup_file = sorted_bam

        runner.add(f"samtools view -bq 1 '{dedup_file}' > '{unique1_file}'")

        index_bam(unique1_file, runner)
        misc.delete_files(list({bam_file, sorted_bam, dedup_file}), runner)

    return unique1_file


def align_chipseq_experiment(fastq_files, rename_function, genome_file, local_run, dry_run, aligned_folder=None):
    first_file = fastq_files[0]
    if isinstance(first_file, list):
        first_file = first_file[0]

    fastq_folder = os.path.dirname(first_file)
    if aligned_folder is None:
        aligned_folder = os.path.join(fastq_folder, "aligned")

    print("Aligned files folder:", aligned_folder)
    logs_folder = os.path.join(aligned_folder, "logs")
    misc.make_sure_folder_exists(logs_folder)

    summary = []
    for i, fastq in enumerate(fastq_files):
        file1 = fastq
        if isinstance(file1, list):
            file1 = file1[0]

        log_file = os.path.join(logs_folder, misc.replace_extension(os.path.basename(file1), "alignment_log"))
        if local_run:
            runner = parallel_run.LocalRunner(log_file)
        else:
            runner = parallel_run.PMACSRunner(log_file, threads=16, memory_gb=100,
                                              name=f"bowtie_{i + 1}_{len(fastq_files)}")

            runner.add("source /etc/profile.d/modules.sh")
            # runner.add("module add samtools-1.1")
            runner.add("module add bowtie2/2.3.4.1")

        bam_file = align_chipseq(fastq, aligned_folder, genome_file, runner)
        good_name = rename_function(bam_file)
        if good_name is None:
            print(bam_file, "is dropped")
            continue
            
        summary.append((os.path.basename(good_name), os.path.basename(bam_file)))

        if os.path.exists(good_name):
            print(f"Skipping {file1} because destination {good_name} exists")
            continue

        misc.rename_file(bam_file, good_name, runner)
        misc.rename_file(bam_file + '.bai', good_name + '.bai', runner)

        runner.run(dry=dry_run)

    with open(os.path.join(aligned_folder, "summary.tsv"), "w") as f:
        f.write("\t".join(["file name", "original sample name"]) + "\n")
        for file, orig_name in sorted(summary):
            f.write(file + "\t" + orig_name + "\n")


def bam_pbc(bam_file):
    return float(parallel_run.capture_stdout(f"R -q -s -e \"cat(encodeChIPqc::PBC('{bam_file}'))\"",
                                             suppress_stderr=True))


def chipseq_alignment_qc(fastq_r1, fastq_r2, aligned_qc_folder, genome_fasta):
    misc.make_sure_folder_exists(aligned_qc_folder)
    base_name = os.path.basename(fastq_r1)
    log_file = os.path.join("logs", f"qc_align_{base_name}.log")
    runner = parallel_run.LocalRunner(log_file, threads=40)
    bam = align_chipseq([fastq_r1, fastq_r2], aligned_qc_folder, genome_fasta, runner, 40, unique_only=False)
    runner.run()
    pbc1 = bam_pbc(bam)
    bam_reads = num_reads_in_bam(bam)
    with open(log_file, "r") as f:
        lines = f.readlines()

    alignment_lines = [l for l in lines if "overall alignment rate" in l]
    alignment_percent = alignment_lines[-1].split(" ")[0]

    num_unique_reads = int(parallel_run.capture_stdout(f"samtools rmdup {bam} - | samtools view -c -f 3 -",
                                                       suppress_stderr=True))

    nrf = num_unique_reads / bam_reads

    return {"file": fastq_r1, "aligned_reads": bam_reads, "alignment_percent": alignment_percent, "PBC1": pbc1,
            "NRF": nrf}
