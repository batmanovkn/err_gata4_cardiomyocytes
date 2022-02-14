import os
import glob

import pandas as pd

from . import parallel_run


def dump_fastq(run_accession, folder, runner=parallel_run.SimpleRunner(), gzip=True):
    fastq_gz = os.path.join(folder, run_accession + '_1.fastq.gz')
    if os.path.exists(fastq_gz):
        print(f"{fastq_gz} exists, skipping")
        return

    fastq_split_gz = [os.path.join(folder, run_accession + f"_{i}" + '.fastq.gz') for i in [1, 2]]
    if all(os.path.exists(f) for f in fastq_split_gz):
        print(f"{fastq_split_gz} exist, skipping")
        return

    print(f"Dumping GEO run {run_accession} into {folder}")
    runner.add(f"fastq-dump {run_accession} -O {folder} --split-files")
    if gzip:
        runner.add(f"pigz {folder}/{run_accession}*")


def read_runinfo_table(table_file, sep="\t"):
    data = pd.read_csv(table_file, sep=sep)
    return zip(data.Run, data.treatment, data.mouse_id, data.Instrument)
