import subprocess


def fastq_first_length(fastq_file):
    first_line = subprocess.check_output(f"gunzip < {fastq_file} | head -n 1", shell=True, encoding="UTF8")
    return int(first_line.split("=")[-1])
