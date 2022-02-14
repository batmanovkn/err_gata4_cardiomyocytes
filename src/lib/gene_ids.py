import pandas as pd


def gene_id_to_name(gene_annotations, name_field="gene_name", allow_missing=False):
    id_to_name = {}
    with open(gene_annotations) as f:
        for line in f:
            if line[0] == "#":
                continue

            info = line.split("\t")[8]
            info_fields = {}
            for field in info.split("; "):
                items = field.strip().split(" ")
                name = items[0]
                value = " ".join(items[1:])
                info_fields[name] = value.replace('"', "")
            
            if name_field not in info_fields and allow_missing:
                id_to_name[info_fields["gene_id"]] = info_fields["gene_id"]
            else:
                id_to_name[info_fields["gene_id"]] = info_fields[name_field]

    return id_to_name


def refflat_gene_id_to_name(refflat_file):
    data = pd.read_csv(refflat_file, sep="\t", header=None)
    return dict(zip(data[1], data[0]))

