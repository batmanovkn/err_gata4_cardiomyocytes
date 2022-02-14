import argparse
import fileinput

import openpyxl
import pandas as pd

from . import gprofiler


libraries = {
    "KEGG": "KEGG pathways",
    "GO:MF": "Molecular function",
    "GO:CC": "Cellular component",
    "GO:BP": "Biological process"
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform gProfiler-based GO term enrichment analysis")
    parser.add_argument("--alpha", required=False, default=0.05, type=float, help="p-value cutoff")
    parser.add_argument("--enrichment_cutoff", required=False, default=2.0, type=float, help="minimum enrichment")
    parser.add_argument("--gene_name_map", required=False, type=str, help=".tsv file with gene names")
    parser.add_argument("--top_n", required=False, default=None, type=int, help="maximum number of results to plot in figures")    

    parser.add_argument("organism", type=str, choices=list(gprofiler.organisms), help="organism DB to use")
    parser.add_argument("tag", type=str, help="name tag for results")
    parser.add_argument("result_tables_file", type=str, help="path to .xlsx file with results")
    parser.add_argument("pics_folder", type=str, help="folder to put pictures in")

    args = parser.parse_args()

    genes = []
    for line in fileinput.input(files="-"):
        genes.append(line.strip())

    if args.gene_name_map:
        gene_ids = pd.read_csv(args.gene_name_map, sep="\t")
        gene_ids = dict(zip(gene_ids["id"], gene_ids.name))
    else:
        gene_ids = None

    res = gprofiler.plot_significant_results(genes, libraries, gprofiler.organisms[args.organism],
                                             args.enrichment_cutoff, args.alpha,
                                             args.pics_folder, None, suffix=args.tag, return_genes=True, top_n_to_plot=args.top_n)

    wb = openpyxl.load_workbook(args.result_tables_file)
    ws = wb.create_sheet(title=args.tag[:26])

    gprofiler.write_res_from_libraries(ws, res, libraries, enrichment_cutoff=args.enrichment_cutoff,
                                       gene_name_map=gene_ids)

    wb.save(args.result_tables_file)
