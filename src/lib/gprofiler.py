import os

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import requests
import pandas as pd

import matplotlib.style
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib_venn import venn3, venn2

# Say, "the default sans-serif font is Arial"
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.size'] = 15
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.loc'] = "upper left"
mpl.rcParams['figure.facecolor'] = "white"

from . import misc, kegg_info


organisms = {"mouse": 'mmusculus', "human": "hsapiens", "dog": "cfamiliaris", "rat": "rnorvegicus", "pig": "sscrofa"}

main_libraries = {
    "KEGG": "KEGG pathways",
    "GO:MF": "Molecular function",
    "GO:CC": "Cellular component",
    "GO:BP": "Biological process"
}


def query_gost(genes, library, cutoff, organism, return_genes=False, electronic_evidence=False,
               ordered=False, background=None):
    params = {
        'organism': organism,
        'query': "\n".join(genes),
        'sources': [library],
        'user_threshold': cutoff,
        'no_evidences': not return_genes,   # skip lookup for evidence codes. Speeds up queries, if there is no
                                            # interest in
                                            # evidence codes.
        'no_iea': not electronic_evidence,  # Ignore electronically annotated GO annotations
        "ordered": ordered                  # ordered query: the input genes are assumed ordered by importance
    }
       
    if background is not None:
        params['domain_scope'] = 'custom'  # use the genes in the probe as the statistical background.
        params['background'] = "\n".join(background)
            
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json=params,
        headers={
            'User-Agent': 'FullPythonRequest'
        }
    )

    res = r.json()
    #print(res)
    if "result" not in res:
        raise RuntimeError("gProfiler returned an error: " + str(res))
        
    if return_genes:
        accepted_genes = np.array(list(res["meta"]["genes_metadata"]["query"]["query_1"]["mapping"].keys()), str)
    else:
        accepted_genes = None

    results = pd.DataFrame(columns=["term_name", "term_id", "parent_ids", "term_description", "p_value", "enrichment",
                                    "num_genes_in_term", "num_genes_in_intersection"])

    if return_genes:
        results["intersected_genes"] = np.nan

    for r in res["result"]:
        N = r["effective_domain_size"]
        B = r["term_size"]
        n = r["query_size"]
        b = r["intersection_size"]
        enrichment = (b / n) / (B / N)
        row = {"term_name": r["name"], "term_id": r["native"], "parent_ids": r["parents"],
               "term_description": r["description"], "p_value": r["p_value"], "enrichment": enrichment,
               "num_genes_in_term": B, "num_genes_in_intersection": b}

        if return_genes:
            intersected_idx = np.array([len(ii) > 0 for ii in r["intersections"]])
            row["intersected_genes"] = list(accepted_genes[intersected_idx])

        results = results.append(row, ignore_index=True)

    return results

                
def query_convert(ids, target_namespace, organism):
    params = {
       'organism': organism,
       'query': "\n".join(ids),
       "target": target_namespace
    }
            
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
        json=params,
        headers={
            'User-Agent': 'FullPythonRequest'
        }
    )

    res = r.json()
    if "result" not in res:
        raise RuntimeError("gProfiler returned an error: " + str(res))
        
    results = [{"incoming": r["incoming"], "converted": r["converted"]} for r in res["result"]]
    return pd.DataFrame.from_records(results)


def query_orthology(ids, from_organism, to_organism, numeric=False):
    #query = "\n".join(ids)
    query = list(ids)

    params = {
        'organism': from_organism,
        'query': query,
        "target": to_organism
    }

    if numeric:
        params["numeric_ns"] = "ENTREZGENE_ACC"

    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/orth/orth/',
        json=params,
        headers={
            'User-Agent': 'FullPythonRequest'
        }
    )

    res = r.json()
    if "result" not in res:
        raise RuntimeError("gProfiler returned an error: " + str(res))

    results = []

    for r in res["result"]:
        results.append((r["incoming"], r["converted"], r["name"]))

    return pd.DataFrame(results, columns=["incoming", "converted_id", "converted_name"])
        

def print_res(res, enrichment_cutoff=None):
    res = list(res)
    seen_children = set()
    for name, desc, p, enrichment, genes, goid, parents in res:
        if enrichment_cutoff is not None and enrichment < enrichment_cutoff:
            continue

        seen_children.update(parents)

    for name, desc, p, enrichment, genes, goid, parents in sorted(res, key=lambda x: x[3], reverse=True):
        if goid in seen_children:
            continue

        if enrichment_cutoff is not None and enrichment < enrichment_cutoff:
            continue

        print("%s: p=%g, enr=%g, num_genes=%d" % (name, p, enrichment, genes))
        print("\t%s" % desc)


def plot_res_enrichr_style(res: pd.DataFrame, top_n=None, fontsize=None):
    seen_children = set()
    for parents in res.parent_ids:
        seen_children.update(parents)

    res = res.loc[~res.term_id.isin(seen_children), :]
    res = res.sort_values(by="enrichment", ascending=False)
    if top_n is not None:
        res = res.iloc[:top_n, :]

    num_results = res.shape[0]

    renderer = plt.gcf().canvas.get_renderer()
    ax = plt.gca()
    inv_transform = ax.transData.inverted()
    bbox = plt.gcf().get_window_extent(renderer=renderer).get_points()

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if num_results == 0:
        ax.spines['bottom'].set_visible(False)
        plt.xticks([])
        plt.yticks([])
        plt.text(*inv_transform.transform(bbox).mean(axis=0), "No significant enrichment", color='black',
                 fontsize=18, verticalalignment='center', ha="center")

        return False

    res["y"] = np.arange(num_results - 1, -1, -1)
    cmap = plt.get_cmap("Wistia")
    c_norm = colors.Normalize(vmin=2, vmax=2 + (res.enrichment.max() - 2) * 1.25)
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cmap)
    plt.barh(res.y, width=res.enrichment, color=scalar_map.to_rgba(res.enrichment))

    if fontsize is None:
        fontsize = 18
        if num_results > 10:
            fontsize = 10

        if num_results > 20:
            fontsize = 8

        if num_results > 30:
            fontsize = 6

    bbox = plt.gcf().get_window_extent(renderer=renderer).get_points()
    for data_row in res.itertuples():
        text = ax.text(0, data_row.y, "  " + data_row.term_name, color='black', fontsize=fontsize,
                       verticalalignment='center', wrap=True)

        bb = text.get_window_extent(renderer=renderer)
        if bb.get_points()[1, 0] > bbox[1, 0] * 0.95:
            text.remove()
            text = ax.text(0, data_row.y, "  " + data_row.term_name, color='black', fontsize=fontsize,
                           verticalalignment='center', wrap=True)

        ax.text(0, data_row.y, str(data_row.num_genes_in_intersection) + " ", color='black',
                fontsize=fontsize, verticalalignment='center', wrap=True, horizontalalignment="right")

    res["p_stars"] = ""
    res.loc[res.p_value < 0.05, "p_stars"] = "*"
    res.loc[res.p_value < 0.01, "p_stars"] = "**"
    res.loc[res.p_value < 0.005, "p_stars"] = "***"

    res.p_stars += "    "

    plt.yticks(res.y, labels=res.p_stars, fontsize=fontsize)
    plt.ylim([-0.8, max(5.0, num_results + 0.8)])

    plt.xlabel("Enrichment")
    return True


def plot_res_table_style(res: pd.DataFrame, top_n=None, fontsize=None, max_str_len=110, remove_parents=True):
    if remove_parents:
        seen_children = set()
        for parents in res.parent_ids:
            seen_children.update(parents)

        res = res.loc[~res.term_id.isin(seen_children), :]

    res = res.sort_values(by="enrichment", ascending=False)
    if top_n is not None:
        res = res.iloc[:top_n, :]

    num_results = res.shape[0]

    renderer = plt.gcf().canvas.get_renderer()
    ax = plt.gca()
    inv_transform = ax.transData.inverted()
    bbox = plt.gcf().get_window_extent(renderer=renderer).get_points()

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if num_results == 0:
        ax.spines['bottom'].set_visible(False)
        plt.xticks([])
        plt.yticks([])
        plt.text(*inv_transform.transform(bbox).mean(axis=0), "No significant enrichment", color='black',
                 fontsize=18, verticalalignment='center', ha="center")

        return 0

    res["y"] = np.arange(num_results - 1, -1, -1)
    cmap = plt.get_cmap("Wistia")
    c_norm = colors.Normalize(vmin=2, vmax=2 + (res.enrichment.max() - 2) * 1.25)
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cmap)
    plt.barh(res.y, width=res.enrichment, color=scalar_map.to_rgba(res.enrichment))

    if fontsize is None:
        fontsize = 18
        if num_results >= 10:
            fontsize = 10

        if num_results > 20:
            fontsize = 8

        if num_results > 30:
            fontsize = 6

    res["p_stars"] = ""
    res.loc[res.p_value < 0.05, "p_stars"] = "*"
    res.loc[res.p_value < 0.01, "p_stars"] = "**"
    res.loc[res.p_value < 0.005, "p_stars"] = "***"

    for data_row in res.itertuples():
        #ax.text(0, data_row.y, f"{data_row.p_stars}  ",
        #        color='black',
        #        fontsize=fontsize, verticalalignment='center', horizontalalignment="left")

        ax.text(data_row.enrichment, data_row.y, f"{data_row.num_genes_in_intersection}",
                color='black',
                fontsize=fontsize, verticalalignment='center', horizontalalignment="left")

    if fontsize < 18:
        max_str_len = int(max_str_len * 18 / fontsize)

    long_terms = res.term_name.str.len() > max_str_len
    res.loc[long_terms, "term_name"] = res.term_name[long_terms].str.slice(0, max_str_len) + "…"

    plt.yticks(res.y, labels=res.term_name, fontsize=fontsize, wrap=True)
    plt.ylim([-0.8, max(5.0, num_results + 0.8)])
    plt.subplots_adjust(left=0.8, right=0.95, top=0.9, bottom=0.2)
    plt.xticks(fontsize=15)

    plt.xlabel("Enrichment", fontsize=18)
    return num_results


def plot_significant_results(genes, libraries: dict, organism, enrichment_cutoff, enrichment_alpha, pics_folder,
                             name_desc=None,
                             suffix="", top_n_to_plot=None, remove_parents=True, plot_style="table", **query_kwargs):
    misc.make_sure_folder_exists(pics_folder)
    results = {}
    for library_id, library_name in libraries.items():
        file = os.path.join(pics_folder, f"{library_id} {suffix}.png")
        if len(genes) == 0:
            print(f"No genes for {file}")
            continue

        print(f"Enrichment analysis for {len(genes)} genes in {library_name}...")
        res = query_gost(genes, library_id, enrichment_alpha, organism=organism, **query_kwargs)
        results[library_id] = res

        fig = plt.figure(figsize=(10, 5))

        if plot_style == "enrichr":
            plot_res = plot_res_enrichr_style
        elif plot_style == "table":
            plot_res = plot_res_table_style
        else:
            assert False, "Unknown plot style " + plot_style

        #print(f"{len(res)} results")

        num_plotted = plot_res(res.loc[res.enrichment >= enrichment_cutoff, ], top_n_to_plot,
                               remove_parents=remove_parents)

        if num_plotted == 0:
            print(f"No data for {file}")
            plt.close(fig)
            continue

        title = library_name
        if top_n_to_plot is not None and num_plotted == top_n_to_plot:
            title += f"\ntop {top_n_to_plot} results by enrichment"

        plt.suptitle(title, fontsize=18)
        plt.savefig(file, dpi=300)
        plt.close(fig)

        if name_desc is not None:
            for data_row in res.itertuples():
                name_desc[library_name][data_row.term_name] = data_row.term_description, data_row.term_id
            
    return results


def plot_significant_by_direction(data, data_name, libraries: dict, organism, de_cutoff, de_alpha, enrichment_cutoff,
                                  enrichment_alpha, pics_folder, name_desc, de_pvalue_column="FDR", top_n_to_plot=None,
                                  plot_style="table",
                                  **query_kwargs):
    misc.make_sure_folder_exists(pics_folder)
    up_genes = data.genes[(data.logFC > de_cutoff) & (data[de_pvalue_column] < de_alpha)].values
    down_genes = data.genes[(data.logFC < -de_cutoff) & (data[de_pvalue_column] < de_alpha)].values
    results = {}
    for dir_data, dir_name in [(up_genes, "up"), (down_genes, "down"),
                               (frozenset(up_genes) | frozenset(down_genes), "both")]:
        if top_n_to_plot is not None and dir_name == "both":
            top_n = top_n_to_plot * 2
        else:
            top_n = top_n_to_plot

        res = plot_significant_results(dir_data, libraries, organism, enrichment_cutoff, enrichment_alpha, pics_folder,
                                       name_desc, suffix=f"_{data_name}_{dir_name}", top_n_to_plot=top_n,
                                       plot_style=plot_style,
                                       **query_kwargs)

        results[dir_name] = res

    return results


def term_genes(terms, namespace, organism):
    res = {t: set() for t in terms}
    for row in query_convert(terms, namespace, organism).itertuples():
        res[row.incoming].add(row.converted)
        
    return res


def write_res_from_libraries(target_worksheet, res_by_library, libraries, enrichment_cutoff, gene_name_map=None,
                             add_columns=None):
    all_libs_res = []
    for lib_id, lib_res in res_by_library.items():
        lib_res["GO library"] = libraries[lib_id]
        lib_res = lib_res.loc[lib_res.enrichment >= enrichment_cutoff]
        if len(lib_res) > 0:
            all_libs_res.append(lib_res)

    if not all_libs_res:
        print("No results to put to worksheet")
        return

    all_libs_res = pd.concat(all_libs_res)
    if add_columns is None:
        add_columns = []

    all_libs_res = all_libs_res.set_index("term_id", drop=False)
    write_res(target_worksheet, all_libs_res, gene_name_map=gene_name_map, add_columns=add_columns + ["GO library"])


def write_res(target_worksheet, res_df, gene_name_map=None, add_columns=None):
    if isinstance(gene_name_map, str):
        gene_name_map = pd.read_csv(gene_name_map, sep="\t")
        gene_name_map = dict(zip(gene_name_map.id, gene_name_map.name))

    res_df = res_df.copy()
    if gene_name_map is not None:
        res_df["affected gene names"] = object()

    for i in range(len(res_df)):
        if res_df.term_id.iloc[i].startswith("KEGG:"):
            res_df.loc[res_df.index[i], "term_description"] = \
                kegg_info.get_kegg_map_description(res_df.term_id.iloc[i][5:])

        if gene_name_map is not None:
            res_df.loc[res_df.index[i], "affected gene names"] = \
                ", ".join((gene_name_map[g] if g in gene_name_map else g) for g in res_df.intersected_genes.iat[i])

        assert isinstance(res_df.intersected_genes.iloc[i], list)
        assert len(res_df.intersected_genes.iloc[i]) == res_df.num_genes_in_intersection.iloc[i]
        res_df.loc[res_df.index[i], "intersected_genes"] = ", ".join(res_df.intersected_genes.iloc[i])

    res_df["% affected genes in term"] = \
        (res_df.num_genes_in_intersection / res_df.num_genes_in_term * 100).astype(float)

    res_df = res_df.round({"enrichment": 2, "% affected genes in term": 2})

    if add_columns is None:
        add_columns = []

    res_df = res_df[add_columns + ["term_name", "term_id", "term_description", "p_value", "enrichment",
                                   "num_genes_in_term",
                                   "num_genes_in_intersection", "% affected genes in term", "intersected_genes"] +
                    ["affected gene names"] * (gene_name_map is not None)]

    res_df = res_df.rename(columns={"term_name": "GO term", "term_id": "GO ID",
                                    "term_description": "GO term description",
                                    "p_value": "P value", "num_genes_in_term": "# genes in term",
                                    "num_genes_in_intersection": "# affected genes in term",
                                    "intersected_genes": "affected genes"})

    if gene_name_map is not None:
        res_df = res_df.rename(columns={"affected genes": "affected gene IDs"})

    misc.write_dataframe_to_worksheet(target_worksheet, res_df)


def map_gene_names(data, gene_id_column, from_organism, to_organism, remove_na=False):
    gene_name_map = query_orthology(data[gene_id_column], organisms[from_organism], organisms[to_organism])
    if remove_na:
        gene_name_map = gene_name_map[gene_name_map.converted_name != "N/A"]

    data = data.merge(gene_name_map, left_on=gene_id_column, right_on="incoming")
    data = data.rename(columns={"converted_name": f"{to_organism}_gene"})
    del data["incoming"]
    del data["converted_id"]
    return data


def go_analysis_for_genes(genes, tag, gene_name_map, wb, alpha, enrichment_cutoff, pics_folder, organism_name,
                          top_n_to_plot=None, remove_parents=True, **query_kwargs):
    res = plot_significant_results(genes, main_libraries, organisms[organism_name],
                                   enrichment_cutoff, alpha,
                                   pics_folder, None, suffix=tag, return_genes=True, top_n_to_plot=top_n_to_plot,
                                   remove_parents=remove_parents, **query_kwargs)

    ws = wb.create_sheet(title=tag[:26] + " GO")
    write_res_from_libraries(ws, res, main_libraries, enrichment_cutoff=enrichment_cutoff, gene_name_map=gene_name_map)
