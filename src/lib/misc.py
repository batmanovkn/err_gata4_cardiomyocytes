import os
from collections import Counter
import tempfile

from openpyxl.utils.dataframe import dataframe_to_rows


def remove_extension(file_name):
    if file_name.endswith(".gz"):
        file_name = file_name.replace(".gz", "")

    return os.path.join(os.path.dirname(file_name), ".".join(os.path.basename(file_name).split(".")[:-1]))


def replace_extension(file_name, new_extension):
    return remove_extension(file_name) + "." + new_extension


def in_folder(file_name, folder):
    return os.path.join(folder, os.path.basename(file_name))


def make_sure_folder_exists(folder):
    folder = os.path.abspath(folder)
    if not os.path.exists(folder):
        make_sure_folder_exists(os.path.dirname(folder))
        os.mkdir(folder)


def write_dataframe_to_worksheet(ws, df):
    for r in dataframe_to_rows(df, index=False, header=True):
        ws.append(r)

    for cell in ws[1]:
        cell.style = 'Headline 4'


def rename_file(from_path, to_path, runner):
    runner.add(f"mv '{from_path}' '{to_path}'")


def copy_file(from_path, to_path, runner):
    runner.add(f"cp '{from_path}' '{to_path}'")


def join_quoted_paths(paths):
    if isinstance(paths, str):
        return f"'{paths}'"

    if isinstance(paths, list):
        return ' '.join(join_quoted_paths(p) for p in paths)

    assert False


def delete_files(paths, runner):
    runner.add(f"rm {join_quoted_paths(paths)}")


def make_values_unique(d):
    counter = Counter(d.values())
    seen = Counter()
    for key in d:
        if counter[d[key]] > 1:
            seen[d[key]] += 1
            d[key] += f"_{seen[d[key]]}"


class TempFileName:
    def __init__(self) -> None:
        with tempfile.NamedTemporaryFile(delete=False) as f:
            self.name = f.name

    def __enter__(self):
        return self.name

    def __exit__(self, *args):
        os.remove(self.name)
