import csv
import json
from Bio import AlignIO
import numpy as np


def write_alignment(alignment, output_file, format="fasta"):
    with open(output_file, "w") as output_handle:
        AlignIO.write(alignment, output_handle, format)


def write_pruned_dict_to_tsv(dictionary, file_name, algorithm):
    headers = [
        "Iteration",
        "Original Index",
        "Global ChiScore",
        "Initial Significant ChiScore",
        "ChiScore Difference",
        "Significant taxa",
        "Algorithm",
    ]
    with open(file_name, "w", newline="") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for col, values in dictionary.items():
            row = {
                "Iteration": values[0],
                "Original Index": values[1],
                "Global ChiScore": values[2],
                "Initial Significant ChiScore": values[3],
                "ChiScore Difference": values[4],
                "Significant taxa": values[5],
                "Algorithm": algorithm,
            }
            writer.writerow(row)


def write_score_dict_to_json(dictionary, file_name):
    def convert_ndarray(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    with open(file_name, "w") as jsonfile:
        json.dump(dictionary, jsonfile, indent=4, default=convert_ndarray)


def write_score_dict_to_tsv(dictionary, file_name):
    sorted_dict = dict(
        sorted(
            dictionary.items(), key=lambda item: abs(item[1]["zscore"]), reverse=True
        )
    )
    with open(file_name, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["Row", "Empirical-Pvalue", "Z-Score"])
        for row, values in sorted_dict.items():
            writer.writerow([row, values["empirical_pvalue"], values["zscore"]])
