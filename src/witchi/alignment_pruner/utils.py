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
        "Null-Pvalue",
        "Algorithm",
    ]
    with open(file_name, "w", newline="") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for col, values in dictionary.items():
            null_p = values[6] if len(values) > 6 else None
            row = {
                "Iteration": values[0],
                "Original Index": values[1],
                "Global ChiScore": values[2],
                "Initial Significant ChiScore": values[3],
                "ChiScore Difference": values[4],
                "Significant taxa": values[5],
                "Null-Pvalue": "" if null_p is None else null_p,
                "Algorithm": algorithm,
            }
            writer.writerow(row)


def write_score_dict_to_json(dictionary, file_name):
    """Write the score dictionary to a JSON file."""

    def convert_ndarray(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    with open(file_name, "w") as jsonfile:
        json.dump(dictionary, jsonfile, indent=4, default=convert_ndarray)


def write_score_dict_to_tsv(dictionary, file_name):
    """Write the score dictionary to a TSV file, ordered by descending absolute z-score."""
    sorted_dict = dict(
        sorted(
            dictionary.items(),
            key=lambda item: abs(item[1]["zscore"]),
            reverse=True,
        )
    )

    with open(file_name, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["Row", "Empirical-Pvalue", "Z-Score"])
        for row, values in sorted_dict.items():
            writer.writerow([row, values["empirical_pvalue"], values["zscore"]])


def write_before_after_score_dict_to_tsv(dictionary, file_name):
    """Write per-taxon before/after Z-scores and empirical p-values to TSV."""
    with open(file_name, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(
            [
                "Row",
                "Before-Z",
                "Before-Pvalue",
                "After-Z",
                "After-Pvalue",
                "Delta-Z",
            ]
        )
        for row, v in dictionary.items():
            writer.writerow(
                [
                    row,
                    v["before_zscore"],
                    v["before_pvalue"],
                    v["after_zscore"],
                    v["after_pvalue"],
                    v["delta_zscore"],
                ]
            )


def _robust_zscore(observed, null_pool):
    """Compute a robust Z-score using mean centering and MAD scaling.

    Centering on the mean (not median) ensures consistency between
    conditional (per-alignment permutation) and marginal (simulation)
    null distributions, because the mean has the tower property while
    the median does not.
    """
    mean = np.mean(null_pool)
    median = np.median(null_pool)
    mad = np.median(np.abs(null_pool - median))
    if mad == 0:
        return 0.0 if observed == mean else np.sign(observed - mean) * np.inf
    return (observed - mean) / (mad / 0.6745)


def make_score_dict(
    per_row_chi2,
    permutated_per_row_chi2,
    empirical_pvalues,
    alignment,
):
    """Make a dictionary of chi-squared scores for each row in the alignment."""
    per_row_chi2 = np.array(per_row_chi2)
    row_names = [record.id for record in alignment]

    zscores = np.array(
        [
            _robust_zscore(per_row_chi2[i], permutated_per_row_chi2)
            for i in range(len(per_row_chi2))
        ]
    )

    row_empirical_pvalue_dict = {}
    for i in range(len(row_names)):
        row_empirical_pvalue_dict[row_names[i]] = {
            "empirical_pvalue": empirical_pvalues[i],
            "zscore": zscores[i],
        }
    # sort row_empirical_pvalue_dict by z-score in descending order
    row_empirical_pvalue_dict = dict(
        sorted(
            row_empirical_pvalue_dict.items(),
            key=lambda item: item[1]["zscore"],
            reverse=True,
        )
    )

    return row_empirical_pvalue_dict


def make_before_after_score_dict(
    before_chi2,
    after_chi2,
    permutated_per_row_chi2,
    before_pvalues,
    after_pvalues,
    alignment,
):
    """Per-taxon Z and empirical p, before and after pruning, with delta-Z.

    Sorted by descending before-Z (most-biased-initially first).
    """
    row_names = [record.id for record in alignment]
    before = np.asarray(before_chi2)
    after = np.asarray(after_chi2)
    before_z = np.array(
        [_robust_zscore(before[i], permutated_per_row_chi2) for i in range(len(before))]
    )
    after_z = np.array(
        [_robust_zscore(after[i], permutated_per_row_chi2) for i in range(len(after))]
    )

    d = {}
    for i, name in enumerate(row_names):
        d[name] = {
            "before_zscore": float(before_z[i]),
            "before_pvalue": float(before_pvalues[i]),
            "after_zscore": float(after_z[i]),
            "after_pvalue": float(after_pvalues[i]),
            "delta_zscore": float(after_z[i] - before_z[i]),
        }
    return dict(sorted(d.items(), key=lambda kv: kv[1]["before_zscore"], reverse=True))
