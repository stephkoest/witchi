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

    has_stratum = any("stratum" in v for v in sorted_dict.values())
    with open(file_name, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        if has_stratum:
            writer.writerow(["Row", "Stratum", "Empirical-Pvalue", "Z-Score"])
            for row, values in sorted_dict.items():
                writer.writerow(
                    [
                        row,
                        values["stratum"],
                        values["empirical_pvalue"],
                        values["zscore"],
                    ]
                )
        else:
            writer.writerow(["Row", "Empirical-Pvalue", "Z-Score"])
            for row, values in sorted_dict.items():
                writer.writerow([row, values["empirical_pvalue"], values["zscore"]])


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
    per_taxon_pools=None,
    name_to_stratum=None,
):
    """Make a dictionary of chi-squared scores for each row in the alignment.

    Parameters
    ----------
    per_taxon_pools : dict or None
        If provided, maps taxon index to its stratum's null pool.
        Each taxon's robust Z-score is computed against its own pool.
        If None, the global pooled null is used for all taxa.
    """
    per_row_chi2 = np.array(per_row_chi2)
    row_names = [record.id for record in alignment]

    zscores = np.array(
        [
            _robust_zscore(
                per_row_chi2[i],
                (
                    per_taxon_pools[i]
                    if per_taxon_pools is not None
                    else permutated_per_row_chi2
                ),
            )
            for i in range(len(per_row_chi2))
        ]
    )

    row_empirical_pvalue_dict = {}
    for i in range(len(row_names)):
        entry = {
            "empirical_pvalue": empirical_pvalues[i],
            "zscore": zscores[i],
        }
        if name_to_stratum is not None:
            entry["stratum"] = name_to_stratum.get(row_names[i], 0)
        row_empirical_pvalue_dict[row_names[i]] = entry
    # sort row_empirical_pvalue_dict by z-score in descending order
    row_empirical_pvalue_dict = dict(
        sorted(
            row_empirical_pvalue_dict.items(),
            key=lambda item: item[1]["zscore"],
            reverse=True,
        )
    )

    return row_empirical_pvalue_dict
