import numpy as np
import os
import time
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from joblib import Parallel, delayed

from .chi_square_calculator import (
    ChiSquareCalculator,
    _count_rows_from_int,
    _encode_alignment_int,
)
from .permutation_test import PermutationTest
from .sequence_type_detector import SequenceTypeDetector
from .alignment_reader import AlignmentReader
from witchi import version_banner

from .utils import (
    write_alignment,
    write_pruned_dict_to_tsv,
    write_score_dict_to_json,
    write_before_after_score_dict_to_tsv,
    make_before_after_score_dict,
)


class AlignmentPruner:
    def __init__(
        self,
        file,
        format="fasta",
        max_residue_pruned=100,
        permutations=100,
        num_workers_chisq=1,
        num_workers_permute=1,
        top_n=10,
        pruning_algorithm="wasserstein",
        delta_null=True,
    ):
        self.file = file
        self.format = format
        self.max_residue_pruned = max_residue_pruned
        self.permutations = permutations
        self.num_workers_chisq = num_workers_chisq
        self.num_workers_permute = num_workers_permute
        self.top_n = top_n
        self.pruning_algorithm = pruning_algorithm
        self.delta_null = delta_null
        self.chi_square_calculator = None
        self.alignment_size = None
        self.initial_top_n = self.top_n
        self._null_max_deltas = None
        self._biased_taxa_threshold = None
        self._stop_reason = None
        self._before_stats = None
        self._after_stats = None
        self.prune_strategies = {
            "squared": self._prune_squared,
            "quartic": self._prune_quartic,
            "wasserstein": self._prune_wasserstein,
        }

    def run(self):
        """Main method to run the recursive chi-square pruning."""
        start_time = time.time()
        print("Starting pruning workflow")
        print(f"Chi² calculator using {self.num_workers_chisq} worker(s)")
        print(f"Permutation test using {self.num_workers_permute} worker(s)")
        print(f"Pruning algorithm: {self.pruning_algorithm}")
        print(f"Initial top_n: {self.top_n}")
        reader = AlignmentReader(self.file, self.format)
        alignment, alignment_array = reader.run()

        is_dna, char_set = SequenceTypeDetector.detect(alignment)

        self.chi_square_calculator = ChiSquareCalculator(
            char_set, self.num_workers_chisq
        )
        self.permutation_test = PermutationTest(
            self.num_workers_permute, self.permutations
        )
        (
            sums,
            maxes,
            upper_box_threshold,
            upper_threshold,
            permutated_per_row_chi2,
        ) = self.permutation_test.compute_null(
            alignment_array,
            self.chi_square_calculator,
        )

        # Precompute Z-score parameters from pooled null (consistent with _robust_zscore)
        self._null_z_mean = float(np.mean(permutated_per_row_chi2))
        _median = float(np.median(permutated_per_row_chi2))
        _mad = float(np.median(np.abs(permutated_per_row_chi2 - _median)))
        self._null_z_scale = _mad / 0.6745

        # Compress null Z-scores to quantiles for Wasserstein pruning.
        null_z = (permutated_per_row_chi2 - self._null_z_mean) / self._null_z_scale
        K = min(200, len(null_z))
        positions = np.linspace(0, 1, K + 2)[1:-1]
        self._null_z_quantiles = np.quantile(null_z, positions)

        if self.delta_null:
            print("Computing null delta distribution for stopping criterion...")
            self._compute_null_deltas(alignment_array)

        pruned_alignment_array, prune_dict, score_dict = self.recursive_prune(
            alignment_array,
            sums,
            upper_threshold,
            permutated_per_row_chi2,
        )

        pruned_sequences = self.update_sequences(alignment, pruned_alignment_array)
        pruned_alignment = MultipleSeqAlignment(pruned_sequences)
        # build output filename suffix
        suffix = f"_{self.pruning_algorithm}_s{self.top_n}_pruned.fasta"
        output_alignment_pruned_file = os.path.splitext(self.file)[0] + suffix
        output_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", ".tsv"
        )
        write_alignment(pruned_alignment, output_alignment_pruned_file)
        write_pruned_dict_to_tsv(prune_dict, output_tsv_file, self.pruning_algorithm)

        output_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", ".tsv"
        )
        empirical_pvalues_after = self.permutation_test.calc_empirical_pvalue(
            score_dict["after_real"],
            permutated_per_row_chi2,
        )
        empirical_pvalues_before = self.permutation_test.calc_empirical_pvalue(
            score_dict["before_real"],
            permutated_per_row_chi2,
        )
        before_after_dict = make_before_after_score_dict(
            score_dict["before_real"],
            score_dict["after_real"],
            permutated_per_row_chi2,
            empirical_pvalues_before,
            empirical_pvalues_after,
            alignment,
        )
        output_score_tsv_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", "_scores.tsv"
        )
        write_before_after_score_dict_to_tsv(before_after_dict, output_score_tsv_file)

        # v2 score_dict: reshape pooled null back to (n_perms, n_taxa) and
        # add provenance + the delta-null distribution that drove stopping.
        before_permuted_2d = np.asarray(permutated_per_row_chi2).reshape(
            self.permutations, -1
        )
        score_dict_out = {
            "schema_version": 2,
            "witchi_version": version_banner().replace("witchi ", ""),
            "algorithm": self.pruning_algorithm,
            "stop_reason": self._stop_reason,
            "taxa": [record.id for record in alignment],
            "before_real": score_dict["before_real"],
            "after_real": score_dict["after_real"],
            "before_permuted": before_permuted_2d,
            "null_max_deltas": self._null_max_deltas,
        }
        output_json_file = os.path.splitext(self.file)[0] + suffix.replace(
            ".fasta", "_score_dict.json"
        )
        write_score_dict_to_json(score_dict_out, output_json_file)

        n_taxa = len(alignment)
        n_input_cols = self.alignment_size
        n_removed = len(prune_dict)
        n_output_cols = n_input_cols - n_removed
        pct = (n_removed / n_input_cols) * 100 if n_input_cols else 0.0
        before = self._before_stats or {
            "z": float("nan"),
            "p": float("nan"),
            "biased": 0,
        }
        after = self._after_stats or before
        print(
            f"Pruned {n_input_cols} → {n_output_cols} columns "
            f"(-{n_removed}, {pct:.2f}%) | "
            f"Z {before['z']:.2f} → {after['z']:.2f} | "
            f"p {before['p']:.2f} → {after['p']:.2f} | "
            f"biased {before['biased']}/{n_taxa} → {after['biased']}/{n_taxa} | "
            f"stop: {self._stop_reason}"
        )
        print(f"Pruned alignment saved to {output_alignment_pruned_file}")
        print(f"Columns pruned in order saved to: {output_tsv_file}")
        print(f"Taxa before/after Z and p-values printed to {output_score_tsv_file}")
        print(f"Chiscore values from permutation testing printed to {output_json_file}")

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Execution time: {elapsed_time:.2f} seconds")

    def prune(
        self,
        alignment_array,
        expected_observed,
        count_rows_array,
        permutated_per_row_chi2,
    ):
        try:
            return self.prune_strategies[self.pruning_algorithm](
                alignment_array,
                expected_observed,
                count_rows_array,
                permutated_per_row_chi2,
            )
        except KeyError:
            raise ValueError(f"Unknown pruning algorithm: {self.pruning_algorithm}")

    def _prune_squared(self, alignment_array, expected_observed, count_rows_array, *_):
        initial_global_chi2 = self.chi_square_calculator.calculate_global_chi2(
            expected_observed, count_rows_array
        )
        chi2_differences = self.chi_square_calculator.calculate_global_chi2_difference(
            count_rows_array, self._alignment_int, initial_global_chi2
        )
        return initial_global_chi2, chi2_differences

    def _prune_quartic(self, alignment_array, expected_observed, count_rows_array, *_):
        initial_global_chi2 = (
            self.chi_square_calculator.calculate_quartic_row_global_chi2(
                expected_observed, count_rows_array
            )
        )
        chi2_differences = self.chi_square_calculator.calculate_quartic_chi2_difference(
            count_rows_array, self._alignment_int, initial_global_chi2
        )
        return initial_global_chi2, chi2_differences

    def _prune_wasserstein(
        self,
        alignment_array,
        expected_observed,
        count_rows_array,
        permutated_per_row_chi2,
    ):
        wasserstein = self.chi_square_calculator.calculate_row_zscore_wasserstein(
            expected_observed,
            count_rows_array,
            self._null_z_quantiles,
            self._null_z_mean,
            self._null_z_scale,
        )
        chi2_differences = (
            self.chi_square_calculator.calculate_wasserstein_zscore_difference(
                count_rows_array,
                self._alignment_int,
                wasserstein,
                self._null_z_quantiles,
                self._null_z_mean,
                self._null_z_scale,
            )
        )
        return wasserstein, chi2_differences

    def _compute_null_deltas(self, alignment_array):
        """Compute null distribution of max per-column delta scores.

        For each of P_null permuted alignments, compute per-column deltas
        using the active pruning algorithm, then record the maximum delta.
        The resulting distribution tests whether observed deltas during
        pruning are distinguishable from chance.
        """
        p_null = self.permutations
        char_set = self.chi_square_calculator.char_set

        def compute_one(i):
            iter_seed = 12345 + self.permutations + i
            rng = np.random.default_rng(iter_seed)
            permuted_array = np.apply_along_axis(rng.permutation, 0, alignment_array)
            count_rows = self.chi_square_calculator.calculate_row_counts(permuted_array)
            expected = self.chi_square_calculator.calculate_expected_observed(
                count_rows
            )
            permuted_int = _encode_alignment_int(permuted_array, char_set)
            if self.pruning_algorithm == "squared":
                initial = self.chi_square_calculator.calculate_global_chi2(
                    expected, count_rows
                )
                deltas = self.chi_square_calculator.calculate_global_chi2_difference(
                    count_rows, permuted_int, initial
                )
            elif self.pruning_algorithm == "quartic":
                initial = self.chi_square_calculator.calculate_quartic_row_global_chi2(
                    expected, count_rows
                )
                deltas = self.chi_square_calculator.calculate_quartic_chi2_difference(
                    count_rows, permuted_int, initial
                )
            elif self.pruning_algorithm == "wasserstein":
                initial = self.chi_square_calculator.calculate_row_zscore_wasserstein(
                    expected,
                    count_rows,
                    self._null_z_quantiles,
                    self._null_z_mean,
                    self._null_z_scale,
                )
                deltas = (
                    self.chi_square_calculator.calculate_wasserstein_zscore_difference(
                        count_rows,
                        permuted_int,
                        initial,
                        self._null_z_quantiles,
                        self._null_z_mean,
                        self._null_z_scale,
                    )
                )
            else:
                raise ValueError(f"Unknown pruning algorithm: {self.pruning_algorithm}")
            return max(deltas.values())

        # Avoid nested joblib: outer loop owns all worker capacity here.
        saved_chisq_workers = self.chi_square_calculator.num_workers
        self.chi_square_calculator.num_workers = 1
        try:
            max_deltas_list = Parallel(n_jobs=self.num_workers_permute)(
                delayed(compute_one)(i) for i in range(p_null)
            )
        finally:
            self.chi_square_calculator.num_workers = saved_chisq_workers
        max_deltas = np.asarray(max_deltas_list, dtype=np.float64)
        self._null_max_deltas = max_deltas
        print(
            f"Null max-delta distribution (n={p_null}): "
            f"mean={np.mean(max_deltas):.6f}, "
            f"p95={np.percentile(max_deltas, 95):.6f}"
        )

    def recursive_prune(
        self,
        alignment_array,
        sums,
        upper_threshold,
        permutated_per_row_chi2,
    ):
        """Recursively prune the alignment array."""
        prune_dict = {}
        score_dict = {}
        iteration = 0
        removed_columns_count = 0
        self.alignment_size = alignment_array.shape[1]
        original_indices = list(range(alignment_array.shape[1]))

        from .utils import _robust_zscore

        self._biased_taxa_threshold = self._compute_biased_taxa_threshold(
            permutated_per_row_chi2, alignment_array.shape[0]
        )

        top_n_indices = []
        initial_global_chi2 = None
        chi2_differences = {}
        null_pvals_for_batch = {}

        # Reactive touchdown rollback state: keep a reference to the
        # original alignment so we can reconstruct it (via prune_dict)
        # when rolling back the overshoot batch. in_touchdown flips True
        # after the first rollback so we don't rollback a second time.
        original_alignment_array = alignment_array
        in_touchdown = False

        # Across-iteration caches for the inner chi^2 path: encoded
        # alignment (int8 codepoints) and un-fudged char counts. Both are
        # maintained incrementally as columns are removed so each
        # iteration avoids re-encoding the alignment and re-counting
        # characters. Rebuilt on touchdown rollback.
        char_set = self.chi_square_calculator.char_set
        n_chars = len(char_set)
        self._alignment_int = _encode_alignment_int(alignment_array, char_set)
        self._count_rows_raw = _count_rows_from_int(self._alignment_int, n_chars)

        while removed_columns_count < self.max_residue_pruned:
            stats = self._calculate_per_row_stats(alignment_array)
            count_rows_array = stats["count_rows"]
            expected_observed = stats["expected_observed"]
            per_row_chi2 = stats["per_row_chi2"]
            global_chi2 = stats["sum"]

            alignment_empirical_p, significant_count = self._calc_empirical_pvals(
                stats, sums, permutated_per_row_chi2
            )
            alignment_z = _robust_zscore(np.sum(per_row_chi2), sums)

            if iteration == 0:
                score_dict["before_permuted"] = permutated_per_row_chi2
                score_dict["before_real"] = per_row_chi2
                self._before_stats = {
                    "z": alignment_z,
                    "p": alignment_empirical_p,
                    "biased": significant_count,
                }
            else:
                # only write to prune_dict after the first pruning iteration
                for col in top_n_indices:
                    pruned_col = original_indices[col]
                    prune_dict[pruned_col] = [
                        iteration,
                        pruned_col,
                        global_chi2,
                        initial_global_chi2,
                        chi2_differences[col],
                        significant_count,
                        null_pvals_for_batch.get(col),
                    ]
                    removed_columns_count += 1

                original_indices = [
                    i for j, i in enumerate(original_indices) if j not in top_n_indices
                ]

            should_stop, stop_reason = self._check_stopping_criteria(
                alignment_empirical_p
            )

            print(
                f"Columns removed: {removed_columns_count}, "
                f"{(removed_columns_count / self.alignment_size) * 100:.2f}% | "
                f"Biased taxa permutation: {significant_count} | "
                f"Alignment Z-score: {alignment_z:.2f} | "
                f"Alignment p-value: {alignment_empirical_p:.2f}"
            )

            # Reactive touchdown rollback: undo the most recent batch
            # and continue with a 10x smaller top_n. Bounds the overshoot
            # from either signal (alignment-level stop OR delta-null
            # partial batch) and gives the adaptive walk room to engage
            # in the final small batches. Fires at most once per run.
            target_top_n = max(1, self.initial_top_n // 10)
            if (
                should_stop
                and not in_touchdown
                and self.top_n > target_top_n
                and len(top_n_indices) > 0
            ):
                alignment_array, original_indices, n_recent = (
                    self._do_touchdown_rollback(
                        top_n_indices,
                        prune_dict,
                        original_alignment_array,
                        target_top_n,
                        char_set,
                        n_chars,
                    )
                )
                removed_columns_count -= n_recent
                in_touchdown = True
                top_n_indices = []

                print(
                    f"Touchdown: alignment-level overshoot "
                    f"(reason='{stop_reason}', "
                    f"alignment_p={alignment_empirical_p:.3f}). "
                    f"Rolled back {n_recent} columns, reduced top_n to "
                    f"{target_top_n}, resuming pruning."
                )
                continue

            if should_stop:
                self._stop_reason = stop_reason
                self._after_stats = {
                    "z": alignment_z,
                    "p": alignment_empirical_p,
                    "biased": significant_count,
                }
                break

            initial_global_chi2, chi2_differences = self.prune(
                alignment_array,
                expected_observed,
                count_rows_array,
                permutated_per_row_chi2,
            )

            # Sort candidates in descending delta order, take up to top_n
            sorted_candidates = sorted(
                chi2_differences.items(), key=lambda x: x[1], reverse=True
            )[: self.top_n]

            # Delta null: walk candidates and remove those above the null max.
            # Stops at first failure within the batch, giving graceful
            # touchdown — batch size naturally shrinks as we approach
            # convergence. Exits the pruning loop only when rank-1 itself
            # fails (zero columns justified in this iteration).
            if self._null_max_deltas is not None:
                justified = []
                null_pvals_for_batch = {}
                last_p = None
                for col, delta in sorted_candidates:
                    last_p = np.sum(self._null_max_deltas >= delta) / len(
                        self._null_max_deltas
                    )
                    if last_p > 0.05:
                        break
                    justified.append(col)
                    null_pvals_for_batch[col] = float(last_p)

                if not justified:
                    max_delta = sorted_candidates[0][1]
                    self._stop_reason = (
                        f"delta not significant (p={last_p:.3f}, "
                        f"max_delta={max_delta:.6f})"
                    )
                    self._after_stats = {
                        "z": alignment_z,
                        "p": alignment_empirical_p,
                        "biased": significant_count,
                    }
                    break

                if len(justified) < len(sorted_candidates):
                    fail_delta = sorted_candidates[len(justified)][1]
                    # Harmonised touchdown: a partial batch means the
                    # current top_n was too aggressive. If the touchdown
                    # budget is unspent and a previous batch exists, roll
                    # it back and resume with reduced top_n (same
                    # mechanism as the alignment-level overshoot path).
                    if (
                        not in_touchdown
                        and self.top_n > target_top_n
                        and len(top_n_indices) > 0
                    ):
                        alignment_array, original_indices, n_recent = (
                            self._do_touchdown_rollback(
                                top_n_indices,
                                prune_dict,
                                original_alignment_array,
                                target_top_n,
                                char_set,
                                n_chars,
                            )
                        )
                        removed_columns_count -= n_recent
                        in_touchdown = True
                        top_n_indices = []
                        print(
                            f"Touchdown: delta-null partial batch "
                            f"(rank {len(justified) + 1} fails at "
                            f"p={last_p:.3f}, delta={fail_delta:.6f}). "
                            f"Rolled back {n_recent} columns, reduced "
                            f"top_n to {target_top_n}, resuming pruning."
                        )
                        continue
                    print(
                        f"  Partial batch: removing {len(justified)}/"
                        f"{len(sorted_candidates)} columns "
                        f"(rank {len(justified) + 1} fails at "
                        f"p={last_p:.3f}, delta={fail_delta:.6f})"
                    )
            else:
                justified = [col for col, _ in sorted_candidates]

            top_n_indices = np.array(justified)
            # Decrement cached raw counts by the removed columns'
            # contributions, then drop the columns from both the string
            # and int8 forms in lockstep.
            self._count_rows_raw -= _count_rows_from_int(
                self._alignment_int[:, top_n_indices], n_chars
            )
            alignment_array = np.delete(alignment_array, top_n_indices, axis=1)
            self._alignment_int = np.delete(self._alignment_int, top_n_indices, axis=1)
            iteration += 1

        # Natural exit (max_residue_pruned cap): the loop top-check fires
        # AFTER the last iteration has already pruned its batch but BEFORE
        # the next iteration's else-branch records it. Recover by recording
        # the unlogged final batch and refreshing per-row stats so
        # after_real matches the on-disk pruned alignment.
        if (
            self._stop_reason is None
            and top_n_indices is not None
            and len(top_n_indices) > 0
        ):
            final_stats = self._calculate_per_row_stats(alignment_array)
            final_p, final_biased = self._calc_empirical_pvals(
                final_stats, sums, permutated_per_row_chi2
            )
            for col in top_n_indices:
                pruned_col = original_indices[col]
                prune_dict[pruned_col] = [
                    iteration,
                    pruned_col,
                    final_stats["sum"],
                    initial_global_chi2,
                    chi2_differences[col],
                    final_biased,
                    null_pvals_for_batch.get(col),
                ]
                removed_columns_count += 1
            per_row_chi2 = final_stats["per_row_chi2"]
            self._stop_reason = (
                f"reached --max_residue_pruned cap ({self.max_residue_pruned})"
            )
            self._after_stats = {
                "z": _robust_zscore(np.sum(per_row_chi2), sums),
                "p": final_p,
                "biased": final_biased,
            }
        elif self._stop_reason is None:
            # Loop never entered (cap=0) or exited with no prior pruning.
            self._stop_reason = (
                f"reached --max_residue_pruned cap ({self.max_residue_pruned})"
            )
            if self._after_stats is None:
                self._after_stats = self._before_stats

        score_dict["after_real"] = per_row_chi2
        self.top_n = self.initial_top_n
        return alignment_array, prune_dict, score_dict

    def _check_stopping_criteria(self, alignment_empirical_p):
        """Stop on alignment p-value: not-biased threshold (0.05) without
        delta-null; over-prune floor (0.95) with it, where the per-column
        walk supplies the evidence to prune past not-biased."""
        threshold = 0.95 if self.delta_null else 0.05
        if alignment_empirical_p > threshold:
            return True, "alignment p-value"
        return False, None

    def _compute_biased_taxa_threshold(self, permuted_chi2, n_taxa):
        """Chi^2 cutoff above which a taxon counts as biased.

        Equivalent to the per-taxon Bonferroni empirical test
        clip(#{null >= obs} / M * n_taxa, 0, 1) <= 0.05, which is monotone
        in obs and so collapses to obs > tau. With M = len(pooled null) and
        k = floor(0.05 * M / n_taxa) the allowed null exceedances, tau is
        the (M - k - 1)-th order statistic of the pooled null: a taxon is
        biased iff its chi^2 strictly exceeds tau (exact in all tie cases).
        Returns +inf when no taxon can qualify.
        """
        pooled = np.asarray(permuted_chi2).ravel()
        M = pooled.size
        if M == 0 or n_taxa <= 0:
            return np.inf
        k = int(np.floor(0.05 * M / n_taxa))
        kth = M - k - 1
        if kth < 0:
            return np.inf
        return float(np.partition(pooled, kth)[kth])

    def _calc_empirical_pvals(self, stats, permuted_sums, permuted_chi2):
        """Alignment-level empirical p-value plus the biased-taxa diagnostic.

        The alignment p-value (scalar, drives stopping) is computed exactly.
        The biased-taxa count uses the precomputed chi^2 cutoff
        (self._biased_taxa_threshold), which is mathematically identical to
        the per-taxon Bonferroni empirical test but O(N) rather than
        O(P*N^2). Falls back to the exact per-taxon computation if the
        threshold was not precomputed.
        """
        alignment_empirical_p = self.permutation_test.calc_empirical_pvalue(
            stats["sum"], permuted_sums
        )[0]
        per_row_chi2 = stats["per_row_chi2"]
        if self._biased_taxa_threshold is not None:
            significant_count = int(
                np.sum(np.asarray(per_row_chi2) > self._biased_taxa_threshold)
            )
        else:
            empirical_pvals = self.permutation_test.calc_empirical_pvalue(
                per_row_chi2,
                permuted_chi2,
            )
            significant_count = sum(p <= 0.05 for p in empirical_pvals)
        return alignment_empirical_p, significant_count

    def _do_touchdown_rollback(
        self,
        top_n_indices,
        prune_dict,
        original_alignment_array,
        target_top_n,
        char_set,
        n_chars,
    ):
        """Roll back the most recent applied batch and reduce top_n.

        Mutates `prune_dict` (deletes the last batch's entries) and
        `self.top_n`, `self._alignment_int`, `self._count_rows_raw`.
        Returns (alignment_array, original_indices, n_recent) for the
        caller to install.
        """
        n_recent = len(top_n_indices)
        recent_cols = list(prune_dict.keys())[-n_recent:]
        for col in recent_cols:
            del prune_dict[col]

        cols_to_remove = sorted(prune_dict.keys())
        if cols_to_remove:
            alignment_array = np.delete(
                original_alignment_array, cols_to_remove, axis=1
            )
        else:
            alignment_array = original_alignment_array
        original_indices = [
            i
            for i in range(original_alignment_array.shape[1])
            if i not in set(cols_to_remove)
        ]
        # Rebuild the cached encoded alignment + raw counts to match
        # the rolled-back string alignment. Touchdown fires at most
        # once per run, so the from-scratch cost is fine.
        self._alignment_int = _encode_alignment_int(alignment_array, char_set)
        self._count_rows_raw = _count_rows_from_int(self._alignment_int, n_chars)
        self.top_n = target_top_n
        return alignment_array, original_indices, n_recent

    def _calculate_per_row_stats(self, alignment_array):
        """Calculate row-wise chi² statistics and summary metrics."""
        count_rows_array = self._cached_count_rows()
        expected_observed = self.chi_square_calculator.calculate_expected_observed(
            count_rows_array
        )
        per_row_chi2 = self.chi_square_calculator.calculate_row_chi2(
            expected_observed, count_rows_array
        )

        return {
            "count_rows": count_rows_array,
            "expected_observed": expected_observed,
            "per_row_chi2": per_row_chi2,
            "median": np.median(per_row_chi2),
            "mean": np.mean(per_row_chi2),
            "sum": np.sum(per_row_chi2),
        }

    def _cached_count_rows(self):
        """Apply the +1 zero-safety fudge fresh to the cached raw counts.

        The fudge is a step function on the raw counts (fires whenever any
        cell is zero), so it is re-evaluated each iteration rather than
        incrementally maintained.
        """
        if (self._count_rows_raw == 0).any():
            return self._count_rows_raw + 1
        return self._count_rows_raw.copy()

    def update_sequences(self, alignment, pruned_alignment_array):
        """Update sequences with pruned alignments."""
        pruned_sequences = []

        for i, record in enumerate(alignment):
            pruned_seq = "".join(pruned_alignment_array[i, :])
            record.seq = Seq(pruned_seq)
            pruned_sequences.append(record)

        return pruned_sequences
