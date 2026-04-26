import numpy as np
from joblib import Parallel, delayed


def _encode_alignment_int(alignment_array, char_set):
    """(T, L) string alignment -> (T, L) int8: char index, or -1 for gap.

    Codepoint LUT: one indexing pass over the uint32 view of the '<U1'
    array replaces 20 string compares per column.
    """
    arr_c = np.ascontiguousarray(alignment_array)
    codepoints = arr_c.view(np.uint32).reshape(arr_c.shape) & 0xFF
    table = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(char_set):
        table[ord(c)] = i
    return table[codepoints]


def _count_rows_from_int(alignment_int, n_chars):
    """(C, T) un-fudged char counts from int-encoded alignment.

    Vectorised bincount over (c * T + t) for nongap entries. Equivalent
    to ChiSquareCalculator.calculate_row_counts without the +1 zero-safety
    fudge — callers that need the fudge should re-apply it.
    """
    T = alignment_int.shape[0]
    nongap = alignment_int >= 0
    chars = alignment_int[nongap].astype(np.int64)
    flat_t = np.broadcast_to(np.arange(T)[:, None], alignment_int.shape)[nongap]
    return (
        np.bincount(chars * T + flat_t, minlength=n_chars * T)
        .reshape(n_chars, T)
        .astype(np.int64)
    )


def _per_row_chi2_block(
    cr,
    alignment_int,
    U,
    total_chars_base,
    N_base,
    col_totals_base,
    j_start,
    j_end,
):
    """Per-row chi^2 for one column block, when each block-column is removed.

    Marginal identity: e[c,t] = gf_j[c] * col_totals_j[t]
        per_row[t] = S[t] / col_totals_j[t] - col_totals_j[t]
        S[t]       = (U.T @ inv_gf)[t] + nongap[t] * (1 - 2*cr_at_t) * inv_gf_at_t
    with U[c,t] = cr[c,t]^2, inv_gf[c] = N_j / total_chars_j[c].

    Returns (per_row_block (T, b), zero_mask (b,)) where zero_mask flags
    block-columns with degenerate totals — caller falls back to the direct
    per-col path for those.
    """
    C, T = cr.shape
    b = j_end - j_start
    block = alignment_int[:, j_start:j_end]
    nongap = block >= 0
    nongap_f = nongap.astype(np.float64)
    chars_s = np.where(nongap, block, np.int8(0))

    flat_b = np.broadcast_to(np.arange(b, dtype=np.int64), (T, b))
    count_at_col = (
        np.bincount(
            chars_s[nongap].astype(np.int64) * b + flat_b[nongap],
            minlength=C * b,
        )
        .reshape(C, b)
        .astype(np.float64)
    )

    n_nongap = nongap_f.sum(axis=0)
    N_j = N_base - n_nongap
    tc_j = total_chars_base[:, None] - count_at_col
    col_totals_j = col_totals_base[:, None] - nongap_f

    zero_mask = (tc_j <= 0).any(axis=0) | (col_totals_j <= 0).any(axis=0)
    safe_tc = np.where(tc_j > 0, tc_j, 1.0)
    inv_gf = N_j[None, :] / safe_tc
    inv_gf_at_tb = inv_gf[chars_s, flat_b]
    cr_at_tb = cr.T[np.arange(T)[:, None], chars_s]

    A = U.T @ inv_gf
    S = A + nongap_f * (1.0 - 2.0 * cr_at_tb) * inv_gf_at_tb
    safe_ct = np.where(col_totals_j > 0, col_totals_j, 1.0)
    return S / safe_ct - col_totals_j, zero_mask


def _reduce_chunk(
    cr,
    alignment_int,
    U,
    total_chars_base,
    N_base,
    col_totals_base,
    chunk_start,
    chunk_end,
    B,
    reducer,
):
    """Process one column chunk by iterating L2-sized blocks within it.

    Decouples the joblib dispatch unit (chunk) from the L2 cache unit
    (block): one chunk = one joblib task, but inside the chunk many
    blocks run serially so per-block ufuncs stay cache-resident. Returns
    ((chunk_size,) reduced, list of zero-col indices in absolute coords).
    """
    out = np.empty(chunk_end - chunk_start, dtype=np.float64)
    zero_indices = []
    for j_start in range(chunk_start, chunk_end, B):
        j_end = min(j_start + B, chunk_end)
        per_row_block, zero_mask = _per_row_chi2_block(
            cr,
            alignment_int,
            U,
            total_chars_base,
            N_base,
            col_totals_base,
            j_start,
            j_end,
        )
        out[j_start - chunk_start : j_end - chunk_start] = reducer(per_row_block)
        if zero_mask.any():
            zero_indices.extend(int(j_start + k) for k in np.where(zero_mask)[0])
    return out, zero_indices


def _reduce_sum(per_row_block):
    return per_row_block.sum(axis=0)


def _reduce_sum_of_squares(per_row_block):
    return (per_row_block**2).sum(axis=0)


class _WassersteinReducer:
    """Picklable per-block reducer for the wasserstein difference."""

    __slots__ = ("positions", "null_q")

    def __init__(self, positions, null_q):
        self.positions = positions
        self.null_q = null_q

    def __call__(self, per_row_block):
        obs_q = np.quantile(per_row_block, self.positions, axis=0)
        return np.mean(np.abs(obs_q - self.null_q), axis=0)


class _WassersteinZReducer:
    """Picklable per-block reducer for the wasserstein-Z difference."""

    __slots__ = ("positions", "null_q", "mean", "scale")

    def __init__(self, positions, null_q, mean, scale):
        self.positions = positions
        self.null_q = null_q
        self.mean = mean
        self.scale = scale

    def __call__(self, per_row_block):
        z = (per_row_block - self.mean) / self.scale
        obs_q = np.quantile(z, self.positions, axis=0)
        return np.mean(np.abs(obs_q - self.null_q), axis=0)


class ChiSquareCalculator:
    def __init__(self, char_set, num_workers=2):
        self.char_set = char_set
        self.num_workers = num_workers

    def calculate_row_counts(self, alignment_array):
        """
        Count occurrences of all characters in the provided char_set for each column
        of the alignment array, ignoring gaps ('-').
        """
        count_rows_array = np.array(
            [(alignment_array == char).sum(axis=1) for char in self.char_set]
        )
        if np.any(count_rows_array == 0):
            # Avoid division by zero
            count_rows_array = count_rows_array + 1

        return count_rows_array

    def calculate_expected_observed(self, count_rows_array):
        """
        Fully vectorized calculation of expected and observed frequencies for each column in
        the alignment using precomputed character counts.
        """
        total_chars = np.sum(count_rows_array, axis=1)
        total_chars_sum = np.sum(total_chars)
        global_frequencies = total_chars / total_chars_sum
        column_totals = np.sum(count_rows_array, axis=0, keepdims=True)
        expected_frequencies = global_frequencies[:, np.newaxis] * column_totals
        if np.any(expected_frequencies == 0):
            # Avoid division by zero
            expected_frequencies = expected_frequencies + 1

        return expected_frequencies

    def calculate_global_chi2(self, expected_values, count_rows_array):
        """Calculate the global chi-squared score for the entire alignment."""
        # Calculate chi-squared values
        chi2_values = (count_rows_array - expected_values) ** 2 / expected_values
        total_chi2 = np.sum(chi2_values)

        return total_chi2

    def calculate_row_chi2(self, expected_values, count_rows_array):
        """Calculate the chi-squared score for all taxa."""
        chi2_values = (count_rows_array - expected_values) ** 2 / expected_values
        chi2_sum = np.sum(chi2_values, axis=0)
        return chi2_sum

    def calculate_row_chi2_wasserstein(
        self, expected_values, count_rows_array, null_quantiles
    ):
        """Wasserstein-1 distance between observed per-taxon chi² and null.

        Parameters
        ----------
        null_quantiles : (K,) array
            Pre-compressed quantiles of the null distribution at K
            uniformly-spaced positions.
        """
        per_row_chi2 = self.calculate_row_chi2(expected_values, count_rows_array)
        K = len(null_quantiles)
        positions = np.linspace(0, 1, K + 2)[1:-1]
        obs_quantiles = np.quantile(per_row_chi2, positions)
        return float(np.mean(np.abs(obs_quantiles - null_quantiles)))

    def calculate_quartic_row_global_chi2(self, expected_values, count_rows_array):
        """Calculate the chi-squared score for all taxa."""
        chi2_values = (count_rows_array - expected_values) ** 2 / expected_values
        chi2_sum = np.sum(chi2_values, axis=0)
        quartic_chi2_sum = np.sum(chi2_sum**2)

        return quartic_chi2_sum

    def _reduce_chi2_per_col_removed(self, count_rows_array, alignment_int, reducer):
        """Stream per-row chi^2 blocks through `reducer`; return (L,) per-col.

        Two-level work splitting:
          * inner block (B, sized for L2 cache) — locality unit, never
            materialises (T, L); each block is reduced to (b,) and dropped.
          * outer chunk (sized for joblib amortisation) — when
            num_workers > 1, the column range is split into n_chunks each
            covering many blocks; one chunk = one joblib task, so the
            per-task overhead amortises across multiple blocks of work.
        Degenerate columns fall back to a direct per-col path on the
        int8 column slice with the same reducer.
        """
        cr = count_rows_array.astype(np.float64, copy=False)
        C, T = cr.shape
        L = alignment_int.shape[1]

        total_chars_base = cr.sum(axis=1)
        N_base = float(total_chars_base.sum())
        col_totals_base = cr.sum(axis=0)
        U = cr * cr

        B = max(16, min(256, 2_000_000 // (8 * T)))
        # Decide whether to dispatch chunks via joblib. Two guards:
        #   * each chunk must hold at least ~4 blocks (cache fit + overhead).
        #   * total kernel work T*L must exceed ~5e6 to amortise the
        #     ~50-100 ms per-task loky cost. On small alignments
        #     (e.g. T*L ~ a few million) the parallel branch is skipped
        #     entirely and we run a single serial chunk.
        if T * L < 5_000_000 or self.num_workers <= 1:
            n_chunks = 1
        else:
            n_chunks = min(self.num_workers, max(1, L // B // 4))

        if n_chunks > 1:
            bounds = np.linspace(0, L, n_chunks + 1, dtype=int)
            chunk_ranges = list(zip(bounds[:-1], bounds[1:]))
            results = Parallel(n_jobs=n_chunks)(
                delayed(_reduce_chunk)(
                    cr,
                    alignment_int,
                    U,
                    total_chars_base,
                    N_base,
                    col_totals_base,
                    cs,
                    ce,
                    B,
                    reducer,
                )
                for cs, ce in chunk_ranges
            )
        else:
            chunk_ranges = [(0, L)]
            results = [
                _reduce_chunk(
                    cr,
                    alignment_int,
                    U,
                    total_chars_base,
                    N_base,
                    col_totals_base,
                    0,
                    L,
                    B,
                    reducer,
                )
            ]

        out = np.empty(L, dtype=np.float64)
        zero_cols = []
        for (cs, ce), (chunk_out, zero_in_chunk) in zip(chunk_ranges, results):
            out[cs:ce] = chunk_out
            zero_cols.extend(zero_in_chunk)

        for j in zero_cols:
            single_col = _count_rows_from_int(alignment_int[:, j : j + 1], C)
            cc = count_rows_array - single_col
            per_row_single = self.calculate_row_chi2(
                self.calculate_expected_observed(cc), cc
            )
            out[j] = float(reducer(per_row_single[:, None])[0])

        return out

    def calculate_global_chi2_difference(
        self, count_rows_array, alignment_int, initial_global_chi2
    ):
        """Chi-squared difference per column: how much the global chi^2 drops if column j is removed."""
        deltas = self._reduce_chi2_per_col_removed(
            count_rows_array, alignment_int, _reduce_sum
        )
        return dict(enumerate(initial_global_chi2 - deltas))

    def calculate_wasserstein_difference(
        self, count_rows_array, alignment_int, wasserstein, null_quantiles
    ):
        """Wasserstein-1 difference per column against the chi^2 null quantiles."""
        K = len(null_quantiles)
        positions = np.linspace(0, 1, K + 2)[1:-1]
        null_q = np.asarray(null_quantiles, dtype=np.float64)[:, None]
        reducer = _WassersteinReducer(positions, null_q)
        deltas = self._reduce_chi2_per_col_removed(
            count_rows_array, alignment_int, reducer
        )
        return dict(enumerate(wasserstein - deltas))

    def calculate_row_zscore_wasserstein(
        self,
        expected_values,
        count_rows_array,
        null_z_quantiles,
        null_z_mean,
        null_z_scale,
    ):
        """Wasserstein-1 distance between observed per-taxon Z-scores and null Z distribution.

        Converts per-taxon chi² to Z-scores using the pooled null parameters,
        then compares quantiles against the null Z quantiles.
        """
        per_row_chi2 = self.calculate_row_chi2(expected_values, count_rows_array)
        z_scores = (per_row_chi2 - null_z_mean) / null_z_scale
        K = len(null_z_quantiles)
        positions = np.linspace(0, 1, K + 2)[1:-1]
        obs_z_quantiles = np.quantile(z_scores, positions)
        return float(np.mean(np.abs(obs_z_quantiles - null_z_quantiles)))

    def calculate_wasserstein_zscore_difference(
        self,
        count_rows_array,
        alignment_int,
        wasserstein,
        null_z_quantiles,
        null_z_mean,
        null_z_scale,
    ):
        """Wasserstein-1 difference per column against the null Z-score quantiles."""
        K = len(null_z_quantiles)
        positions = np.linspace(0, 1, K + 2)[1:-1]
        null_q = np.asarray(null_z_quantiles, dtype=np.float64)[:, None]
        reducer = _WassersteinZReducer(positions, null_q, null_z_mean, null_z_scale)
        deltas = self._reduce_chi2_per_col_removed(
            count_rows_array, alignment_int, reducer
        )
        return dict(enumerate(wasserstein - deltas))

    def calculate_quartic_chi2_difference(
        self, count_rows_array, alignment_int, initial_global_chi2
    ):
        """Chi-squared difference per column on the quartic (sum of per-row chi^2 squared) score."""
        deltas = self._reduce_chi2_per_col_removed(
            count_rows_array, alignment_int, _reduce_sum_of_squares
        )
        return dict(enumerate(initial_global_chi2 - deltas))
