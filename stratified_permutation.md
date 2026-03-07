# Similarity-Stratified Permutation with Freedman-Lane Recentering

## Method

Standard permutation treats all taxa as exchangeable, which can misrepresent the expected compositional variation on uneven trees. Similarity-stratified permutation addresses this in three steps:

### 1. Stratification

Taxa are assigned to strata based on evolutionary isolation. A seed-projection maps each taxon into a low-dimensional distance space, and mean k-nearest-neighbour distances serve as an isolation metric. Fisher-Jenks natural breaks with Calinski-Harabasz scoring determine the number of strata (capped at K=3). A minimum stratum size of ceil(20N/P)+1 is enforced to guarantee that Bonferroni-corrected significance remains achievable.

### 2. Within-stratum permutation

For each of P permutations, alignment columns are shuffled independently among taxa within each stratum. This preserves the correlation structure imposed by shared evolutionary history while breaking taxon-specific compositional signal.

### 3. Freedman-Lane recentering

After permutation, the null chi-squared for each taxon is computed using within-stratum residuals weighted by global character frequencies:

```
chi2_FL[i] = sum_c ( (n_ci - f_s[c] * L_i)^2 / (f_g[c] * L_i) )
```

where n_ci is the observed count of character c for taxon i, f_s[c] is the stratum frequency, f_g[c] is the global frequency, and L_i is the total character count for taxon i.

Standard (uncorrected) stratified permutation computes chi-squared against the global mean. When a stratum differs compositionally from the global average, this inflates the null — the stratum-level offset is conflated with within-stratum variation. The Freedman-Lane correction replaces the global reference in the numerator with the stratum mean, eliminating both the constant inflation term and its cross-term, while retaining global frequencies in the denominator so that observed and null chi-squared remain on the same scale.

This approach is a direct application of the Freedman-Lane residual permutation procedure (Freedman & Lane, 1983) to compositional data on a phylogenetic tree.

## References

- Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of reported significance levels. *Journal of Business & Economic Statistics*, 1(4), 292-298.
- Anderson, M.J. and Legendre, P. (1999). An empirical comparison of permutation methods for tests of partial regression coefficients in a linear model. *Journal of Statistical Computation and Simulation*, 62, 271-303.
