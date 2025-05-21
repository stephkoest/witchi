# WitChi: A Compositional Bias Pruning Tool for Multiple Sequence Alignments using Chi-squared Statistics
![CI](https://github.com/stephkoest/witchi/actions/workflows/ci.yml/badge.svg)

WitChi is a sophisticated alignment pruning tool designed to reduce compositional bias in multiple sequence alignments (MSAs). It utilizes statistical methods, permutation tests, and parallel processing to iteratively prune alignment columns, optimizing for unbiased taxonomic composition.

![Logo_v2](https://github.com/user-attachments/assets/2af4f0ff-cbbe-48be-a50f-b45de3061b40)

## Key Features
* **Recursive Chi-Square Pruning**: Iteratively removes biased columns based on Chi-square statistics.
* **Multiple Pruning Algorithms**: Supports squared, Wasserstein, and quartic delta Chi-square pruning.
* **Permutation-Based Thresholding**: Empirically estimates unbiased Chi-square score distributions using permutation tests.
* **Parallel Processing**: Leverages multi-threading for computational efficiency.
* **Modular Design**: Perform only pruning or permutation testing depending on your analysis goals.
* **Flexible Configuration**: Adjustable parameters such as pruning depth, top N columns to prune, and permutation count.

![FigO_v1](https://github.com/user-attachments/assets/c68faa7d-07da-4ec6-9382-2d533e63cb64)
Overview of the WitChi workflows for detecting and reducing compositional bias in multiple sequence alignments. (A) The TEST workflow (gray background) computes taxon-specific χ² scores and establishes an empirical null distribution by column permutation (100x), allowing the identification of biased taxa. The PRUNE workflow (green background) iteratively removes alignment columns with the highest Δχ², followed by a convergence check to determine whether pruning should continue or stop. (B) Example observed MSA and one corresponding column- permuted MSA, illustrating how taxon-specific biases are homogenised while preserving global composition. The bar plot on the right compares taxon χ² scores between the observedl and permuted alignments. (C) Left: Δχ² scores per alignment column, with the most biased column flagged for removal (dashed box). Right: Density distributions of taxon-specific χ² scores before and during the pruning loop, showing how pruning shifts scores toward the null distribution. Once no further biased taxa are detected, pruning converges to an unbiased alignment (right panel).

## Installation
**1. Clone the repository:**
```bash
git clone https://github.com/stephkoest/witchi.git
cd witchi
```

**2. Create the Conda Environment:**

```bash
conda env create -f environment.yml
conda activate witchi-env
```
**3. Install locally:**
```bash
pip install -e .
```

**4. Running Tests:**

As an example we use the compositionally biased 5 taxon dataset from Foster et al. 2022 (PMID: 36083446). The dataset is available in the `tests` directory.
Run the tests using unittest:
```bash
python -m unittest discover -s tests -p 'test_witchi.py'
```

## Usage
### Pruning Alignment
Prune alignment columns recursively based on Chi-square test:

```bash
witchi prune --file alignment.fasta --format fasta --max_residue_pruned 100 --permutations 100 --num_workers 4 --top_n 2 --pruning_algorithm quartic
```

#### Options:
- `--file`: Path to the alignment file.
- `--format`: Alignment file format (default: fasta).
- `--max_residue_pruned`: Maximum columns to prune (default: 100).
- `--permutations`: Number of permutations for empirical distribution (default: 100).
- `--num_workers`: Number of CPU threads for parallelization (default: 2).
- `--top_n`: Number of top biased columns to prune per iteration (default: 1).
- `--pruning_algorithm`: Pruning algorithm to use (squared, wasserstein, quartic).

### Permutation Testing
Run permutation tests to establish empirical Chi-square distributions:
```bash
witchi test --file alignment.fasta --format fasta --num_workers 4 --permutations 100
```

#### Options:
- `--file`: Path to the alignment file.
- `--format`: Alignment format (default: fasta).
- `--num_workers`: Number of CPU threads (default: 2).
- `--permutations`: Number of permutations (default: 100).
- `--create_output`: Flag to create output file with z-scores and pseudo p-values per taxon.

## Pruning Algorithms
- **Squared Pruning**: Prioritizes columns with the highest delta Chi-square score.
- **Wasserstein Pruning**: Guides pruning by minimizing Wasserstein distance to the unbiased distribution.
- **Quartic Pruning**: Targets columns that maximize squared taxon delta Chi-square score differences.

## Output
- Pruned Alignment File: A new alignment file with reduced bias.
- Statistical Reports: Summary of the pruning process, including Chi-square distributions and convergence metrics.
- Future: Plots: Visual representation of the pruning impact on Chi-square distributions.

## How It Works
**1. Read Alignment:**

  * The alignment file is parsed and converted to a NumPy array.

**2. Empirical Thresholds:**

  * Permutation tests generate expected Chi-square score distributions.

**3. Pruning Loop:**

  * Iteratively removes the most biased columns based on the selected algorithm.
  * Monitors progress using metrics like Wasserstein distance.

**4. Final Output:**

  * Produces a pruned alignment and statistical summaries.

## Example Workflow
1. Run permutation test:
```bash
witchi test --file example.nex --format "nexus"
```
You can run this on multiple subsets of the dataset, with different taxon sampling for example.
Once tyou have chose datasets you want to prune, go to the next step.

2. Prune up to 50 residues of the alignment with Wasserstein distance guidance:
```bash
witchi prune --file example.nex --max_residue_pruned 50 --pruning_algorithm wasserstein
```

## Future Features
Real-Time Visualization: Live plotting of distribution shifts in the terminal during pruning.

## License
witchi is licensed under the MIT License.

## Publication
Preprint coming soon...

## Contributing
We welcome contributions! Please fork the repository, create a branch, and submit a pull request.

## Acknowledgments
