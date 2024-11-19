# WitChi: A Compositional Bias Pruning Tool for Multiple Sequence Alignments
WitChi is a sophisticated alignment pruning tool designed to reduce compositional bias in multiple sequence alignments (MSAs). It utilizes statistical methods, permutation tests, and parallel processing to iteratively prune alignment columns, optimizing for unbiased taxonomic composition.
![Logo_v2](https://github.com/user-attachments/assets/e8f43248-0904-4414-afeb-ffaa8cf01f27)


## Key Features
* **Recursive Chi-Square Pruning**: Iteratively removes biased columns based on Chi-square statistics.
* **Multiple Pruning Algorithms**: Supports global, outlyingness, Wasserstein, and squared delta Chi-square pruning.
* **Permutation-Based Thresholding**: Empirically estimates unbiased Chi-square score distributions using permutation tests.
* **Parallel Processing**: Leverages multi-threading for computational efficiency.
* **Modular Design**: Perform only pruning or permutation testing depending on your analysis goals.
* **Flexible Configuration**: Adjustable parameters such as pruning depth, top N columns to prune, and permutation count.

## Installation
Clone the repository and install dependencies:

```bash
git clone https://github.com/username/bewitchi.git
cd bewitchi
pip install -r requirements.txt
```

## Usage
### Pruning Alignment
Prune alignment columns recursively based on Chi-square test:

```bash
python witchi.py prune --file alignment.fasta --format fasta --max_residue_pruned 100 --permutations 500 --num_workers 4 --top_n 3 --pruning_algorithm wasserstein
```
#### Options:
- `--file`: Path to the alignment file.
- `--format`: Alignment file format (default: fasta).
- `--max_residue_pruned`: Maximum columns to prune (default: 100).
- `--permutations`: Number of permutations for empirical distribution (default: 100).
- `--num_workers`: Number of CPU threads for parallelization (default: 2).
- `--top_n`: Number of top biased columns to prune per iteration (default: 1).
- `--pruning_algorithm`: Pruning algorithm to use (global, outlyingness, wasserstein, squared).

### Permutation Testing
Run permutation tests to establish empirical Chi-square distributions:

```bash
python witchi.py test --file alignment.fasta --format fasta --num_workers 4 --permutations 100
```

#### Options:
- `--file`: Path to the alignment file.
- `--format`: Alignment format (default: fasta).
- `--num_workers`: Number of CPU threads (default: 2).
- `--permutations`: Number of permutations (default: 100).

## Pruning Algorithms
- **Global Pruning**: Prioritizes columns with the highest global delta Chi-square score.
- **Outlyingness Pruning**: Focuses on removing columns contributing to outlier taxa.
- **Wasserstein Pruning**: Guides pruning by minimizing Wasserstein distance to the unbiased distribution.
- **Squared Pruning**: Targets columns that maximize squared Chi-square score differences.

## Output
- Pruned Alignment File: A new alignment file with reduced bias.
- Statistical Reports: Summary of the pruning process, including Chi-square distributions and convergence metrics.
- Future: Plots: Visual representation of the pruning impact on Chi-square distributions.

## How It Works
1. Read Alignment:
  * The alignment file is parsed and converted to a NumPy array.
2. Empirical Thresholds:
  * Permutation tests generate expected Chi-square score distributions.
3. Pruning Loop:
  * Iteratively removes the most biased columns based on the selected algorithm.
  * Monitors progress using metrics like Wasserstein distance.
4. Final Output:
  * Produces a pruned alignment and statistical summaries.

## Example Workflow
1. Run permutation test:
```bash
python witchi.py test --file example.fasta --permutations 100
```
2. Prune alignment with Wasserstein distance guidance:
```bash
python witchi.py prune --file example.fasta --max_residue_pruned 50 --pruning_algorithm wasserstein
```

## Future Features
Real-Time Visualization: Live plotting of distribution shifts during pruning.

## License
Bewitchi is licensed under the MIT License.

Contributing
We welcome contributions! Please fork the repository, create a branch, and submit a pull request.

Acknowledgments
