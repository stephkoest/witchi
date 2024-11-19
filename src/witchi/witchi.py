import argparse
from alignment_pruner.alignment_pruner import AlignmentPruner
from alignment_pruner.permutation_test import PermutationTest

def main():
    parser = argparse.ArgumentParser(description='Recursively prune alignment based on chi-square test.')
    subparsers = parser.add_subparsers(dest='command', required=True)

    prune_parser = subparsers.add_parser('prune', help='Prune alignment based on chi-square test.')
    prune_parser.add_argument('--file', required=True, help='File containing the alignment')
    prune_parser.add_argument('--format', default='fasta', help='Alignment format, default is fasta')
    prune_parser.add_argument('--max_residue_pruned', type=int, default=100,
                              help='Maximum number of columns to prune')
    prune_parser.add_argument('--permutations', type=int, default=100, help='Number of alignment permutations.')
    prune_parser.add_argument('--num_workers', type=int, default=2,
                              help='Number of parallel workers (cores) to use, default is 2.\nScales poorly with more than 4 workers.')
    prune_parser.add_argument('--top_n', type=int, default=1,
                              help='Number of top delta-chi2 columns to remove. Default is 1.')
    prune_parser.add_argument('--pruning_algorithm', default='squared',
                              help='Pruning algorithm to use: global, outlyingness, wasserstein, squared')

    test_parser = subparsers.add_parser('test', help='Run permutation test on alignment.')
    test_parser.add_argument('--file', required=True, help='File containing the alignment')
    test_parser.add_argument('--format', default='fasta', help='Alignment format. (fasta)')
    test_parser.add_argument('--num_workers', type=int, default=2,
                             help='Number of parallel workers (cores) to use. (2)')
    test_parser.add_argument('--permutations', type=int, default=100, help='Number of alignment permutations. (100)')

    args = parser.parse_args()

    if args.command == 'prune':
        pruner = AlignmentPruner(file=args.file, format=args.format, max_residue_pruned=args.max_residue_pruned,
                                 permutations=args.permutations, num_workers=args.num_workers, top_n=args.top_n,
                                 pruning_algorithm=args.pruning_algorithm)
        pruner.run()
    elif args.command == 'test':
        tester = PermutationTest(num_workers=args.num_workers, permutations=args.permutations)
        tester.run_test(alignment_file=args.file, alignment_format=args.format)

if __name__ == '__main__':
    main()