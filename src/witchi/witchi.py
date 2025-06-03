import argparse
from witchi.alignment_pruner.alignment_pruner import AlignmentPruner
from witchi.alignment_pruner.permutation_test import PermutationTest


def main():
    parser = argparse.ArgumentParser(
        description="Recursively prune alignment based on chi-square test."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prune_parser = subparsers.add_parser(
        "prune", help="Prune alignment based on chi-square test."
    )
    prune_parser.add_argument(
        "--file", required=True, help="File containing the alignment"
    )
    prune_parser.add_argument(
        "--format", default="fasta", help="Alignment format, default is fasta"
    )
    prune_parser.add_argument(
        "--max_residue_pruned",
        type=int,
        default=100,
        help="Maximum number of columns to prune",
    )
    prune_parser.add_argument(
        "--permutations",
        type=int,
        default=100,
        help="Number of alignment permutations.",
    )
    prune_parser.add_argument(
        "--num_workers",
        type=int,
        default=2,
        help="Number of parallel workers (cores) to use, default is 2.\nScales poorly with more than 4 workers.",
    )
    prune_parser.add_argument(
        "--top_n",
        type=int,
        default=1,
        help="Number of top delta-chi2 columns to remove. Default is 1.",
    )
    prune_parser.add_argument(
        "--pruning_algorithm",
        default="quartic",
        help="Pruning algorithm to use: squared, quartic, wasserstein",
    )
    # add touchdown
    prune_parser.add_argument(
        "--touchdown",
        action="store_true",
        help="If flag is set, Touchdown mode is activated (experimental!). "
        "Consider setting top_n to around 1-2%% of the original alignment length.\n"
        "The number of "
        "columns to prune will be reduced to 0.1%% of the original alignment length "
        "per iteration after a threshold of 99.0%% is reached. ",
    )

    test_parser = subparsers.add_parser(
        "test", help="Run permutation test on alignment."
    )
    test_parser.add_argument(
        "--file", required=True, help="File containing the alignment"
    )
    test_parser.add_argument(
        "--format", default="fasta", help="Alignment format. (fasta)"
    )
    test_parser.add_argument(
        "--num_workers",
        type=int,
        default=2,
        help="Number of parallel workers (cores) to use. (2)",
    )
    test_parser.add_argument(
        "--permutations",
        type=int,
        default=100,
        help="Number of alignment permutations. (100)",
    )
    test_parser.add_argument(
        "--create_output",
        action="store_true",
        help="Flag to create output file with scores.",
    )

    args = parser.parse_args()

    if args.command == "prune":
        pruner = AlignmentPruner(
            file=args.file,
            format=args.format,
            max_residue_pruned=args.max_residue_pruned,
            permutations=args.permutations,
            num_workers=args.num_workers,
            top_n=args.top_n,
            pruning_algorithm=args.pruning_algorithm,
            touchdown=args.touchdown,
        )
        pruner.run()
    elif args.command == "test":
        tester = PermutationTest(
            num_workers=args.num_workers, permutations=args.permutations
        )
        tester.run_test(
            alignment_file=args.file,
            alignment_format=args.format,
            create_output=args.create_output,
        )


if __name__ == "__main__":
    main()
