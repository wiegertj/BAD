import argparse
from rich_argparse import RichHelpFormatter
from prediction.predictor import Predictor


def main():
    parser = argparse.ArgumentParser(description="Bold Assertor of Difficulty (BAD)", formatter_class=RichHelpFormatter)
    parser.add_argument("-msa", required=True, type=str, help=f"{'ABSOLUTE PATH |':<15}{'absolute path to the MSA file in fasta format':^50}")
    parser.add_argument("-tree", required=True, type=str,
                        help=f"{'ABSOLUTE PATH |':<15}{'absolute path to the tree file in newick format (.bestTree)':^50}")
    parser.add_argument("-model", required=True, type=str,
                        help=f"{'ABSOLUTE PATH |':<15}{'absolute path to the raxml-ng model file (.bestModel)':^50}")
    parser.add_argument("-query", required=True, type=str,
                        help=f"{'ABSOLUTE PATH |':<15}{'absolute path to the query file in fasta format':^50}")
    parser.add_argument("-o", required=False, type=str, help=f"{'VALUE |':<15}{'output folder name':^50}{'| default: BAD_output':>20}")
    parser.add_argument("-raxmlng", required=False, type=str, default="raxml-ng",
                        help=f"{'ABSOLUTE PATH |':<15}{'path to RAxML-NG':^50}{'| default: raxml-ng':>20}")
    parser.add_argument("-redo", action='store_true', required=False,
                        help=f"{'':<15}{'if set, all existing results will be overwritten':^50}")
    parser.add_argument("-shap", action='store_false', required=False,
                        help=f"{'':<15}{'if set, Shapley explanations of the prediction will be calculated':^50}")
    parser.add_argument("-threads", required=False, type=str, default='auto', help=f"{'NUMBER |':<15}{'number of threads to be used by RAxML-NG':^50}{'| default: auto':>20}")

    args = parser.parse_args()
    if args.o is None:
        output = "BAD_output"
    else:
        output = args.o
    predictor = Predictor(args.msa, args.tree, args.model, args.query, output, args.raxmlng, args.redo, args.shap, args.threads)
    predictor.predict()


if __name__ == "__main__":
    main()
