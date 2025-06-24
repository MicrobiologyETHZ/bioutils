
import argparse
from bioutils.sequence import filter_pseudo_sequences, rename_protein_sequences


def welcome():
    print("Welcome to bioinformatic util package")


def main():
    parser = argparse.ArgumentParser(
        description="Bioinformatics utility pacakge.")
    subparsers = parser.add_subparsers(
        title="Commands", dest="command", required=True)

    welcome_parser = subparsers.add_parser(
        "welcome", help="Printe welcome message")
    welcome_parser.set_defaults(func=lambda args: welcome())

    # Define the 'filter-pseudo' command
    nopseudo_parser = subparsers.add_parser(
        "filter-pseudo", help="Filter pseudo sequences from aa FASTA file")
    nopseudo_parser.add_argument("input_file", help="Path to the input file")
    nopseudo_parser.add_argument("output_file", help="Path to the output file")
    nopseudo_parser.set_defaults(
        func=lambda args: filter_pseudo_sequences(args.input_file, args.output_file))

    # Define the 'rename-headers' command
    newheaders_parser = subparsers.add_parser(
        "rename-headers", help="Rename headers in a FASTA file from a field in gff file")
    newheaders_parser.add_argument(
        "input_file", help="Path to the input fasta")
    newheaders_parser.add_argument(
        "output_file", help="Path to the output fasta file")
    newheaders_parser.add_argument("gff_file", help="Path to the gff file")
    newheaders_parser.add_argument(
        "gff_descriptor", help="ID to extract from gff", default='ID', nargs="?")
    newheaders_parser.add_argument(
        "match_in_gff", help="Field in gff file that matches the current header", default='protein_id', nargs="?")
    newheaders_parser.add_argument(
        "header_sep", help="Split header by", default=' ', nargs="?")
    newheaders_parser.add_argument(
        "match_in_header", help="Index of split header that would match a field in gff file", default=2, type=int, nargs="?")
    newheaders_parser.set_defaults(
        func=lambda args: rename_protein_sequences(args.input_file, args.output_file, args.gff_file,
                                                   args.gff_descriptor, args.match_in_gff, args.header_sep, args.match_in_header))

    # Parse the arguments
    args = parser.parse_args()

    # Call the appropriate function for the chosen command
    if hasattr(args, "func"):
        args.func(args)


if __name__ == "__main__":
    main()
