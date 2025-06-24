from Bio import SeqIO
from pathlib import Path
import pyranges as pr


def filter_pseudo_sequences(input_fasta, output_fasta):
    """
    Filters out sequences with 'pseudo=true' in their headers from a multi-FASTA file.

    Parameters
    ----------
    input_fasta : str
        The path to the input multi-FASTA file.
    output_fasta : str
        The path to the output multi-FASTA file.
    """
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        # Iterate over each sequence in the input multi-FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Check if 'pseudo=true' is not in the header
            if 'pseudo=true' not in record.description:
                # Write the sequence to the output file
                SeqIO.write(record, output_handle, "fasta")


def rename_protein_sequences(input_fasta, output_fasta, input_gff,
                             gff_col_replace='ID', gff_col_match='protein_id', sep=' ', match_term_index=2):
    gff = pr.read_gff3(input_gff)[[gff_col_replace, gff_col_match]].as_df()
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        # Iterate over each sequence in the input multi-FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            to_match = record.id.split(sep)[match_term_index]
            record.description = ''
            record.id = gff[gff[gff_col_match] ==
                            to_match][gff_col_replace].values[0]
            SeqIO.write(record, output_handle, "fasta")


if __name__ == "__main__":
    ifile = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/grab/dbcan_output/faa_files/GCA_001688725.2_ASM168872v2_protein.faa"
    ofile = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/grab/dbcan_output/faa_files/CP015401.2.edited.faa"
    gff = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/grab/annotations/CP015401.2.gff3"
    rename_protein_sequences(ifile, ofile, gff, sep=' ', match_term_index=0)
