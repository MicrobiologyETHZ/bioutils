from Bio import SeqIO

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

