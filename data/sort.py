from Bio import SeqIO
import re

# This script takes a fasta file and sorts the records by the numeric ID in the header.
# This was needed because seqkit was not sorting the fasta file correctly, i.e. record 2
# would be placed after all records starting with 1 such as 10, 100 etc.

# The script filters out "uncultured bacterium" records.

# Input and output file paths
input_file = "muscle.fasta"
output_file = "muscle-sorted.fasta"

def extract_id(record):
    """Extracts the numeric ID from the FASTA header (assumes ">id ...")."""
    match = re.match(r"^>(\d+)", record.description)
    if match:
        return int(match.group(1))
    # Fallback: try to extract from the first word
    try:
        return int(record.id)
    except ValueError:
        return float('inf')  # Put unparseable IDs at the end

# Read all records from the input FASTA file
records = list(SeqIO.parse(input_file, "fasta"))

# Sort records by extracted numeric ID
records.sort(key=extract_id)

# Write sorted records to the output FASTA file
SeqIO.write(records, output_file, "fasta")

print(f"Sorted {len(records)} records and saved to {output_file}")
