import argparse
from pathlib import Path
from Bio import SeqIO


def get_fasta_headers(sequence_file):
    if not Path(sequence_file).is_file():
        raise FileNotFoundError(f"Sequence file {sequence_file} does not exist.")
    
    with open(sequence_file, "r") as f:
        return [record.description for record in SeqIO.parse(f, "fasta")]

def check_embedding_exist(esm_embedding_dir, sequence_file):

    fasta_headers = get_fasta_headers(sequence_file)

    for protein_id in fasta_headers:
        embedding_file = Path(esm_embedding_dir) / f"{protein_id}.pt"
        if not embedding_file.exists():
            return False

    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--esm_embedding_dir', type=str, help='Output directory for the ESM embeddings')
    parser.add_argument('--sequence_file', type=str, help='Input Path to the sequence file of protein target')
    args = parser.parse_args()

    if check_embedding_exist(args.esm_embedding_dir, args.sequence_file):
        print("exists")
    else:
        print("not exists")

