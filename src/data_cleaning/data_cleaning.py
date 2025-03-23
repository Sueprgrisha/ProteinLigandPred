import sys
import numpy as np

valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")


def is_valid_protein_sequences(sequences):

    # Remove whitespace and convert to uppercase for all sequences
    sequences = np.char.replace(sequences, " ", "")
    sequences = np.char.replace(sequences, "\n", "")
    sequences = np.char.upper(sequences)

    # Check for valid characters using vectorized functions
    valid_chars = np.vectorize(lambda seq: all(char in valid_amino_acids for char in seq))
    has_no_stop_codons = np.vectorize(lambda seq: '*' not in seq and 'X' not in seq)
    has_min_length = np.vectorize(lambda seq: len(seq) >= 5)

    # Combine conditions
    valid_mask = valid_chars(sequences) & has_no_stop_codons(sequences) & has_min_length(sequences)

    return valid_mask

def check_prot_seq(exp_val, prot_seq_array):
    mask_non_prot = is_valid_protein_sequences(prot_seq_array)
    non_prot_array = prot_seq_array[~mask_non_prot]
    return list(np.setdiff1d(non_prot_array, np.array(exp_val)))
