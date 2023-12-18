import numpy as np
from Bio import SeqIO
from scipy.fftpack import dct
from collections import defaultdict
from scipy.stats import kurtosis


def generate_k_mers(input_string, k) -> list:
    """
    Generates a list of k-mers
    :param input_string: string to create k-mers from
    :param k: length of the k-mers
    :return: list of k-mers
    """
    k_mers = []

    if len(input_string) < k:
        return k_mers

    for i in range(len(input_string) - k + 1):
        k_mer = input_string[i:i + k]
        k_mers.append(k_mer)

    return k_mers


def fraction_shared_kmers(sequence1, sequence2, k) -> float:
    """
    Calculates the fraction of shared k-mers
    Parameters
        :param binary_string1: first sequence
        :param binary_string2: second sequence
        :param k: length of the k-mers
    Returns
        :return: fraction of shared k-mers between both inputs
    """
    kmers1 = generate_k_mers(sequence1, k)
    kmers2 = generate_k_mers(sequence2, k)

    set_kmers1 = set(kmers1)
    set_kmers2 = set(kmers2)

    shared_kmers = set_kmers1.intersection(set_kmers2)

    try:

        fraction_shared = len(shared_kmers) / (len(set_kmers1) + len(set_kmers2) - len(shared_kmers))
    except ZeroDivisionError:
        fraction_shared = 0

    return fraction_shared


def dna_to_numeric(sequence, isAA) -> np.array:
    """
    Linearly encodes the sequence
    Parameters
        :param sequence: sequence to generate numeric representation of
        :param isAA: indicates if sequence is DNA (False) or AA (True)
    Returns
        :return: numerical encoded sequence as numpy array
    """
    if not isAA:
        mapping = {'A': 63, 'C': 127, 'G': 191, 'T': 255, '-': 0, 'N': 0} # linear mapping for DNA
        mapping = defaultdict(lambda: 0, mapping)
    else:
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        step = 256 // len(amino_acids)
        mapping = {aa: i * step for i, aa in enumerate(amino_acids)} # linear mapping for AA
        mapping = defaultdict(lambda: 0, mapping)
    numeric_sequence = [mapping[base] for base in sequence]
    return np.array(numeric_sequence)


def encode_dna_as_image(sequence) -> np.array:
    """
    Transforms sequence to square image
    Parameters
        :param sequence: sequence to transform in square image
    Returns
        :return: square image of sequence as numpy array
    """
    side_length = int(np.ceil(np.sqrt(len(sequence))))
    image = np.resize(sequence, (side_length, side_length))
    return image


def compute_hamming_distance(hash_value_1, hash_value_2) -> int:
    """
    Calculates the hamming distance between two hash values
    Parameters
        :param hash_value_1: hash of first sequence
        :param hash_value_2: hash of first sequence
    Returns
        :return: hamming distance
    """
    distance = sum(c1 != c2 for c1, c2 in zip(hash_value_1, hash_value_2))
    return distance


def compute_dct_sign_only_hash(sequence, isAA) -> str:
    """
    Computes the discrete cosine transformation sign only hash value of a sequence
    Parameters
        :param sequence: DNA or AA sequence to compute sign only hash of
        :param isAA: indicates of sequence isAA or not
    Returns
        :return: hash value of sequence
    """
    sequence = sequence.upper()
    numeric_sequence = dna_to_numeric(sequence, isAA)
    image = encode_dna_as_image(numeric_sequence)

    dct_coeffs = dct(dct(image, axis=0), axis=1)
    sign_only_sequence = np.sign(dct_coeffs)
    size_ = 16

    sign_only_sequence[sign_only_sequence >= 0] = 1
    sign_only_sequence[sign_only_sequence < 0] = 0
    binary_sequence = sign_only_sequence[:size_, :size_]

    hash_value = "".join([str(int(bit)) for bit in binary_sequence.flatten()])
    return hash_value


def compute_perceptual_kmer_similarity(msa_file, query_file, isAA) -> list:
    """
    Computes the perceptual hash 25-mer similarity between MSA and query sequences
    Parameters
        :param msa_file: MSA file for distance computation
        :param query_file: query file for distance computation
        :param isAA: defines if data is DNA or AA
    Returns
        :return: list with an entry for each query sequence and its kurtosis of the 25-mer similarities
    """
    results = []

    for record_query in SeqIO.parse(query_file, 'fasta'):
        kmer_sims25 = []
        hash_query = compute_dct_sign_only_hash(record_query.seq, isAA)
        for record_msa in SeqIO.parse(msa_file, 'fasta'):
            if record_msa.id != record_query.id:
                hash_msa = compute_dct_sign_only_hash(record_msa.seq, isAA)
                if hash_msa != 0:
                    kmer_sim25 = fraction_shared_kmers(hash_msa, hash_query, 25)
                    kmer_sims25.append(kmer_sim25)

        kur_kmer_sim25 = kurtosis(kmer_sims25, fisher=True)

        results.append((record_query.id, kur_kmer_sim25))
    return results
