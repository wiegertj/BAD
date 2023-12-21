import statistics
import pandas as pd
from itertools import product
from probables import BloomFilter
from scipy.stats import skew, kurtosis
from Bio import SeqIO


class KmerComputer:

    def __init__(self, msa_filepath, tree_filepath, model_filepath, query_filepath, isAA):
        self.msa_filepath = msa_filepath
        self.tree_filepath = tree_filepath
        self.model_filepath = model_filepath
        self.query_filepath = query_filepath
        self.isAA = isAA
        self.k_mer_length = 15
        self.k_mer_max_gap = 0.3
        self.k_mer_max_n = 0.3
        self.fp_rate_bf = 0.01
        self.bloom_filters_MSA = self.get_bloom_filters_MSA()

    def filter_gapped_kmers(self, sequence) -> list:
        """
        Returns a list of k-mers for the given sequence considering a max_gap_percentage.
        Ambiguity code gets resolved on the fly by considering each possible k-mer.

                Parameters:
                        :param sequence:  DNA sequence

                Returns:
                        :return kmers: list of k-mers

        """
        sequence = sequence.upper()
        kmer_list = []
        nucleotide_ambiguity_code = {
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['G', 'C'],
            'W': ['A', 'T'],
            'K': ['G', 'T'],
            'M': ['A', 'C'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T']
        }

        if not self.isAA:

            for i in range(len(sequence) - int(self.k_mer_length) + 1):

                kmer = sequence[i:i + int(self.k_mer_length)]
                gap_count = kmer.count('-')
                if (kmer != ('N' * int(self.k_mer_length))) and (
                        kmer.count('N') / int(self.k_mer_length) <= self.k_mer_max_n):
                    if gap_count / int(self.k_mer_length) <= self.k_mer_max_gap:

                        ambiguous_positions = [i for i, char in enumerate(kmer) if char in nucleotide_ambiguity_code]

                        expanded_kmers = []
                        if ambiguous_positions:
                            combinations = product(
                                *(nucleotide_ambiguity_code[char] for char in kmer if
                                  char in nucleotide_ambiguity_code))
                            for combination in combinations:
                                expanded_kmer = list(kmer)
                                for position, nucleotide in zip(ambiguous_positions, combination):
                                    expanded_kmer[position] = nucleotide
                                expanded_kmers.append(''.join(expanded_kmer))
                            kmer_list.extend(expanded_kmers)
                        else:
                            kmer_list.append(kmer)
        else:
            for i in range(len(sequence) - int(self.k_mer_length) + 1):
                kmer = sequence[i:i + int(self.k_mer_length)]
                gap_count = kmer.count('-')
                if gap_count / self.k_mer_length <= self.k_mer_max_gap:
                    kmer_list.append(kmer)

        return kmer_list

    def compute_string_similarity_statistics(self, query) -> tuple:
        """
        Computes k-mer similarity using a bloom filter of the query and all the bloom filters of the MSA sequences.
        Then it computes summary statistics over those k-mer similarities.
        Ambiguity code gets resolved on the fly by considering each possible k-mer for DNA data.

                Parameters:
                        :param query: DNA sequence
                        :param k: k-mer length
                        :param max_gap_percent: maximum percentage of a k-mer to be valid

                Returns:
                         :return tuple: (dataset, sampleId, mean_similarity, std_similarity, sk_similarity, kur_similarity)
        """
        kmers_query = self.filter_gapped_kmers(str(query.seq))
        query_bf = self.bloom_filter(set(kmers_query), len(kmers_query), self.fp_rate_bf)

        result_similarities = []
        for bloom_filter_ref in self.bloom_filters_MSA:

            if not (bloom_filter_ref[0] == query.id):
                similarity = 0
                for kmer in set(kmers_query):
                    similarity += bloom_filter_ref[1].check(kmer) * query_bf.check(kmer)

                result_similarities.append(
                    similarity / len(set(kmers_query)))

        mean_similarity = sum(result_similarities) / len(result_similarities)
        std_similarity = statistics.stdev(result_similarities)
        kur_similarity = kurtosis(result_similarities, fisher=True)
        sk_similarity = skew(result_similarities)

        return query.id,  mean_similarity, std_similarity, sk_similarity, kur_similarity

    def bloom_filter(self, filtered_kmers, size, fp_rate) -> BloomFilter:
        """
        Returns a bloomfilter with the k-mers added

                Parameters:
                        :param filtered_kmers: list of k-mers to be put into a bloom filter
                        :param size: size of the bloom filter
                        :param fp_rate: false positive rate of the bloom filter

                Returns:
                        :return bf: bloom filter with k-mers
        """

        bf_ = BloomFilter(size, fp_rate)

        for kmer in filtered_kmers:
            bf_.add(kmer)
        return bf_

    def get_bloom_filters_MSA(self) -> list:
        """
        Creates a list of Bloom filters for each sequence in the MSA
        Returns
            :return bloom_filters_MSA: list of Bloom filters
        """
        bloom_filters_MSA = []
        for record in SeqIO.parse(self.msa_filepath, 'fasta'):
            kmers = self.filter_gapped_kmers(str(record.seq))
            kmers = set(kmers)
            bf = self.bloom_filter(kmers, len(kmers), self.fp_rate_bf)
            bloom_filters_MSA.append((record.id, bf))
        return bloom_filters_MSA

    def compute_kmer_similarity(self, logger) -> pd.DataFrame:
        results = []
        counter = 0
        for query_record in SeqIO.parse(self.query_filepath, 'fasta'):
            counter += 1
            results.append(self.compute_string_similarity_statistics(query_record))
            if counter != 0 and counter % 20 == 0:
                logger.info(f"Finished computing k-mer features for {counter} queries ... ")

        return pd.DataFrame(results, columns=["queryId", "mean_15mer_similarity", "std_15mer_similarity",
                                              "skewness_15mer_similarity", "kurtosis_15mer_similarity"])
