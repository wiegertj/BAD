import os
import re
import sys
import subprocess
import ete3
import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter
from scipy.stats import skew
from Bio import AlignIO, SeqRecord, Seq, SeqIO
from features.kmer_processing import KmerComputer
from features.perceptual_hash_processing import compute_perceptual_kmer_similarity
from features.substitution_processing import count_subst_freqs
from utils.utils import setup_logger, check_raxml_availability, normalize_branch_lengths, is_AA


class FeatureComputer:
    """
           Performs all feature computations necessary for BAD.

           Attributes
           ----------
           logger : logger
               for communcation with command line
           msa_filepath : str
               absolute path to the .fasta file
           model_filepath : str
               absolute path to the RAxML-NG model file
           tree_filepath : str
               absolute path to the .newick tree file
           output_prefix : str
                name of the ouput BAD produces
           raxml_ng_path : str
               path to the RAxML-NG executable

           Methods
           -------
           compute_parsimony_bootstrap_features(self):
               Computes the number of unique topologies of 200 parsimony bootstrap starting trees
           compute_parsimony_features(self):
               Computes mean of normalized RF distance for 1000 parsimony starting trees
           compute_kmer_similarity_features(self):
               Calls k-mer similarity computation
           compute_tree_features(self):
               Computes standard deviation of normalized branch lengths and skewness of closeness centrality
           compute_pars_subs_features(self):
               Computes parsimony substitution frequency features
           compute_perceptual_hash_features(self):
               Computes perceptual hash k-mer similarity features
           compute_site_composition_features(self):
               Computes features about the how the query matches invariant sites of the MSA

           """

    def __init__(self, msa_filepath, tree_filepath, model_filepath, query_filepath, output_prefix, raxml_ng_path, log_path, threads):
        self.msa_filepath = msa_filepath
        self.threads = threads
        self.tree_filepath = tree_filepath
        self.model_filepath = model_filepath
        self.query_filepath = query_filepath
        self.output_prefix = output_prefix
        self.raxml_ng_path = raxml_ng_path
        self.logger = setup_logger("FeatureComputer", log_path)
        self.isAA = is_AA(self.msa_filepath)
        if self.isAA:
            self.logger.info("Automatic data type recognition has recognized AA data")
        else:
            self.logger.info("Automatic data type recognition has recognized DNA data")
        if raxml_ng_path == "raxml-ng":
            self.raxml_path = check_raxml_availability()
            if self.raxml_path == None:
                self.logger.error(
                    "RAxML-NG not found in path variable. Please make sure it exists or pass the full path as input parameter (-raxmlng PATH)")
        else:
            self.raxml_path = raxml_ng_path
            if not os.path.isfile(self.raxml_path):
                self.logger.error(
                    f"RAxML-NG not found in the provided path {raxml_ng_path}. Please make sure it exists and provide the absolute path. Exiting BAD.")
                sys.exit()

    def compute_parsimony_bootstrap_features(self) -> float:

        """
            Performs the parsimony bootstrap. Samples alignment sites with replacement 200 times. For each samples MSA
            RAxML-NG is called to infer the corresponding parsimony starting tree. Stores those trees in a final file.
            Calculates the number of topologies of the set of trees.

                    Returns:
                            :return float: number of unique topologies

            """
        alignment = AlignIO.read(self.msa_filepath, "fasta")
        sequence_data = [list(record.seq) for record in alignment]
        alignment_array = np.array(sequence_data)
        original_ids = [record.id for record in alignment]

        trees_path = os.path.join(os.curdir, "parsimony_bootstraps_tmp.txt")

        for x in range(200):
            if (x % 40 == 0) and (x != 0):
                self.logger.info(f"Finished computing {x} from 200 parsimony bootstraps ... ")

            sampled_columns = np.random.choice(alignment_array.shape[1], size=alignment_array.shape[1],
                                               replace=True)

            replicate_alignment = alignment_array[:, sampled_columns]

            seq_records = [SeqRecord.SeqRecord(Seq.Seq(''.join(seq)), id=original_ids[i], description="") for i, seq
                           in
                           enumerate(replicate_alignment)]

            msa_new = AlignIO.MultipleSeqAlignment(seq_records)

            new_msa_path = os.path.join(os.curdir, "parsimony_bootstrap_tmp_" + str(x) + ".fasta")

            output_prefix = "parsimony_bootstrap_tmp_" + str(
                x)

            SeqIO.write(msa_new, new_msa_path, "fasta")

            raxml_command = [
                self.raxml_path,
                "--start",
                "--model", self.model_filepath,
                "--tree", "pars{1}",
                "--msa", new_msa_path,
                "--redo",
                "--prefix", output_prefix,
                "--threads", f"{self.threads}",
                "--log", "ERROR"
            ]
            try:
                subprocess.run(raxml_command, shell=False)
            except:
                self.logger.error()

            result_tree_path = os.path.join(os.curdir, output_prefix + ".raxml.startTree")

            try:
                with open(os.path.abspath(result_tree_path), 'r') as tree_file:
                    newick_tree = tree_file.read()
            except FileNotFoundError:
                self.logger.error(f"Starting trees not found: {result_tree_path}")
                continue

            with open(trees_path, 'a') as trees_file:
                if not os.path.exists(os.path.abspath(trees_path)):
                    trees_file.write(newick_tree)
                else:
                    trees_file.write(newick_tree)

        raxml_command = [self.raxml_path,
                         "--rfdist",
                         "--tree", os.path.abspath(trees_path),
                         "--redo",
                         "--prefix", output_prefix,
                         "--threads", "auto{60}"]
        result = subprocess.run(" ".join(raxml_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                shell=True)

        numbers = re.findall(r'set:\s+(-?[\d.]+)', result.stdout)
        numbers = [int(num) if num.isdigit() else float(num) for num in numbers]
        no_top_boot = numbers[2]
        self.logger.info("Finished computing parsimony bootstrap features!")
        return no_top_boot

    def compute_parsimony_features(self) -> float:
        """
        Calls RAxML-NG for creating 1000 parsimony starting trees, then it computes the normalized RF distance among them.

                Returns:
                        :return float: normalized RF distance

        """
        output_prefix = "parsimony_tmp_1000"

        raxml_command = [
            self.raxml_path,
            "--start",
            "--model", self.model_filepath,
            "--tree", "pars{1000}",
            "--msa", self.msa_filepath,
            "--redo",
            "--prefix", output_prefix,
            "--threads", f"{self.threads}",
            "--log", "ERROR",
        ]

        subprocess.run(raxml_command, shell=False)
        trees_path = os.path.join(os.curdir, "parsimony_tmp_1000.raxml.startTree")

        raxml_command = [self.raxml_path,
                         "--rfdist",
                         "--tree", os.path.abspath(trees_path),
                         "--redo",
                         "--prefix", output_prefix,
                         "--threads", f"{self.threads}"
                         ]
        result = subprocess.run(" ".join(raxml_command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                shell=True)
        numbers = re.findall(r'set:\s+(-?[\d.]+)', result.stdout)

        numbers = [int(num) if num.isdigit() else float(num) for num in numbers]
        avg_rel_rf_no_boot = numbers[1]
        self.logger.info("Finished computing parsimony starting tree features!")
        return avg_rel_rf_no_boot

    def compute_tree_features(self) -> (float, float):
        """
        Computes summary statistics for the tree.

                Returns:
                        :return (float, float): standard deviation of normalized branch lengths, skewness of closeness centrality
        """
        with open(self.tree_filepath, 'r') as file:
            newick_tree = file.read()
        tree = ete3.Tree(newick_tree)

        tree = normalize_branch_lengths(tree)
        branch_lengths = [node.dist for node in tree.traverse() if not node.is_root()]
        std_length = np.std(branch_lengths)

        G = nx.DiGraph()

        def traverse_and_add_edges(node):
            for child in node.children:
                edge_weight = node.get_distance(child)
                G.add_edge(node.name, child.name, weight=edge_weight)
                traverse_and_add_edges(child)

        traverse_and_add_edges(tree)
        closeness_centrality = nx.closeness_centrality(G, distance='weight')
        values = list(closeness_centrality.values())
        self.logger.info("Finished computing reference tree features!")

        return std_length, skew(values)

    def compute_site_composition_features(self) -> list:
        """
        Computes summary statistics of the matches of the query characters on invariant sites of the MSA.

                Returns:
                        :return list: list with an entry for all queries including:
                                      - queryId
                                      - transversion count on invariant sites (threshold: 0.5)
                                      - transversion count on invariant sites (threshold: 0.7)
                                      - standard deviation of query char fraction in invariant sites (threshold: 0.7)
                                      - fraction of matches between query and invariant sites (threshold: 0.9)
                                      - minimum of query char fraction in invariant sites (threshold: 0.5)

        """
        alignment = AlignIO.read(self.msa_filepath, 'fasta')
        analyzed_sites_9 = []
        analyzed_sites_8 = []
        analyzed_sites_7 = []
        analyzed_sites_5 = []

        # Calculate invariant sites with different thresholds
        for position in range(len(alignment[0])):
            residues_at_position = [str(record.seq[position]) for record in alignment]
            char_counts = Counter(residues_at_position)
            most_common_char, most_common_count = char_counts.most_common(1)[0]
            total_count = sum(count for char, count in char_counts.items())
            proportion_most_common = most_common_count / total_count if total_count > 0 else 0

            if proportion_most_common < 0.5 or most_common_char in ['-', 'N']:
                analyzed_sites_5.append((0, most_common_char))
            else:
                analyzed_sites_5.append((1, most_common_char))

            if proportion_most_common < 0.8 or most_common_char in ['-', 'N']:
                analyzed_sites_8.append((0, most_common_char))
            else:
                analyzed_sites_8.append((1, most_common_char))

            if proportion_most_common < 0.7 or most_common_char in ['-', 'N']:
                analyzed_sites_7.append((0, most_common_char))
            else:
                analyzed_sites_7.append((1, most_common_char))

            if proportion_most_common < 0.9 or most_common_char in ['-', 'N']:
                analyzed_sites_9.append((0, most_common_char))
            else:
                analyzed_sites_9.append((1, most_common_char))

        # Check how query is matching invariant sites
        results = []
        for record in SeqIO.parse(self.query_filepath, 'fasta'):

            match_counter_9 = 0
            total_inv_sites_9 = 0
            for i, (flag, char) in enumerate(analyzed_sites_9):
                # Check if the corresponding site in the query has a 1 and if the characters are equal
                if flag == 1:
                    total_inv_sites_9 += 1
                if flag == 1 and str(record.seq)[i] == char and char not in ['-', 'N']:
                    match_counter_9 += 1

            match_counter_5 = 0
            total_inv_sites_5 = 0
            for i, (flag, char) in enumerate(analyzed_sites_5):
                # Check if the corresponding site in the query has a 1 (is invariant w.r.t. t=0.5) and if the characters are equal
                if flag == 1:
                    total_inv_sites_5 += 1
                if flag == 1 and str(record.seq)[i] == char and char not in ['-', 'N']:
                    match_counter_5 += 1

            match_counter_7 = 0
            total_inv_sites_7 = 0
            for i, (flag, char) in enumerate(analyzed_sites_7):
                # Check if the corresponding site in the query has a 1 (is invariant w.r.t. t=0.7) and if the characters are equal
                if flag == 1:
                    total_inv_sites_7 += 1
                if flag == 1 and str(record.seq)[i] == char and char not in ['-', 'N']:
                    match_counter_7 += 1

            match_counter_8 = 0
            total_inv_sites_8 = 0
            for i, (flag, char) in enumerate(analyzed_sites_8):
                # Check if the corresponding site in the query has a 1 and if the characters are equal
                if flag == 1:
                    total_inv_sites_8 += 1
                if flag == 1 and str(record.seq)[i] == char and char not in ['-', 'N']:
                    match_counter_8 += 1

            # Get transitions/transversion fraction on all induced mutations
            transition_count7 = 0  # transitions induced by query on invariant sites with t=0.7
            transversion_count7 = 0  # transversions induced by query on invariant sites with t=0.7
            mut_count7 = 0  # all mutations induced by query on invariant sites with t=0.7
            fraction_char_rests7 = []  # fractions of the query residues in the invariant sites with t=0.7
            for i, (flag, char) in enumerate(analyzed_sites_7):
                # Check if the corresponding site in the query has a 1 and if the characters are equal
                if flag == 1:
                    total_inv_sites_7 += 1
                if flag == 1 and str(record.seq)[
                    i] != char:  # if site is invariant AND sequence does not match, check which mutation
                    mut_count7 += 1
                    if char in ["C", "T", "U"]:
                        if str(record.seq)[i] in ["A", "G"]:
                            transversion_count7 += 1
                    elif char in ["A", "G"]:
                        if str(record.seq)[i] in ["C", "T", "U"]:
                            transversion_count7 += 1
                    else:
                        transition_count7 += 1

                    residues_at_position = [str(record.seq[i]) for record in alignment]
                    residue_counts = Counter(residues_at_position)
                    most_common_residue, most_common_count = residue_counts.most_common(1)[0]
                    residues_at_position_del_most_common = [r for r in residues_at_position if r != most_common_residue]
                    if str(record.seq)[i] in residues_at_position_del_most_common:
                        count_char = residues_at_position_del_most_common.count(str(record.seq)[i])
                        fraction_char_rest = count_char / len(residues_at_position_del_most_common)
                    else:
                        fraction_char_rest = 0
                    fraction_char_rests7.append(fraction_char_rest)

            # Repeat for t=0.5
            transition_count5 = 0
            transversion_count5 = 0
            mut_count5 = 0
            fraction_char_rests5 = []
            for i, (flag, char) in enumerate(analyzed_sites_5):
                if flag == 1:
                    total_inv_sites_5 += 1
                if flag == 1 and str(record.seq)[i] != char:
                    mut_count5 += 1
                    if char in ["C", "T", "U"]:
                        if str(record.seq)[i] in ["A", "G"]:
                            transversion_count5 += 1
                    elif char in ["A", "G"]:
                        if str(record.seq)[i] in ["C", "T", "U"]:
                            transversion_count5 += 1
                    else:
                        transition_count5 += 1

                    residues_at_position = [str(record.seq[i]) for record in alignment]
                    residue_counts = Counter(residues_at_position)
                    most_common_residue, most_common_count = residue_counts.most_common(1)[0]
                    residues_at_position_del_most_common = [r for r in residues_at_position if r != most_common_residue]
                    if str(record.seq)[i] in residues_at_position_del_most_common:
                        count_char = residues_at_position_del_most_common.count(str(record.seq)[i])
                        fraction_char_rest = count_char / len(residues_at_position_del_most_common)
                    else:
                        fraction_char_rest = 0
                    fraction_char_rests5.append(fraction_char_rest)

            # Normalize results
            if mut_count7 > 0:
                transversion_count_rel7 = transversion_count7 / mut_count7
            else:
                transversion_count_rel7 = 0

            #########

            if mut_count5 > 0:
                transversion_count_rel5 = transversion_count5 / mut_count5
            else:
                transversion_count_rel5 = 0

            if len(fraction_char_rests7) > 0:
                std_fraction_char_rests7 = np.std(fraction_char_rests7)
            else:
                std_fraction_char_rests7 = 0.0

            if len(fraction_char_rests5) > 0:
                min_fraction_char_rests5 = np.min(fraction_char_rests5)
            else:
                min_fraction_char_rests5 = -1

            seq_length = len(str(record.seq))

            if self.isAA:
                transversion_count_rel5 = -1
                transversion_count_rel7 = -1

            results.append((record.id, transversion_count_rel5, transversion_count_rel7, std_fraction_char_rests7,
                            match_counter_9 / seq_length, min_fraction_char_rests5))
        self.logger.info("Finished computing invariant site features!")
        return results

    def compute_kmer_similarity_features(self) -> pd.DataFrame:
        kmer_computer = KmerComputer(self.msa_filepath, self.tree_filepath, self.model_filepath, self.query_filepath,
                                     self.isAA)
        kmer_similarities_df = kmer_computer.compute_kmer_similarity(self.logger)
        self.logger.info("Finished computing k-mer similarity features!")
        return kmer_similarities_df

    def compute_pars_subs_features(self) -> float:
        result = count_subst_freqs(self.tree_filepath, self.msa_filepath)
        with open(self.tree_filepath, "r") as tree_file:
            tree_data = tree_file.read()
            tree = ete3.Tree(tree_data)
            max_subst_freq = max(result) / len(tree.get_leaves())
        self.logger.info("Finished computing parsimony substitution features!")
        return max_subst_freq

    def compute_perceptual_hash_features(self) -> list:
        perc_hash = compute_perceptual_kmer_similarity(self.msa_filepath, self.query_filepath, self.isAA)
        self.logger.info("Finished perceptual hash features!")
        return perc_hash
