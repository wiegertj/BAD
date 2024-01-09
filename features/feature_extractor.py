import time
import os
import pandas as pd
from features.feature_computer import FeatureComputer
from utils.utils import setup_logger


class FeatureExtractor:
    """
           Coordinates the feature computation.

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
           current_directory : os.path
               current working directory
           feature_computer : FeatureComputer
               object that encapsulates all feature computations

           Methods
           -------
           extract_features(self):
               triggers feature computation in the right order, merges result and tracks execution time
           """
    def __init__(self, msa_file_path, tree_file_path, model_file_path, query_file, output_prefix, raxml_ng_path, redo, log_file_path, threads):
        self.logger = setup_logger("FeatureExtractor", log_file_path)
        self.msa_file_path = msa_file_path
        self.tree_file_path = tree_file_path
        self.model_file_path = model_file_path
        self.feature_computer = FeatureComputer(msa_file_path, tree_file_path, model_file_path, query_file, output_prefix, raxml_ng_path, log_file_path, threads)
        self.current_directory = os.path.abspath(os.curdir)

    def extract_features(self) -> pd.DataFrame:
        """
        Triggers feature computation in the right order, merges result and tracks execution time

                Returns:
                        :return pd.DataFrame: dataframe with all features

        """
        self.logger.info("Starting feature extraction ... ")
        start_time = time.time()
        mean_nrf_parsimony_trees = self.feature_computer.compute_parsimony_features()
        no_topologies_parsimony_bootstrap = self.feature_computer.compute_parsimony_bootstrap_features()
        std_branch_length, skw_clo_sim = self.feature_computer.compute_tree_features()
        site_comp_features = self.feature_computer.compute_site_composition_features()
        df_site_comp_features = pd.DataFrame(site_comp_features, columns=["queryId", "transversion_frac_query_msa_t5",
                                                                          "transversion_frac_query_msa_t7",
                                                                          "inv_site_std_frac_query_msa_t7",
                                                                          "inv_site_matches_query_msa_t9",
                                                                          "min_frac_query_msa_t5"])
        pars_sub_features = self.feature_computer.compute_pars_subs_features()



        perc_hash_features = self.feature_computer.compute_perceptual_hash_features()
        df_perc_hash_features = pd.DataFrame(perc_hash_features,
                                             columns=["queryId", "kurtosis_25mer_similarity_perc_hash"])
        df_kmer_features = self.feature_computer.compute_kmer_similarity_features()

        df_merged = df_site_comp_features.merge(df_perc_hash_features, on=["queryId"], how="inner")
        df_merged = df_merged.merge(df_kmer_features, on=["queryId"], how="inner")
        df_merged["std_branch_length"] = std_branch_length
        df_merged["skewness_closeness_centrality"] = skw_clo_sim
        df_merged["mean_nrf_parsimony_trees"] = mean_nrf_parsimony_trees
        df_merged["no_topologies_parsimony_bootstrap"] = no_topologies_parsimony_bootstrap
        df_merged["max_parsimony_subst_freq"] = pars_sub_features

        df_merged["kurtosis_25mer_similarity_perc_hash"].fillna(-1, inplace=True)

        df_merged = df_merged[["queryId", "mean_15mer_similarity",
                               "std_15mer_similarity",	"skewness_15mer_similarity",	"kurtosis_15mer_similarity",	"inv_site_matches_query_msa_t9",	"transversion_frac_query_msa_t7",	"inv_site_std_frac_query_msa_t7",	"transversion_frac_query_msa_t5",	"min_frac_query_msa_t5",	"std_branch_length",	"skewness_closeness_centrality",	"kurtosis_25mer_similarity_perc_hash",	"max_parsimony_subst_freq",	"mean_nrf_parsimony_trees",	"no_topologies_parsimony_bootstrap"]]

        self.logger.info("Feature extraction finished!")
        elapsed_time = time.time() - start_time
        self.logger.info(f"Elapsed time: {round(elapsed_time, 2)} seconds")

        return df_merged
