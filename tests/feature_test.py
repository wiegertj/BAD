import os
import unittest
import pandas as pd
import shutil

from features.feature_extractor import FeatureExtractor

test_msa_path = os.path.abspath(os.path.curdir + "/data/raw/test_19473_0_reference.fasta")
test_model_path = os.path.abspath(os.path.curdir + "/data/raw/test_19473_0_msa_model.txt")
test_tree_path = os.path.abspath(os.path.curdir + "/data/raw/test_19473_0_taxon5.newick")
test_query_path = os.path.abspath(os.path.curdir + "/data/raw/test_19473_0_query.fasta")


class FeatureExtractorTest(unittest.TestCase):
    def setUp(self) -> None:
        os.chdir(os.path.dirname(__file__))
        tmp_folder_path = os.path.abspath(os.path.join(os.curdir, "test"))
        os.makedirs(tmp_folder_path)
        tmp_folder_path_work = os.path.abspath(os.path.join(os.curdir, "test", "tmp"))
        os.makedirs(tmp_folder_path_work)
        os.chdir(tmp_folder_path_work)
        self.feature_computer = FeatureExtractor(test_msa_path,
                                                 test_tree_path,
                                                 test_model_path, test_query_path,
                                                 "test", "raxml-ng", True, os.path.join(os.curdir, "bad.log"), threads='auto')
        self.features = self.feature_computer.extract_features()
        self.features_ground_truth = pd.read_csv(os.path.abspath(os.path.dirname(__file__) + "/data/comparison/ground_truth.csv"))

    def tearDown(self) -> None:
        shutil.rmtree(os.path.abspath(os.path.dirname(__file__) + "/test"))

    def test_features(self):
        self.assertEqual(round(self.features["mean_15mer_similarity"].iloc[0], 3),
                         round(self.features_ground_truth["mean_15mer_similarity"].iloc[0], 3))
        self.assertEqual(round(self.features["std_15mer_similarity"].iloc[0], 3),
                         round(self.features_ground_truth["std_15mer_similarity"].iloc[0], 3))
        self.assertEqual(round(self.features["skewness_15mer_similarity"].iloc[0],3),
                         round(self.features_ground_truth["skewness_15mer_similarity"].iloc[0], 3))
        self.assertEqual(round(self.features["kurtosis_15mer_similarity"].iloc[0], 3),
                         round(self.features_ground_truth["kurtosis_15mer_similarity"].iloc[0], 3))

        self.assertEqual(round(self.features["inv_site_matches_query_msa_t9"].iloc[0], 3),
                         round(self.features_ground_truth["inv_site_matches_query_msa_t9"].iloc[0], 3))
        self.assertEqual(round(self.features["transversion_frac_query_msa_t7"].iloc[0], 3),
                         round(self.features_ground_truth["transversion_frac_query_msa_t7"].iloc[0], 3))
        self.assertEqual(round(self.features["inv_site_std_frac_query_msa_t7"].iloc[0], 3),
                         round(self.features_ground_truth["inv_site_std_frac_query_msa_t7"].iloc[0], 3))
        self.assertEqual(round(self.features["transversion_frac_query_msa_t5"].iloc[0], 3),
                         round(self.features_ground_truth["transversion_frac_query_msa_t5"].iloc[0], 3))
        self.assertEqual(round(self.features["min_frac_query_msa_t5"].iloc[0], 3),
                         round(self.features_ground_truth["min_frac_query_msa_t5"].iloc[0], 3))

        self.assertEqual(round(self.features["std_branch_length"].iloc[0], 3),
                         round(self.features_ground_truth["std_branch_length"].iloc[0], 3))
        self.assertEqual(round(self.features["skewness_closeness_centrality"].iloc[0], 3),
                         round(self.features_ground_truth["skewness_closeness_centrality"].iloc[0], 3))

        self.assertEqual(round(self.features["kurtosis_25mer_similarity_perc_hash"].iloc[0], 3),
                         round(self.features_ground_truth["kurtosis_25mer_similarity_perc_hash"].iloc[0], 3))
        self.assertEqual(round(self.features["max_parsimony_subst_freq"].iloc[0], 3),
                         round(self.features_ground_truth["max_parsimony_subst_freq"].iloc[0], 3))

        self.assertAlmostEqual(self.features["mean_nrf_parsimony_trees"].iloc[0],
                               self.features_ground_truth["mean_nrf_parsimony_trees"].iloc[0], delta=0.1)
        self.assertAlmostEqual(self.features["no_topologies_parsimony_bootstrap"].iloc[0],
                               self.features_ground_truth["no_topologies_parsimony_bootstrap"].iloc[0], delta=20)


if __name__ == '__main__':
    unittest.main()
