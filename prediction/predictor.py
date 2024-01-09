import os
import sys
import time
import shap
import pandas as pd
import matplotlib.pyplot as plt
from features.feature_extractor import FeatureExtractor
from models.get_model import get_model
from utils.utils import setup_logger, workdir_is_clean, check_file_exists


class Predictor:
    """
        Class for making predictions and writing the results.

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
        query_filepath : str
            absolute path to the .fasta query file
        current_directory : os.path
            current working directory
        output_prefix : str
             name of the ouput BAD produces
        prediction_model : LightGBM.Booster
            model predicting the median of the bootstrap
        prediction : pd.DataFrame
            dataframe with all prediction results

    """

    def __init__(self, msa_filepath, tree_filepath, model_filepath, query_filepath, o="BAD_output",
                 raxml_ng_path="raxml-ng", redo=False, shapley_calc=True, threads='auto'):
        self.logger = setup_logger("Predictor", os.path.join(os.pardir, "bad.log"))
        check_file_exists(msa_filepath, "Fasta MSA file", self.logger)
        check_file_exists(tree_filepath, "Newick tree file", self.logger)
        check_file_exists(model_filepath, "RAxML-NG model file", self.logger)
        self.current_directory = os.path.abspath(os.curdir)
        self.tree_filepath = tree_filepath
        self.prediction_model = get_model("bad")
        self.redo = redo
        self.output_prefix = o
        self.prediction = None
        self.shapley_calc = shapley_calc

        tmp_folder_path = os.path.abspath(os.path.join(os.curdir, self.output_prefix))
        if not (workdir_is_clean(tmp_folder_path, self.redo, o)):
            self.logger.error(
                f"Found exisiting result folder: {tmp_folder_path}, please rename the folder or run in -redo mode. Exiting BAD.")
            sys.exit()

        self.logger = setup_logger("Predictor", os.path.join(tmp_folder_path, "bad.log"))
        self.feature_extractor = FeatureExtractor(msa_filepath,
                                                  tree_filepath,
                                                  model_filepath, query_filepath,
                                                  o, raxml_ng_path, redo, os.path.join(tmp_folder_path, "bad.log"), threads)


    def print_result(self) -> None:
        self.logger.info("Feature values stored to: " + os.path.abspath(
            os.path.join(os.pardir, f"{self.output_prefix}_features.csv")))
        self.logger.info(
            "Difficulties stored to: " + os.path.abspath(os.path.join(os.pardir, f"{self.output_prefix}_result.csv")))
        if self.shapley_calc:
            self.logger.warning(f""" Be aware on how to interpret the results
                   If you want to interpret the Shapley value explanation, please make sure you understand how Shapley values are interpreted. 
                   Easy introduction to Shapley values: https://christophm.github.io/interpretable-ml-book/shapley.html
    
                   For further details on the prediction explanation, please have a look at the Shapley waterfall plots stored at: 
                   {os.path.abspath(os.path.join(os.pardir, f'{self.output_prefix}_QUERYID_shapley_plot.png'))}
               """)

    def predict(self) -> None:
        start_time = time.time()
        self.predict_shapley() if self.shapley_calc else self.predict_no_shapley()
        self.logger.info(f"Total elapsed time: {round(time.time() - start_time, 2)} seconds")

    def predict_no_shapley(self) -> None:
        """
        Function that performs the prediction without the shapley calculations
        """
        features = self.feature_extractor.extract_features()
        self.prediction = features
        self.prediction["placement_difficulty_prediction"] = self.prediction_model.predict(features.drop(columns=["queryId"]))
        features.to_csv(os.path.abspath(os.path.join(os.pardir, f"{self.output_prefix}_features.csv")))
        results = features[["queryId", "placement_difficulty_prediction"]]
        results.to_csv(os.path.abspath(os.path.join(os.pardir, f"{self.output_prefix}_result.csv")))
        for i in range(0, results.shape[0]):
            self.logger.info(f"Query: {results.iloc[i]['queryId']} | Difficulty:{round(results.iloc[i]['placement_difficulty_prediction'], 3) : .3f}")
        self.print_result()

    def predict_shapley(self) -> None:
        """
        Function that performs the prediction and shapley calculations
        """

        # Groups for shapley explanation summary
        group_tree_space = ["mean_nrf_parsimony_trees", "no_topologies_parsimony_bootstrap"]
        group_inv_sites = ["inv_site_std_frac_query_msa_t7", "transversion_frac_query_msa_t7",
                           "transversion_frac_query_msa_t5",
                           "min_frac_query_msa_t5", "inv_site_matches_query_msa_t9"]
        group_sim_qs_msa = ["kurtosis_15mer_similarity", "skewness_15mer_similarity", "std_15mer_similarity",
                            "mean_15mer_similarity", "kurtosis_25mer_similarity_perc_hash"]
        group_tree_msa = ["max_parsimony_subst_freq", "std_branch_length", "skewness_closeness_centrality"]

        features = self.feature_extractor.extract_features()
        self.prediction = self.prediction_model.predict(features.drop(columns=["queryId"]))

        explainer = shap.Explainer(self.prediction_model)
        self.prediction = []

        for queryId in features["queryId"].values.tolist():
            df_tmp = features[features["queryId"] == queryId]
            shap_values = explainer.shap_values(df_tmp.drop(columns=["queryId"]))
            diff_pred = self.prediction_model.predict(df_tmp.drop(columns=["queryId"]))
            self.prediction.append(diff_pred[0])

            shap_values_group_tree_space = shap_values[:, [feature in group_tree_space for feature in df_tmp.drop(
                columns=["queryId"]).columns]].sum(
                axis=1)
            shap_values_group_inv_sites = shap_values[:, [feature in group_inv_sites for feature in df_tmp.drop(
                columns=["queryId"]).columns]].sum(
                axis=1)
            shap_values_group_sim_qs_msa = shap_values[:, [feature in group_sim_qs_msa for feature in df_tmp.drop(
                columns=["queryId"]).columns]].sum(
                axis=1)
            shap_values_group_tree_msa = shap_values[:, [feature in group_tree_msa for feature in df_tmp.drop(
                columns=["queryId"]).columns]].sum(
                axis=1)

            def add_plus(float_value):
                if float_value >= 0:
                    return f"+{round(float_value, 3) : .3f}"
                else:
                    return f"-{abs(round(float_value, 3)) : .3f}"

            sum_shap = add_plus(round(shap_values_group_tree_space[0], 3) + round(shap_values_group_inv_sites[0], 3) + round(shap_values_group_sim_qs_msa[0], 3)
                                + round(shap_values_group_tree_msa[0], 3))

            self.logger.info(f"""Query: {queryId} | Difficulty: {round(diff_pred[0], 3) : .3f}\n 
        {"Explanation overview" : ^55}\n
        {'Feature group' : <30}{'Sum of Shapley values' : ^30} 
        {'-'*30 : <30}{'-'*30 : ^30} 
        {'Tree space features:' : <30}{add_plus(shap_values_group_tree_space[0]) : ^30} 
        {'Invariant site features:' : <30}{add_plus(shap_values_group_inv_sites[0]) : ^30} 
        {'Query/MSA similarity features:' : <30}{add_plus(shap_values_group_sim_qs_msa[0]) : ^30} 
        {'Tree and MSA features:' : <30}{add_plus(shap_values_group_tree_msa[0]) : ^30}
        {'-'*30 : <30}{'-'*30 : ^30} 
        {'All: ' : <30}{sum_shap : ^30}\n
        {f'{round(explainer.expected_value, 3)} (Base value) {sum_shap} (Shapley values) ={round(diff_pred[0], 3) : .3f}'}
    """)

            base_values = explainer.expected_value
            shap.plots.waterfall(
                shap.Explanation(
                    values=shap_values[0],
                    base_values=base_values,
                    data=df_tmp.drop(columns=["queryId"]).iloc[0]
                ),
                show=False
            )
            plt.tight_layout()
            plt.savefig(os.path.abspath(os.path.join(os.pardir, f"{self.output_prefix}_{queryId}_shapley_plot.png")))
            plt.close()

        features["placement_difficulty_prediction"] = self.prediction
        features.to_csv(os.path.abspath(os.path.join(os.pardir, f"{self.output_prefix}_features.csv")))
        results = features[["queryId", "placement_difficulty_prediction"]]
        results.to_csv(os.path.abspath(os.path.join(os.pardir, f"{self.output_prefix}_result.csv")))
        self.print_result()
