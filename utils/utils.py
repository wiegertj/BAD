from Bio import SeqIO
import logging
import os
import subprocess
import sys
import shutil

def check_valid_characters(sequence) -> bool:
    """
    Checks whether all characters are valid DNA + ambiguity chars
    Parameters
        :param sequence: sequence to check
    Returns
        :return: True if DNA, False else
    """
    valid_characters = set("ATUCGMRWSYKVHDBN-_")
    if set(sequence.upper()) <= valid_characters:
        return True
    else:
        return False


def is_AA(msa_filepath):
    """
    Checks if MSA is AA data
    Parameters
        :param msa_filepath: MSA filepath to check
    Returns
        :return: True if AA, False else
    """
    with open(msa_filepath, "r") as fasta_file:
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        for seq_record in sequences:
            isDNA = check_valid_characters(str(seq_record.seq))
            if not isDNA:
                return True
        return False


def setup_logger(instance_name, log_file_path):
    logger = logging.getLogger(instance_name)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.INFO)  # Adjust the level as needed
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger


def check_file_exists(file_path, description, logger):
    if not os.path.exists(file_path):
        logger.error(
            f"MSA file not found at: {description}. Please check if the file exists and provide the absolute path.")
        sys.exit()
    return file_path


def check_raxml_availability():
    if os.name == "posix":
        try:
            return subprocess.check_output(["which", "raxml-ng"], text=True).strip()
        except subprocess.CalledProcessError:
            return None
    elif os.name == "nt":
        try:
            return subprocess.check_output(["where", "raxml-ng"], text=True).strip()
        except subprocess.CalledProcessError:
            return None
    return None


def workdir_is_clean(tmp_folder_path, redo, o):
    if os.path.exists(tmp_folder_path):
        if redo:
            shutil.rmtree(tmp_folder_path)
        else:
            return False

    os.makedirs(tmp_folder_path)
    tmp_folder_path_work = os.path.abspath(os.path.join(os.curdir, o, "tmp"))
    os.makedirs(tmp_folder_path_work)
    os.chdir(tmp_folder_path_work)
    return True


def normalize_branch_lengths(tree):
    """
    Normalizes branch lengths according to the sum of all branch lengths
    :param tree: ete3 tree to normalize the branch lengths of
    :return: ete3 tree with normalized branch lengths
    """
    total_length = 0.0

    for node in tree.traverse():
        if not node.is_root():
            total_length += node.dist

    for node in tree.traverse():
        if node.up:
            node.dist /= total_length
    return tree
