import os
import shutil
import subprocess
import logging
from scipy import stats
import fnmatch
import glob


def run_command(cmd, check=True):
    logging.info(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True, check=check)


def list_files_with_extension(folder, extension):
    return fnmatch.filter(os.listdir(folder), f"*.{extension}")


def file_exists(path):
    return os.path.lexists(path)


def fdr(p_vals):
    """FDR correction for multiple comparisons."""
    ranked_p_values = stats.rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr


def move_to_root_folder(root_path, cur_path):
    """Move downloaded NCBI genomes to the top folder level."""
    for filename in os.listdir(cur_path):
        if os.path.isfile(os.path.join(cur_path, filename)):
            shutil.move(os.path.join(cur_path, filename), os.path.join(root_path, filename))
        elif os.path.isdir(os.path.join(cur_path, filename)):
            move_to_root_folder(root_path, os.path.join(cur_path, filename))
        else:
            sys.exit("Should never reach here.")
    # remove empty folders
    if cur_path != root_path:
        os.rmdir(cur_path)


def gunzip_all(folder):
    """Decompress all .gz files in the specified folder."""
    gzip_pattern = os.path.join(folder, "*.gz")
    run_command(f"gunzip {gzip_pattern} || true")


def rename_file_given_extension(folder, extension):
    """Find the first file in the folder with given extension and rename it to reference.x"""
    all_files = glob.glob(os.path.join(folder, "*."+extension))
    first_file = all_files[0]  # Take the first file found
    os.rename(first_file, os.path.join(folder, "reference."+extension))  # Rename the file


def setup_logging(log_file):
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            # logging.StreamHandler(sys.stdout)
            logging.FileHandler(log_file, mode='a', encoding='utf-8')
        ]
    )


def check_dependencies(tools):
    missing_tools = []
    for tool in tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)
    return missing_tools

