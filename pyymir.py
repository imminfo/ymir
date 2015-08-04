import argparse
import gz
import json
import os


# columns name for the Ymir input format
OUTPUT_NUC = "nuc"
OUTPUT_AA = "aa"
OUTPUT_V = "v"
OUTPUT_D = "d"
OUTPUT_J = "j"
OUTPUT_VEND = "vend"
OUTPUT_DSTART = "dstart"
OUTPUT_DEND = "dend"
OUTPUT_JSTART = "jstart"
OUTPUT_COLUMN_SEP = "\t"
OUTPUT_GENE_SEP = ","


def convert(filename, format_json_path):
	outfile = filename + ".ymir.txt"


def default_ymir_ap():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", help = "input file, list of comma-separated files or a name of a folder with input (text or gzipped text) files", type = str)
    ap.add_argument("-f", "--format", help = "format of input files: tcR, MiTCR, MiGEC, etc. For a list of possible input formats in this Ymir distribution run $python3 pyymir.py formats", type = str)
    ap.add_argument("-m", "--model", help = "path to a folder with a model's .json file. For a list of available models in this Ymir distribution run $python3 pyymir.py models", type = str)
    return ap


def extract_info(sid):
    pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("algos", help = "List of available algorithms for the statistical inference with their parameters in this Ymir distribution", action = "store_true")
    ap.add_argument("formats", help = "list of possible input formats in this Ymir distribution", action = "store_true")
    ap.add_argument("models", help = "List of available models in this Ymir distribution", action = "store_true")
    ap.parse_args()

    for arg in ap:
        extract_info(arg)