import argparse
import gzip
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
    ap.add_argument("-m", "--model", help = "either an alias of the one from available models in Ymir or a path to a folder with a model's .json file. For a list of available models in this Ymir distribution run $python3 pyymir.py models", type = str)
    return ap


def check_filepath(filepath):
    # check if there is exist a file in the given filepath
    pass


def parse_input(args_input, filelist):
    # parse args.input, check for each file and update the given list of filepaths
    pass


def parse_model(alias_or_path):
    # check if the given string is an alias of the one of available models
    # or is a path to a model.
    # If it's an alias, than return a path to this model.
    # If it's a path then check for a model.json file in this folder.
    pass


def parse_format(args_input):
    pass


def parse_output_files(args_input, args_output):
    pass


def parse_output_models(args_input, args_output):
    pass


def reformat(filepath, format):
    pass


def extract_info(arg, jsdata):
    def _pretty_dict(key, d, shift):
        if type(d) is list:
            for key, val in sorted(d):
                _pretty_dict(key, val, shift + 1)
        elif type(d) is dict:
            for key, val in sorted(d.items()):
                _pretty_dict(key, val, shift + 1)
        else:
            if key == "alias":
                print("  " * shift, '"', d, '" :\t', sep = '', end = '')
            elif key == "comment":
                print(d)
            else:
                print("  " * shift, key, ": ", d, sep = '')

    print(jsdata[arg]["title"])
    for val in jsdata[arg]["values"]:
        _pretty_dict(val, jsdata[arg]["values"][val], 1)
    print()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--version", help = "Display an information about the current Ymir distribution", action = "store_true", type = str)
    ap.add_argument("-a", "--algos", help = "List of available algorithms for the statistical inference with their parameters in this Ymir distribution", action = "store_true", type = str)
    ap.add_argument("-f", "--formats", help = "list of possible input formats in this Ymir distribution", action = "store_true", type = str)
    ap.add_argument("-m", "--models", help = "List of available models in this Ymir distribution", action = "store_true", type = str)
    ap.add_argument("-s", "--scripts", help = "List of available ready-to-use scripts", action = "store_true", type = str)
    args = ap.parse_args()

    jsdata = None
    with open("./.info.json") as f:
        jsdata = json.load(f)
        extract_info("ymir", jsdata)
        some_true = False
        for arg, val in vars(args).items():
            if val and arg != "version":
                extract_info(arg, jsdata)
                some_true = True
        if not some_true:
            for arg, val in vars(args).items():
                if arg != "version":
                    extract_info(arg, jsdata)