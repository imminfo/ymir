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

INFO_JSON = "./.info.json"


def default_ymir_ap():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", nargs = "+", help = "input file (text or gzipped) or a folder with input files (of the same format) or a list of space-separated files and/or folders in any combinations", type = str)
    ap.add_argument("-f", "--format", help = "format of input files: tcR, MiTCR, MiGEC, etc. For a list of possible input formats in this Ymir distribution run $python3 pyymir.py formats", type = str)
    ap.add_argument("-m", "--model", help = "either an alias of the one from available models in Ymir or a path to a folder with a model's .json file. For a list of available models in this Ymir distribution run $python3 pyymir.py models", type = str)
    return ap


def parse_input(args):
    """
    Parse the input argument and return a list of files (possibly found in the input folders).
    """

    args_input = args.input

    def _get_files_from_folder(folderpath):
        return [folderpath + "/" + f for f in os.listdir(folderpath) if f[0] != "."]

    res = []
    flag = True
    for x in args_input:
        if os.path.exists(x):
            if os.path.isdir(x):
                res.extend(_get_files_from_folder(x))
            else:
                res.append(x)
        else:
            print('\tERROR: file "', x, '" not found )\':', sep = "")
            flag = False

    return res, flag


def parse_model(args):
    """
    Check if the given string is an alias of the one of available models or
    it is a path to a folder with a model.
    If it's an alias, than return a path to this model.
    If it's a path then check for a model.json file in this folder.
    """

    args_model = args.model

    res = ""
    flag = True
    cur_models = {}

    with open(INFO_JSON) as f:
        jsdata = json.load(f)["models"]["values"]
        cur_models = {jsdata[key]["alias"] : jsdata[key]["path"] for key in jsdata}

    if args_model in cur_models:
        res = cur_models[args_model]
    else:
        if os.path.exists(args_model + "/model.json"):
            res = args_model
        else:
            print('\tERROR: model "', args_model, '" not found o:<', sep = "")

    return res, flag


def parse_format(args_input):
    return "", False


def parse_output_files(args_input, args_output):
    return "", False


def parse_output_models(args_input, args_output):
    return "", False


def convert(filepath, format):
    return filepath + ".ymir_in.txt"


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
    with open(INFO_JSON) as f:
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