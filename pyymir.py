import argparse
import json
import os
import sys

import converters


DIRNAME = os.path.dirname(sys.argv[0])
DIRNAME = DIRNAME if DIRNAME else "./"
INFO_JSON = DIRNAME + "/.info.json"


def make_workpath(script_run_path):
    if script_run_path.rfind('/') != -1:
        return script_run_path[:script_run_path.rfind('/') + 1]
    else:
        return './'


def default_ymir_ap():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", nargs = "+", help = "input file (text or gzipped) or a folder with input files (of the same format) or a list of space-separated files and/or folders in any combinations", type = str, default = "")
    ap.add_argument("-f", "--format", help = "format of input files (tcR, MiTCR, MiGEC, etc.) as an alias or as a Python 3 class from your module linked to the '$YMIR_HOME/converters' package. For a list of possible input formats and their aliases in this Ymir distribution run $python3 pyymir.py -f", type = str, default = "")
    ap.add_argument("-m", "--model", help = "either an alias of the one from available models in Ymir or a path to a folder with a model's .json file. For a list of available models in this Ymir distribution run $python3 pyymir.py -m", type = str, default = "")
    return ap


def parse_input(args):
    """
    Parse the input argument and return a list of files (possibly found in the input folders).
    """

    def _get_files_from_folder(folderpath):
        # return [folderpath + "/" + f for f in os.listdir(folderpath) if f[0] != "." and "ymir_in" not in f]
        return [folderpath + "/" + f for f in os.listdir(folderpath) if f[0] != "."]

    args_input = args.input

    res = []
    flag = True

    if len(args_input):
        for x in args_input:
            if os.path.exists(x):
                if os.path.isdir(x):
                    res.extend(_get_files_from_folder(x))
                else:
                    res.append(x)
            else:
                print('\tERROR: file "', x, '" not found )\':', sep = "")
                flag = False
    else:
        print("\tERROR: please specify at least one input file ]:")
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

    if args_model:
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
    else:
        print("\tERROR: please specify a generation model >:")
        flag = False

    return res, flag


def parse_format(args):
    """
    Parse format and return a class for convert one repertoire format to Ymir's format.

    !!! // TODO: Also checks if all ok with number of related pairs of genes-alignments !!!
    """

    args_format = args.format

    res = None
    flag = True

    if len(args_format):
        with open(INFO_JSON) as f:
            jsdata = json.load(f)["formats"]["values"]
            # print(jsdata)
            if args_format in jsdata:
                res = jsdata[args_format]["pyclass"]
            else:
                # if format is custom than search for it
                print("\tERROR: unknown format or python converter class o:")

            res = getattr(converters, res)
    else:
        print("\tERROR: please specify a format for input files :3")
        flag = False

    return res, flag


def parse_output_files(infiles, args, seq_type):
    """
    Create the output folder for output files.
    """

    args_out = args.output
    if not os.path.exists(args_out):
        os.makedirs(args_out)
    postfix = "nuc" if seq_type == "n" else "aa";
    return [args_out + "/" + x[x.rfind("/") + 1:] + ".ymir_probs_" + postfix + ".txt" for x in infiles], True


def parse_output_models(infiles, args):
    """
    Create the output folder and all inner folders for each model for input files.
    """

    args_out = args.output
    res = []
    for x in infiles:
        new_x = args_out + "/" + x[x.rfind("/") + 1:x.rfind(".")] + ".model/"
        if not os.path.exists(new_x):
            os.makedirs(new_x)
        res.append(new_x)
    return res, True


def convert(filepath, converter):
    out_file = ""
    flag = True
    if converter is converters.YmirConverter:
        out_file = filepath
        flag = True
    else:
        out_file = ''

        if filepath.endswith('.txt.gz'):
            out_file = filepath[:filepath.rfind('.txt.gz')] + ".ymir_in.txt"
        elif filepath.endswith('.gz'):
            out_file = filepath[:filepath.rfind('.gz')] + ".ymir_in.txt"
        else:
            out_file = filepath[:filepath.rfind(".")] + ".ymir_in" + filepath[filepath.rfind("."):]

        print("Converting '", filepath, "' to '", out_file, "'", "...", end = "\t", sep = "")
        c = converter()
        flag = c.convert(file_in=filepath, file_out=out_file)
        if flag: print("Done.")
    return out_file, flag


def extract_info(arg, jsdata):
    """
    Extract information from Ymir's general info .json file.
    """

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


def add_prealign(ap):
    ap.add_argument("--prealign", help = "[NOT IMPLEMENTED YET] if supplied then pre-align all gene segments for each input sequence", action = "store_true")
    ap.add_argument("--threshold", help = "[NOT IMPLEMENTED YET] if '--pre-align' is supplied then this is a minimal alignment score for each gene class to align in form '--threshold=<Vscore>:<Dscore>:<Jscore>' (default = 10:5:10)", type = str)
    ap.add_argument("--scoreV", help = "[NOT IMPLEMENTED YET] if '--pre-align' is supplied then this is a score for each operation in V alignment in form '--scoreV=match:mismatch:ins:del' (default = 1:-1:-3:-3)", type = str)
    ap.add_argument("--scoreD", help = "[NOT IMPLEMENTED YET] if '--pre-align' is supplied then this is a score for each operation in D alignment in form '--scoreD=match:mismatch:ins:del' (default = 1:-2:-4:-4)", type = str)
    ap.add_argument("--scoreJ", help = "[NOT IMPLEMENTED YET] if '--pre-align' is supplied then this is a score for each operation in J alignment in form '--scoreJ=match:mismatch:ins:del' (default = 1:-1:-3:-3)", type = str)
    return ap


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--version", help = "Display an information about the current Ymir distribution", action = "store_true")
    ap.add_argument("-a", "--algos", help = "List of available algorithms for the statistical inference with their parameters in this Ymir distribution", action = "store_true")
    ap.add_argument("-f", "--formats", help = "list of possible input formats in this Ymir distribution", action = "store_true")
    ap.add_argument("-m", "--models", help = "List of available models in this Ymir distribution", action = "store_true")
    ap.add_argument("-s", "--scripts", help = "List of available ready-to-use scripts", action = "store_true")
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
