import os
import sys

from pyymir import *


DIRNAME = os.path.dirname(sys.argv[0])
DIRNAME = DIRNAME if DIRNAME else "./"

if __name__ == "__main__":
    ap = default_ymir_ap()
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", nargs = "+", help = "input file (text or gzipped) or a folder with input files (of the same format) or a list of space-separated files and/or folders in any combinations", type = str, default = "")
    ap.add_argument("-f", "--format", help = "format of input files (tcR, MiTCR, MiGEC, etc.) as an alias or as a Python 3 class from your module linked to the '$YMIR_HOME/converters' package. For a list of possible input formats and their aliases in this Ymir distribution run $python3 pyymir.py -f", type = str, default = "")
    args = ap.parse_args()

    files, input_check = parse_input(args)
    converter, format_check = parse_format(args)

    print()
    if input_check and format_check:
        for i in range(len(files)):
            print(i + 1, " / ", len(files), ":")
            conv_file, convert_flag = convert(files[i], converter)
            if convert_flag:
                print("Successfully converted.")
            else:
                print("Convertation failed, skipping the file.")
    else:
        print("Can't process further, too many errors for me! x_x")
