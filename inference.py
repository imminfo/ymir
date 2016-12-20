import os
import sys

from pyymir import *


DIRNAME = os.path.dirname(sys.argv[0])
DIRNAME = DIRNAME if DIRNAME else "./"


if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-a", "--algorithm", help = "name of the algorithm for the statistical inference (default is 'em'). For a list of available algorithms with their parameters in this Ymir distribution run $python3 pyymir.py -a", default = "em", type = str)
    ap.add_argument("-o", "--output", help = "path to the folder for output model folders (default is './ymir_models/')", type = str, default = "./ymir_models/")
    ap.add_argument("-p", "--parameters", help = "statistical inference algorithm parameters in form '<param>=<value>,<param>=<value>'", default = "niter=10,memory-safe=0", type = str)
    ap.add_argument("-w", "--workmode", help = "[NOT IMPLEMENTED YET] either build a vector of graphs once and recompute their parameters at each step ('-w time', faster execution time), or rebuild graphs at each step ('-w mem', slow execution time, but consumes almost no RAM)", default = "time", type = str)
    ap = add_prealign(ap)
    args = ap.parse_args()

    files, input_check = parse_input(args)
    model, model_check = parse_model(args)
    converter, format_check = parse_format(args)
    out_models, out_check = parse_output_models(files, args)
    algo_params = dict(map(lambda x: x.split("="), args.parameters.split(",")))

    print()
    if model_check and input_check and format_check and out_check:
        for i in range(len(files)):
            print(i + 1, ":")
            conv_file, convert_flag = convert(files[i], converter)
            if convert_flag:
                print()
                os.system(" ".join([DIRNAME + "/build/scripts/Inference", conv_file, model, out_models[i], args.algorithm, "niter", algo_params["niter"], "memory-safe", algo_params["memory-safe"]]))
    else:
        print("Can't process further, too many errors for me! T_T")
