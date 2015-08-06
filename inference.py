

from pyymir import *

if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-a", "--algorithm", help = "name of the algorithm for the statistical inference. For a list of available algorithms with their parameters in this Ymir distribution run $python3 pyymir.py -a", default = "em", type = str)
    ap.add_argument("-o", "--output", help = "output folder, list of comma-separated output folders or a name of an output folder for output model folders", type = str, default = "")
    args = ap.parse_args()

    model, model_check = parse_model(args)
    files, input_check = parse_input(args)
    format, format_check = parse_format(args)
    out_models, out_check = parse_output_models(args.input, args.output)

    if model_check and input_check and format_check and out_check:
        pass

os.system("./build/inference $1 $2 ...")