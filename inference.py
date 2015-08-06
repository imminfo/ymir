from pyymir import *


if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-a", "--algorithm", help = "name of the algorithm for the statistical inference (default is 'em'). For a list of available algorithms with their parameters in this Ymir distribution run $python3 pyymir.py -a", default = "em", type = str)
    ap.add_argument("-o", "--output", help = "path to the folder for output model folders (default is './ymir_models/')", type = str, default = "./ymir_models/")
    args = ap.parse_args()

    files, input_check = parse_input(args)
    model, model_check = parse_model(args)
    format, format_check = parse_format(args)
    out_models, out_check = parse_output_models(files, args)

    if model_check and input_check and format_check and out_check:
        files, convert_flag = convert_files(files, format)
        if convert_flag:
            for inp, out in zip(files, out_models):
                print(inp, out)
                os.system("./build/inference $1 $2 ...")
    else:
        print("Can't process further, too many errors for me! T_T")