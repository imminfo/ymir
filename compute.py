from pyymir import *


if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-w", "--workmode", help = "Either add computed probabilities to the input file ('add'), or print them to the separated output file as a vector ('sep')", default = "add", type = str)
    ap.add_argument("-o", "--output", help = "output file, list of comma-separated output files or a name of an output folder for output files", type = str, default = "")
    args = ap.parse_args()

    files, input_check = parse_input(args)
    model, model_check = parse_model(args)
    format, format_check = parse_format(args)
    out_files, out_check = parse_output_files(args.input, args.output)

    if model_check and input_check and format_check and out_check:
        os.system("./build/compute $1 $2 ...")
    else:
        print("Can't process further, too many errors for me! T_T")