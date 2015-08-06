from pyymir import *


if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-w", "--workmode", help = "either add computed probabilities to the input file ('add', default), or print them to the separated output file as a vector ('sep')", default = "add", type = str)
    ap.add_argument("-o", "--output", help = "path to the output folder for output files (default is './ymir_genprob/')", type = str, default = "./ymir_genprob/")
    args = ap.parse_args()

    files, input_check = parse_input(args)
    model, model_check = parse_model(args)
    converter, format_check = parse_format(args)
    out_files, out_check = parse_output_files(files, args)

    if model_check and input_check and format_check and out_check:
        for i in range(len(files)):
            print(i + 1, ":")
            conv_file, convert_flag = convert(files[i], converter)
            if convert_flag:
                print(" ".join(["./build/compute", "-i", conv_file, "-m", model, "-o", out_files[i]]), sep = "")
    else:
        print("Can't process further, too many errors for me! T_T")