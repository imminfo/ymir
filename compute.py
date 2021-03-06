import os
import sys

from pyymir import *


DIRNAME = os.path.dirname(sys.argv[0])
DIRNAME = DIRNAME if DIRNAME else "./"

if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-w", "--workmode", help = "either add computed probabilities to the input file ('add', default), or print them to the separated output file as a vector ('sep')", default = "add", type = str)
    ap.add_argument("-o", "--output", help = "path to the output folder for output files (default is './ymir_genprob/')", type = str, default = "./ymir_genprob/")
    ap.add_argument("-p", "--predefined", help = "recompute or use predefined model's gene usage (default is to recompute, to change it add '-p' to your script call)", action = "store_true")
    ap.add_argument("-l", "--leave", help = "add this to leave converted files (default is to remove converted files)", action = "store_true")
    ap.add_argument("-e", "--error", help = "error rate in the input data (default 0)", type = float, default = 0)
    ap.add_argument("-n", "--nuc", help = "compute nucleotide sequence probabilities (by default)", default = True, action = "store_true")
    ap.add_argument("-a", "--aa", help = "compute amino acid sequence probabilities", action = "store_true")
    ap = add_prealign(ap)
    args = ap.parse_args()

    files, input_check = parse_input(args)
    model, model_check = parse_model(args)
    converter, format_check = parse_format(args)
    out_files, out_check = parse_output_files(files, args, "n" if not args.aa else "a")

    print()
    if model_check and input_check and format_check and out_check:
        for i in range(len(files)):
            print(i + 1, " / ", len(files), ":")
            conv_file, convert_flag = convert(files[i], converter)
            if convert_flag:
                print()
                os.system(" ".join([DIRNAME + "/build/scripts/Compute", conv_file, model, out_files[i], "0" if args.predefined else "1", "n" if not args.aa else "a"]))

                if not args.leave and args.format != 'ymir':
                    os.remove(conv_file)
            print()
    else:
        print("Can't process further, too many errors for me! T_T")
