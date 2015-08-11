from pyymir import *


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--model", help = "either an alias of the one from available models in Ymir or a path to a folder with a model's .json file. For a list of available models in this Ymir distribution run $python3 pyymir.py -m", type = str, default = "")
    ap.add_argument("-c", "--count", help = "number of clonotypes to generate", default = 10000, type = int)
    ap.add_argument("-o", "--output", help = "output file", type = str, default = "./out.ymir_gen.txt")
    args = ap.parse_args()

    model, model_check = parse_model(args)
    converter, format_check = parse_format(args)

    print()
    if model_check and input_check and format_check and out_check:
        for i in range(len(files)):
            print(i + 1, ":")
            conv_file, convert_flag = convert(files[i], converter)
            if convert_flag:
                print()
                os.system(" ".join(["./build/Compute", conv_file, model, out_files[i], "0" if args.predefined else "1"]))
            print()
    else:
        print("Can't process further, too many errors for me! T_T")