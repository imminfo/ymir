from pyymir import *


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--model", help = "either an alias of the one from available models in Ymir or a path to a folder with a model's .json file. For a list of available models in this Ymir distribution run $python3 pyymir.py -m", type = str, default = "")
    ap.add_argument("-c", "--count", help = "number of clonotypes to generate", default = 10000, type = int)
    ap.add_argument("-o", "--output", help = "output file", type = str, default = "./ymir_artif/out.ymir_gen.txt")
    ap.add_argument("--coding", help = "generate only coding sequences", action = "store_true")
    ap.add_argument("--noncoding", help = "generate only noncoding sequence", action = "store_true")
    ap.add_argument("--oof", help = "generate only out-of-frame sequences", action = "store_true")
    ap.add_argument("--stopcodon", help = "generate only sequences with stop-codon", action = "store_true")
    args = ap.parse_args()

    model, model_check = parse_model(args)
    if not os.path.exists(args.output[:args.output.rfind("/")]):
        os.makedirs(args.output[:args.output.rfind("/")])

    if args.model[-1] == '/':
        args.model = args.model[:-1]
    args.model = args.model[args.model.rfind('/') + 1:]

    out_file = args.output
    if out_file == "./ymir_artif/out.ymir_gen.txt":
        out_file = "./ymir_artif/out.ymir_gen." + args.model + ".txt"

    cod_nonc = '0'
    if args.coding:
        cod_nonc = '1'
    elif args.noncoding:
        cod_nonc = '2'
    elif args.oof:
        cod_nonc = '3'
    elif args.stopcodon:
        cod_nonc = '4'

    if model_check:
        print()
        os.system(" ".join(["./build/scripts/Generate", model, str(args.count), out_file, cod_nonc]))
        print()
    else:
        print("Can't process further, too many errors for me! T_T")
