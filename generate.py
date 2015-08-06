from pyymir import *


if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-c", "--count", help = "number of clonotypes to generate", default = 10000, type = int)
    ap.add_argument("-o", "--output", help = "output file", type = str, default = "out.ymir.txt")
    args = ap.parse_args()

    model, model_check = parse_model(args)

    if model_check:
        os.system("./build/generate $1 $2 ...")
    else:
        print("Can't process further, too many errors for me! T_T")