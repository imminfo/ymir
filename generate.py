

from pyymir import *

if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-c", "--count", help = "number of clonotypes to generate", default = 10000, type = int)
    ap.parse_args()

os.system("./build/generate $1 $2 ...")