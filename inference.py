

from pyymir import *

if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-a", "--algorithm", help = "name of the algorithm for the statistical inference. For a list of available algorithms with their parameters in this Ymir distribution run $python3 pyymir.py -a", default = "em", type = str)
    ap.parse_args()

os.system("./build/inference $1 $2 ...")