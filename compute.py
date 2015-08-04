

from pyymir import *

if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument("-w", "--workmode", help = "Either add computed probabilities to the input file ('add'), or print them to the separated output file as a vector ('sep')", default = "add", type = str)
    ap.parse_args()

os.system("./build/compute $1 $2 ...")