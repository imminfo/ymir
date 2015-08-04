

from pyymir import *

if __name__ == "__main__":
    ap = default_ymir_ap()
    ap.add_argument()
    ap.parse_args()

os.system("./build/compute $1 $2 ...")