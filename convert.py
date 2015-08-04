import json
import argparse
import gz


# columns name for the Ymir input format
OUTPUT_NUC = "nuc"
OUTPUT_AA = "aa"
OUTPUT_V = "v"
OUTPUT_D = "d"
OUTPUT_J = "j"
OUTPUT_VEND = "vend"
OUTPUT_DSTART = "dstart"
OUTPUT_DEND = "dend"
OUTPUT_JSTART = "jstart"
OUTPUT_COLUMN_SEP = "\t"
OUTPUT_GENE_SEP = ","


def convert(filename, format_json_path):
	outfile = filename + ".ymir"


if __name__ == "__main__":
    # "input" - name of input file "f1.txt" or files with colons "f1.txt;f2.txt;f3.txt" or a folder with files with equal format.
    # output will be "f1.txt.ymir"
	# if mode == "convert"
	# if mode == "merge" then merge Ymir output with the input file
	# "format" - name of .json file w/o ".json" from the parsers folder
	# "parser" - name of the .py file with a parser class
	# "filetype" - either human-readable tab-separated file (".ymirt") or binary file (".ymirb")
	ap = ArgumentParser("-i, --input", "-m, --mode", "-f, --format")
	a = 1 + 1;