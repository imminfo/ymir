class RepertoireConverter:
    def __init__(self):
        self.init_ymir_names()

        self.col_sep = "\t"
        self.gene_sep = "\t"
        self.column_i = {"nuc": -1, "aa": -1,
                         "v": -1, "d": -1, "j": -1,
                         "vend": -1, "dstart": -1, "dend": -1, "jstart": -1}


    def init_ymir_names(self):
        # column names
        self.ymir_nuc = "nuc"
        self.ymir_aa = "aa"
        self.ymir_vgene = "vgene"
        self.ymir_dgene = "dgene"
        self.ymir_jgene = "jgene"
        self.ymir_vend = "vend"
        self.ymir_dstart = "dstart"
        self.ymir_dend = "dend"
        self.ymir_jstart = "jstart"
        self.ymir_col_sep = "\t"
        self.ymir_gene_sep = ","


    def convert(self, file_in, file_out):
        """
        Parse and convert input file to the default Ymir's format file.

        :param file_in: path to the input file
        :param file_out: path to the output file
        :return:
        """
        skip = 0
        with open(file_in) as fin:
            skip = self.compute_skip(fin)

        with open(file_in) as fin:
            with open(file_out) as fout:
                for line in fin:
                    skip -= 1
                    if skip == 0: break

                self.parse_header(fin.readline())
                self.write_header(fout)
                for line in fin:
                    fout.write(self.parse_line(line) + "\n")


    def parse_line(self, line):
        raise NotImplementedError


    def parse_header(self, header):
        pass


    def write_header(self, fout):
        fout.write(self.ymir_col_sep.join([self.ymir_nuc, self.ymir_aa,
                                           self.ymir_vgene, self.ymir_dgene, self.ymir_jgene,
                                           self.ymir_vend, self.ymir_dstart, self.ymir_dend, self.ymir_jstart]) + "\n")


    def compute_skip(self, fin):
        return 0


class tcRConverter (RepertoireConverter):
    def __init__(self):
        self = RepertoireConverter.__init__()


    def parse_line(self):
        pass


class MiTCRConverter (RepertoireConverter):
    def __init__(self):
        self = RepertoireConverter.__init__()


    def parse_line(self):
        pass


    def compute_skip(self, fin):
        # check for levels of MiTCR output
        pass