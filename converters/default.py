class RepertoireConverter:
    """
    Parent class for all other converters. Just rewrite functions to make them do what you want to.
    """

    def __init__(self):
        self.ymir_columns = {}  # dict <ymir col name> : <index of related column in the input file>
        self.init_ymir_names()

        self.input_columns = {}  # dict <input col name> : <ymir col name>
        self.init_input_columns()


    def init_ymir_names(self):
        """
        Initialise Ymir's format column names and separator characters.
        """

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

        self.ymir_columns = {self.ymir_nuc: -1, self.ymir_aa: -1,
                         self.ymir_v: -1, self.ymir_d: -1, self.ymir_j: -1,
                         self.ymir_vend: -1, self.ymir_dstart: -1, self.ymir_dend: -1, self.ymir_jstart: -1}


    def init_input_columns(self):
        """
        Initialise input files' column names and separator characters.
        """

        self.input_nuc = ""
        self.input_aa = ""
        self.input_vgene = ""
        self.input_dgene = ""
        self.input_jgene = ""
        self.input_vend = ""
        self.input_dstart = ""
        self.input_dend = ""
        self.input_jstart = ""
        self.input_col_sep = ""
        self.input_gene_sep = ""
        
        self.input_columns = {self.input_nuc: self.ymir_nuc, self.input_aa: self.ymir_aa,
                              self.input_v: self.ymir_v, self.input_d: self.ymir_d, self.input_j: self.ymir_j,
                              self.input_vend: self.ymir_vend, self.input_dstart: self.ymir_dstart, self.input_dend: self.ymir_dend, self.input_jstart: self.ymir_jstart}


    def convert(self, file_in, file_out):
        """
        Parse and convert input file to the default Ymir's format file.

        :param file_in: path to the input file
        :param file_out: path to the output file
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
        """
        Parse the input line and return new line in the Ymir format to write
        to the output file.

        :param line: a line from the input file
        :return: line in Ymir's format
        """

        words = line.strip().split(self.input_col_sep)
        out_words = []

        for key, val in self.ymir_columns.items():
            pass

        return self.ymir_col_sep.join(out_words)


    def parse_header(self, header):
        words = header.strip().split(self.col_sep)
        for i, w in enumerate(words):
            if w in self.input_columns:
                self.ymir_column_i[self.input_columns[w]] = i


    def write_header(self, fout):
        fout.write(self.ymir_col_sep.join([self.ymir_nuc, self.ymir_aa,
                                           self.ymir_vgene, self.ymir_dgene, self.ymir_jgene,
                                           self.ymir_vend, self.ymir_dstart, self.ymir_dend, self.ymir_jstart]) + "\n")


    def compute_skip(self, fin):
        return 0


class tcRConverter (RepertoireConverter):
    def __init__(self):
        self = RepertoireConverter.__init__()


    def init_input_columns(self):
        """
        Initialise input files' column names and separator characters.
        """

        self.input_nuc = "CDR3.nucleotide.sequence"
        self.input_aa = "CDR3.amino.acid.sequence"
        self.input_vgene = "V.gene"
        self.input_dgene = "D.gene"
        self.input_jgene = "J.gene"
        self.input_vend = "V.end"
        self.input_dstart = "D5.end"
        self.input_dend = "D3.end"
        self.input_jstart = "J.start"
        self.input_col_sep = "\t"
        self.input_gene_sep = ", "

        self.input_columns = {self.input_nuc: self.ymir_nuc, self.input_aa: self.ymir_aa,
                              self.input_v: self.ymir_v, self.input_d: self.ymir_d, self.input_j: self.ymir_j,
                              self.input_vend: self.ymir_vend, self.input_dstart: self.ymir_dstart, self.input_dend: self.ymir_dend, self.input_jstart: self.ymir_jstart}


class MiTCRConverter (RepertoireConverter):
    def __init__(self):
        self = RepertoireConverter.__init__()


    def init_input_columns(self):
        """
        Initialise input files' column names and separator characters.
        """

        self.input_nuc = "CDR3 nucleotide sequence"
        self.input_aa = "CDR3 amino acid sequence"
        self.input_vgene = "V segments"
        self.input_dgene = "D segments"
        self.input_jgene = "J segments"
        self.input_vend = "V last nucleotide position"
        self.input_dstart = "D first nucleotide position"
        self.input_dend = "D last nucleotide position"
        self.input_jstart = "J first nucleotide position"
        self.input_col_sep = "\t"
        self.input_gene_sep = ", "

        self.input_columns = {self.input_nuc: self.ymir_nuc, self.input_aa: self.ymir_aa,
                              self.input_v: self.ymir_v, self.input_d: self.ymir_d, self.input_j: self.ymir_j,
                              self.input_vend: self.ymir_vend, self.input_dstart: self.ymir_dstart, self.input_dend: self.ymir_dend, self.input_jstart: self.ymir_jstart}


    def compute_skip(self, fin):
        # check for levels of MiTCR output
        pass