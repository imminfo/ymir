import gzip


class RepertoireConverter:
    """
    Parent class for all other converters. Just rewrite functions to make them do what you want to.
    """

    def __init__(self):
        self.ymir_columns = {}  # dict <ymir col name> : <index of related column in the input file>
        self.init_ymir_names()

        self.input_columns = {}  # dict <input col name> : <ymir col name>
        self.init_input_column_names()
        self.init_input_columns()

        self.set_base()


    def init_ymir_names(self):
        """
        Initialise Ymir's format column names and separator characters.
        """

        self.ymir_nuc = "Nucleotide sequence"
        self.ymir_aa = "Amino acid sequence"
        self.ymir_vgene = "Variable"
        self.ymir_dgene = "Diversity"
        self.ymir_jgene = "Joining"
        self.ymir_vend = "V end"
        self.ymir_dstart = "D start"
        self.ymir_dend = "D end"
        self.ymir_jstart = "J start"
        self.ymir_col_sep = "\t"
        self.ymir_gene_sep = ","

        self.ymir_columns = {self.ymir_nuc: -1, self.ymir_aa: -1,
                             self.ymir_vgene: -1, self.ymir_dgene: -1, self.ymir_jgene: -1,
                             self.ymir_vend: -1, self.ymir_dstart: -1, self.ymir_dend: -1, self.ymir_jstart: -1}

        self.ymir_columns_sorted = [self.ymir_nuc, self.ymir_aa,
                                    self.ymir_vgene, self.ymir_dgene, self.ymir_jgene,
                                    self.ymir_vend, self.ymir_dstart, self.ymir_dend, self.ymir_jstart]


    def init_input_column_names(self):
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


    def init_input_columns(self):
        self.input_columns = {self.input_nuc: self.ymir_nuc, self.input_aa: self.ymir_aa,
                              self.input_vgene: self.ymir_vgene, self.input_dgene: self.ymir_dgene, self.input_jgene: self.ymir_jgene,
                              self.input_vend: self.ymir_vend, self.input_dstart: self.ymir_dstart, self.input_dend: self.ymir_dend, self.input_jstart: self.ymir_jstart}


    def convert(self, file_in, file_out):
        """
        Parse and convert input file to the default Ymir's format file.

        :param file_in: path to the input file
        :param file_out: path to the output file
        """

        skip = 0
        self.gzip = False
        open_fun = lambda x: gzip.open(x) if x.endswith(".gz") else open(x)
        if (file_in.endswith(".gz")):
            self.gzip = True

        with open_fun(file_in) as fin:
            skip = self.compute_skip(fin)

        with open_fun(file_in) as fin:
            with open(file_out, "w") as fout:
                if skip:
                    for line in fin:
                        skip -= 1
                        if skip == 0: break

                self.parse_header(fin.readline())
                self.write_header(fout)
                for line in fin:
                    if line.strip():
                        fout.write(self.parse_line(line) + "\n")


        return True


    def parse_line(self, line):
        """
        Parse the input line and return new line in the Ymir format to write
        to the output file.

        :param line: a line from the input file
        :return: line in Ymir's format
        """

        if self.gzip: line = line.decode()

        words = line.strip().split(self.input_col_sep)
        out_words = ["" for x in range(len(self.ymir_columns))]

        for i, key in enumerate(self.ymir_columns_sorted):
            out_words[i] = words[self.ymir_columns[key]]

        if self._base == 0:
            out_words[5] = str(int(out_words[5]) + 1)
            out_words[6] = str(int(out_words[6]) + 1)
            out_words[7] = str(int(out_words[7]) + 1)
            out_words[8] = str(int(out_words[8]) + 1)

        return self.ymir_col_sep.join(out_words)


    def parse_header(self, header):
        if self.gzip: header = header.decode()

        words = header.strip().split(self.input_col_sep)
        for i, w in enumerate(words):
            if w in self.input_columns:
                self.ymir_columns[self.input_columns[w]] = i


    def write_header(self, fout):
        fout.write(self.ymir_col_sep.join([self.ymir_nuc, self.ymir_aa,
                                           self.ymir_vgene, self.ymir_dgene, self.ymir_jgene,
                                           self.ymir_vend, self.ymir_dstart, self.ymir_dend, self.ymir_jstart]) + "\n")


    def compute_skip(self, fin):
        return 0


    def set_base(self):
        self._base = 0


class YmirConverter:
    def __init__(self):
        pass


    def convert(self, file_in, file_out):
        pass


    def set_base(self):
        self._base = 1


class tcRConverter (RepertoireConverter):
    def __init__(self):
        RepertoireConverter.__init__(self)


    def init_input_column_names(self):
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


    def set_base(self):
        self._base = 0


class MiTCRConverter (RepertoireConverter):
    def __init__(self):
        RepertoireConverter.__init__(self)


    def init_input_column_names(self):
        """
        Initialise input files' column names and separator characters.
        """

        # self.input_nuc = "CDR3 nucleotide sequence"
        # self.input_aa = "CDR3 amino acid sequence"
        # self.input_vgene = "V segments"
        # self.input_dgene = "D segments"
        # self.input_jgene = "J segments"
        # self.input_vend = "Last V nucleotide position"
        # self.input_dstart = "First D nucleotide position"
        # self.input_dend = "Last D nucleotide position"
        # self.input_jstart = "First J nucleotide position"
        # self.input_col_sep = "\t"
        # self.input_gene_sep = ", "
        self.input_nuc = "CDR3.nucleotide.sequence"
        self.input_aa = "CDR3.amino.acid.sequence"
        self.input_vgene = "V.segments"
        self.input_dgene = "D.segments"
        self.input_jgene = "J.segments"
        self.input_vend = "Last.V.nucleotide.position"
        self.input_dstart = "First.D.nucleotide.position"
        self.input_dend = "Last.D.nucleotide.position"
        self.input_jstart = "First.J.nucleotide.position"
        self.input_col_sep = "\t"
        self.input_gene_sep = ", "


    # def compute_skip(self, fin):
    #     # check for levels of MiTCR output
    #     pass


    def set_base(self):
        self._base = 0


class MiGECConverter (RepertoireConverter):
    def __init__(self):
        RepertoireConverter.__init__(self)


    def init_input_column_names(self):
        """
        Initialise input files' column names and separator characters.
        """

        self.input_nuc = "CDR3 nucleotide sequence"
        self.input_aa = "CDR3 amino acid sequence"
        self.input_vgene = "V segments"
        self.input_dgene = "D segments"
        self.input_jgene = "J segments"
        self.input_vend = "Last V nucleotide position"
        self.input_dstart = "First D nucleotide position"
        self.input_dend = "Last D nucleotide position"
        self.input_jstart = "First J nucleotide position"
        self.input_col_sep = "\t"
        self.input_gene_sep = ","


    def set_base(self):
        self._base = 0


class VDJtoolsConverter (RepertoireConverter):
    def __init__(self):
        RepertoireConverter.__init__(self)


    def init_input_column_names(self):
        """
        Initialise input files' column names and separator characters.
        """

        self.input_nuc = "cdr3nt"
        self.input_aa = "cdr3aa"
        self.input_vgene = "v"
        self.input_dgene = "d"
        self.input_jgene = "j"
        self.input_vend = "VEnd"
        self.input_dstart = "DStart"
        self.input_dend = "DEnd"
        self.input_jstart = "JStart"
        self.input_col_sep = "\t"
        self.input_gene_sep = ","


    def set_base(self):
        self._base = 0
