# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys
import logging
import numpy as np
from records import read_cont, read_tab1, read_list, read_text
import matplotlib.pyplot as plt
import pandas as pd


def split(file):
    """
    Split ``ENDF-6`` file  into MFMT sections.
    """
    import re
    pattern = ".{74}0.{5}\n?"
    text = open(file).read()
    U = re.split(pattern, text)
    return list(filter(None, U)) # remove empty lines



class Chunk:
    """
    MFMT section
    """

    @property
    def n(self):
        r"""
        Number of lines of the `ENDF-6` file
        """
        return len(self.text)

    @property
    def line(self):
        """
        Return current line.
        """
        return self.text[self.i]

    def __init__(self, text, mat=None, mf=None, mt=None):
        from io import StringIO
        import pandas as pd
        self.mat = int(text[66:70])
        self.mf = int(text[70:72])
        self.mt = int(text[72:75])
        # first argument must be any object with a read() method (such as a StringIO)
        self.text = pd.read_fwf(StringIO(text), widths=6*[11]+[4,2,3,5],
                                names=("c1","c2","l1","l2","n1","n2","mat","mf","mt","ns"),
                                header=None)
        self.i = 0

    def __iter__(self):
        """
        Make the file an iterator.
        """
        return self

    def next(self):
        """
        Yield a new line of the `ENDF-6` file.
        """
        line = self.line
        self.i += 1
        return line

    def move(self, steps):
        """
        Move backward/forward a number of lines in the ``ENDF-6`` text.

        Inputs:
        - :``steps``: :
            (scalar integer) number of lines to jump
        """
        self.i += int(steps)

    def read_cont(self):
        """
        Read ``ENDF-6`` ``CONT`` record in formatted fortran.

        Outputs:
            - :``out``: :
                (tuple) content of ``CONT`` record

        Found error in:
            - n-17-Cl-035.jeff32
            - n-3-Li-007.jeff32
            - n-63-Eu-152.jeff32
            - n-63-Eu-153.jeff32
            - n-64-Gd-155.jeff32
            - n-77-Ir-193.jeff32
            - n-90-Th-229.jeff32
            - n-94-Pu-238.jeff32
            - n-94-Pu-241.jeff32
            - n-94-Pu-242.jeff32
            - n-97-Bk-250.jeff32
            - n-98-Cf-254.jeff32
        """
        try:
            out = (float(self.line.c1), float(self.line.c2), float(self.line.l1),
                   float(self.line.l2), float(self.line.n1), float(self.line.n2))
            self.i += 1
            return out
        except:
            sys.exit("ERROR: line is not CONT.\n{}".format(self.line))


def list2dict(chunks):
    return 1

def process_section(text):
    mf = int(text[70:72])
    mt = int(text[72:75])
    if mf == 3:
        return read_mf3_mt(text)
    elif mf == 33:
        return read_mf33_mt(text)
    else:
        return None

def read_mf3_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2})
    T, i = read_tab1(str_list, i)
    out.update({"QM" : T.C1, "QI" : T.C2, "LR" : T.L2, "E" : T.x, "XS" : T.y})
    return out


def read_mf33_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    # Subsections are given as dictionary values.
    # Keys are MAT1*100+MT1
    out.update({"ZA" : C.C1, "AWR" : C.C2, "MTL" : C.L2, "SUB" : {}})
    for j in range(C.N2):
        sub = {}
        C, i = read_cont(str_list, i)
        sub.update({"XMF1" : C.C1, "XLFS1" : C.C2, "MAT1" : C.L1, "MT1" : C.L2})
        NC = C.N1
        NI = C.N2
        NCLIST = []
        for k in range(NC):
            C, i = read_cont(str_list, i)
            subsub = {"LTY" : C.L2}
            if subsub["LTY"] == 0:
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2,
                               "CI" : L.B[:L.N2], "XMTI" : L.B[L.N2:]})
                NCLIST.append(subsub)
            elif subsub["LTY"] in (1,2,3):
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2, "MATS" : L.L1, "MTS" : L.L2,
                               "XMFS" : L.B[0], "XLFSS" : L.B[1],
                               "EI" : L.B[2:2+L.N2], "WEI" : L.B[2+L.N2:]})
                NCLIST.append(subsub)
            NCLIST.append(subsub)
        sub.update({"NC" : NCLIST})
        NILIST = []
        for k in range(NI):
            L, i = read_list(str_list, i)
            subsub = {"LB" : L.L2}
            if subsub["LB"] in range(5):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                if subsub["LT"] == 0:
                    subsub.update({"Ek" : L.B[:subsub["NP"]], "Fk" : L.B[subsub["NP"]:]})
                else:
                    Nk = subsub["NP"] - subsub["LT"]
                    ARRk = L.B[:Nk]
                    ARRl = L.B[Nk:]
                    subsub.update({"Ek" : ARRk[:Nk/2], "Fk" : ARRk[Nk/2:],
                                   "El" : ARRl[:subsub["LT"]], "Fl" : ARRl[subsub["LT"]:]})
            elif subsub["LB"] == 5:
                subsub.update({"LS" : L.L1, "NT" : L.NPL, "NE" : L.N2,
                               "Ek" : L.B[:L.N2], "Fkk" : L.B[L.N2:]})
            elif subsub["LB"] == 6:
                subsub.update({"NT" : L.NPL, "NER" : L.N2})
                subsub.update({"NEC": (subsub["NT"]-1)/subsub["NER"]})
                subsub.update({"Ek" : L.B[:(subsub["NER"]+subsub["NEC"])],
                               "Fkk" : L.B[(subsub["NER"]+subsub["NEC"]):]})
            elif subsub["LB"] in (8,9):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                subsub.update({"Ek" : L.B[:subsub["NP"]], "Fk" : L.B[subsub["NP"]:]})
            NILIST.append(subsub)
        sub.update({"NI" : NILIST})
        out["SUB"].update({sub["MAT1"]*1000+sub["MT1"] : sub})
    return out

def write_mf33_mt(MF33):
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "MTL" : C.L2, "SUB" : []})
    for j in range(C.N2):
        sub = {}
        C, i = read_cont(str_list, i)
        sub.update({"XMF1" : C.C1, "XLFS1" : C.C2, "MAT1" : C.L1, "MT1" : C.L2})
        NC = C.N1
        NI = C.N2
        NCLIST = []
        for k in range(NC):
            C, i = read_cont(str_list, i)
            subsub = {"LTY" : C.L2}
            if subsub["LTY"] == 0:
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2,
                               "CI" : L.B[:L.N2], "XMTI" : L.B[L.N2:]})
                NCLIST.append(subsub)
            elif subsub["LTY"] in (1,2,3):
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2, "MATS" : L.L1, "MTS" : L.L2,
                               "XMFS" : L.B[0], "XLFSS" : L.B[1],
                               "EI" : L.B[2:2+L.N2], "WEI" : L.B[2+L.N2:]})
                NCLIST.append(subsub)
            NCLIST.append(subsub)
        sub.update({"NC" : NCLIST})
        NILIST = []
        for k in range(NI):
            L, i = read_list(str_list, i)
            subsub = {"LB" : L.L2}
            if subsub["LB"] in range(5):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                if subsub["LT"] == 0:
                    subsub.update({"Ek" : L.B[:subsub["NP"]], "Fk" : L.B[subsub["NP"]:]})
                else:
                    Nk = subsub["NP"] - subsub["LT"]
                    ARRk = L.B[:Nk]
                    ARRl = L.B[Nk:]
                    subsub.update({"Ek" : ARRk[:Nk/2], "Fk" : ARRk[Nk/2:],
                                   "El" : ARRl[:subsub["LT"]], "Fl" : ARRl[subsub["LT"]:]})
            elif subsub["LB"] == 5:
                subsub.update({"LS" : L.L1, "NT" : L.NPL, "NE" : L.N2,
                               "Ek" : L.B[:L.N2], "Fkk" : L.B[L.N2:]})
            elif subsub["LB"] == 6:
                subsub.update({"NT" : L.NPL, "NER" : L.N2})
                subsub.update({"NEC": (subsub["NT"]-1)/subsub["NER"]})
                subsub.update({"Ek" : L.B[:(subsub["NER"]+subsub["NEC"])],
                               "Fkk" : L.B[(subsub["NER"]+subsub["NEC"]):]})
            elif subsub["LB"] in (8,9):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                subsub.update({"Ek" : L.B[:subsub["NP"]], "Fk" : L.B[subsub["NP"]:]})
            NILIST.append(subsub)
        sub.update({"NI" : NILIST})
        out["SUB"].append(sub)
    return out

def pandas_interpolate(df, interp_column, method='zero'):
    # interp_column is a list
    ug = np.unique(list(df.index) + interp_column)
    df = df.reindex(ug)
    df = df.interpolate(method=method)
    df = df.reindex(interp_column)
    # Now transpose columns to rows and reinterpolate
    df = df.transpose()
    ug = np.unique(list(df.index) + interp_column)
    df = df.reindex(ug)
    df = df.interpolate(method=method)
    df = df.reindex(interp_column)
    df = df.transpose()
    return df

def extract_cov_mf33(MF33, mt=[102]):
    from sandy.cov import triu_matrix
    from sandy.functions import union_grid
    columns = ('MAT', 'MT', 'MAT1', 'MT1', 'COV')
    df = pd.DataFrame(columns=columns)
    for sub in MF33["SUB"].values():
        covs = []
        if len(sub["NI"]) == 0:
            continue
        for nisec in sub["NI"]:
            if nisec["LB"] == 5:
                Fkk = np.array(nisec["Fkk"])
                if nisec["LS"] == 0: # to be tested
                    cov = Fkk.reshape(nisec["NE"]-1, nisec["NE"]-1)
                else:
                    cov = triu_matrix(Fkk, nisec["NE"]-1)
                # add zero row and column at the end of the matrix
                cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                covs.append(pd.DataFrame(cov, index=nisec["Ek"], columns=nisec["Ek"]))
        if len(covs) == 0:
            continue
        uxx = union_grid(*[[list(cov.index) + list(cov.columns)] for cov in covs])
        cov = np.sum([ pandas_interpolate(cov,list(uxx)).as_matrix() for cov in covs ], 0)
        df = df.append({
                "MAT" : MF33['MAT'],
                "MT" : MF33['MT'],
                "MAT1" : sub['MAT1'],
                "MT1" : sub['MT1'],
                "COV" : cov
                }, ignore_index=True)
    return df

def merge_covs():

    pass



#columns = ('MAT', 'MT', 'MAT1', 'MT1', 'COV')
#tape = pd.DataFrame(columns=('MAT', 'MF', 'MT','SEC'))
#for chunk in split("H1.txt"):
#    tape = tape.append({
#        "MAT" : int(chunk[66:70]),
#        "MF" : int(chunk[70:72]),
#        "MT" : int(chunk[72:75]),
#        }, ignore_index=True)
tape = pd.DataFrame([[int(x[66:70]), int(x[70:72]), int(x[72:75]), x, int(x[66:70])*100000+int(x[70:72])*1000+int(x[72:75])] for x in split("H1.txt")],
        columns=('MAT', 'MF', 'MT','TEXT', 'ID'))
tape.set_index('ID', inplace=True)
tape = tape.set_index(['MAT','MF','MT']).sort_index() # Multi-indexing
tape['DATA'] = tape['TEXT'].apply(process_section)

#tape as a dictionary
#tape = {}
#for chunk in split("H1.txt"):
#    mat = int(chunk[66:70])
#    mf = int(chunk[70:72])
#    mt = int(chunk[72:75])
#    if mat not in tape:
#        tape.update({mat : {}})
#    if mf not in tape[mat]:
#        tape[mat].update({mf : {}})
#    tape[mat][mf].update({mt : chunk})

#process sections in dictionary
#for mat in tape:
#    for mf in tape[mat]:
#        for mt in tape[mat][mf]:
#            if mf == 3:
#                tape[mat][mf][mt] = read_mf3_mt(tape[mat][mf][mt])
#            elif mf ==33:
#                tape[mat][mf][mt] = read_mf33_mt(tape[mat][mf][mt])
O=read_mf3_mt(A[2])
P=read_mf33_mt(tape[125][33])
extract_cov_mf33(P)
write_mf33_mt(P)
XS=A[2].splitlines()
ii=0

CC, ii = read_cont(XS,ii)
read_tab1(XS,ii)
B=Chunk(A[-1])
import re
line = re.sub('\n', '', A[1])
AAA=re.findall('.{11}', 'line')
import pandas as pd
line = open('012503001').read()
line = re.sub('\n', '', line)
blocks = re.findall('.{11}', line)
df = pd.read_fwf('012503001', widths=6*[11]+[4,2,3,5],
                 names=("C1","C2","L1","L2","N1","N2","MAT","MF","MT","NS"),
                 header=None)
sys.exit()
df.fillna(0, inplace=True)

class File:
    """ General text file """

    @staticmethod
    def isfile(filename):
        """
        Check if file exists.

        Inputs:
            - filename :
                (string) path+name of the file
        """
        from os.path import isfile
        if not isfile(filename):
            logging.error("FILE : File '{}' does not exist".format(filename))
            sys.exit()

    def __init__(self, filename):
        """
        Initialize file.

        Inputs:
            - filename :
                (string) path+name of the file

        ..Important::
            Use encoding="ascii" and errors="surrogateescape" when opening the
            file in order to process files in ASCII compatible encoding such as
            `utf-8` and `latin-1`.
        """
        self.filename = filename
        self.f = open(self.filename, encoding="ascii", errors="surrogateescape")

    @property
    def filename(self):
        """
        Absolute file path + name.
        """
        return self._filename

    @filename.setter
    def filename(self, filename):
        from os.path import expanduser, abspath
        _f = abspath(expanduser(filename))
        File.isfile(_f)
        self._filename = _f

    @property
    def path(self):
        """
        Absolute file path.
        """
        from os.path import split
        return split(self.filename)[0]

    @property
    def name(self):
        """
        File basename.
        """
        from os.path import split
        return split(self.filename)[1]

    def read(self):
        """
        Load the content of input object `stream`.

        Outputs:
            - text :
                (list of strings) text in the file
        """
        try:
            text = self.f.readlines()
            return text
        except:
            logging.error("FILE : Cannot read file '{}'".format(self.filename))
            sys.exit()

    def load_yaml(self, msg="ERROR! Cannot yaml.load input file"):
        """
        Description
        ===========
        Load the content of input object *stream*.

        Outputs
        ======
         - *load*: text in yaml format
        """
        import yaml
        try:
            load = yaml.load(self.f)
            return load
        except yaml.YAMLError as exc:
            sys.exit(msg + " '{}'".format(self.filename))

#import pytest
#
#class TestFile:
#
#    fileyaml = '../test_objects/input.yaml'
#    filenotyaml = '../test_objects/input.not_yaml'
#    filenotexists = 'nonexistingfile'
#
#    def test_is_not_file(self):
#        from files import File
#        with pytest.raises(SystemExit):
#            File.isfile(self.filenotexists)
#
#    def test_is_file(self):
#        from files import File
#        File.isfile(self.fileyaml)
#
#    def test_init_file(self):
#        from files import File
#        F = File(self.fileyaml)
#        assert F.name == 'input.yaml'
#
#    def test_read_file(self):
#        from files import File
#        F = File(self.fileyaml)
#        assert F.read() == open(self.fileyaml).readlines()
#
#    def test_load_yaml_file(self):
#        from files import File
#        F = File(self.fileyaml)
#        load = F.load_yaml()
#        assert 'replace' in load
#        assert 'mat' in load['replace']
#        assert load['replace']['mat'] == 2631
#        assert 'endf1' in load['replace']
#        assert load['replace']['endf1'] == 'file1'
#        assert 'endf2' in load['replace']
#        assert load['replace']['endf2'] == 'file2'
#        assert 'mf' in load['replace']
#        assert load['replace']['mf'] == 3
#        assert 'mt' in load['replace']
#        assert load['replace']['mt'] == 102
#        assert 'emin' in load['replace']
#        assert load['replace']['emin'] == 1e6
#        assert 'emax' in load['replace']
#        assert load['replace']['emax'] == 10e6
#
#    def test_load_notyaml_file(self):
#        from files import File
#        F = File(self.filenotyaml)
#        with pytest.raises(SystemExit):
#            F.load_yaml()