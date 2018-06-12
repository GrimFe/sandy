#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:58:11 2018

@author: fiorito_l
"""
import pandas as pd
import os, re, sys, time, shutil, argparse, pdb

sab = pd.DataFrame.from_records([[48,9237,1,1,241,'uuo2'],
                                  [42,125,0,8,221,'tol'],
                                  [59,1425,0,1,221,'si'],
                                  [37,125,11,2,223,'pol'],
                                  [2,125,0,2,221,'ph'],
                                  [12,128,0,2,221,'pd'],
                                  [75,825,1,1,239,'ouo2'],
                                  [48,825,0,3,221,'osap'],
                                  [51,825,0,1,222,'ohw'],
                                  [46,825,3,1,237,'obeo'],
                                  [3,125,0,2,221,'oh2'],
                                  [13,128,0,2,221,'od2'],
                                  [52,1225,0,1,249,'mg'],
                                  [38,125,0,12,221,'mesi'],
                                  [10,125,1,2,221,'ice'],
                                  [7,125,12,1,225,'hzr'],
                                  [1,125,0,2,222,'lw'],
                                  [8,125,0,2,249,'hca'],
                                  [31,600,1,1,229,'gr'],
                                  [11,128,0,2,221,'dhw'],
                                  [59,2025,0,1,249,'cah'],
                                  [27,425,3,1,233,'bbeo'],
                                  [26,425,2,1,231,'be'],
                                  [60,1325,0,2,221,'asap']],
columns = ['matde','matdp','icoh','natom','mtref','ext'])



class PyNjoyError(Exception):
    """Exception indicating an error in PyNjoy."""



def TimeDecorator(foo):
    """
    Output the time a function takes
    to execute.
    """
    def wrapper(*args, **kwargs):
        t1 = time.time()
        out = foo(*args, **kwargs)
        t2 = time.time()
        print("Time to run function {}: {} sec".format(foo, t2-t1))
        return out
    return wrapper



class evalFile:

    def __init__(self, file):
        def read_float(x):
            try:
                return float(x[0] + x[1:].replace('+', 'E+').replace('-', 'E-'))
            except:
                return x
        widths = [11,11,11,11,11,11,4,2,3]
        columns = ["C1", "C2", "L1", "L2", "N1", "N2","MAT", "MF", "MT"]
        converters = dict(zip(columns[:6],[read_float]*6))
        tape = pd.read_fwf(file, widths=widths, names=columns, converters=converters).query("MAT>0 & MF>0 & MT>0")
        self.filename = os.path.basename(file)
        self.tape = self.endfFile = os.path.abspath(os.path.realpath(file))
        self.endf = self.mat = tape.MAT.iloc[0] # ENDF-6 MAT number
        info = tape.query("MF==1 & MT==451")
        self.za = int(info.C1.iloc[0])
        self.awr = info.C2.iloc[0]
        self.lrp = int(info.L1.iloc[0])
        if info.L2.iloc[0] == 1:
            self.fission = "yes"
        self.nlib = int(info.N1.iloc[0])
        self.mod = int(info.N2.iloc[0])
        self.lis = int(info.L1.iloc[1])
        self.liso = int(info.L2.iloc[1])
        self.metaStable = "true" if self.liso > 0 else "false"
        self.version = int(info.N2.iloc[1]) # Library format
        self.awi = info.C1.iloc[2]
        emax = info.C2.iloc[2]
        self.rel = int(info.L1.iloc[2])
        self.rev = int(info.N2.iloc[2])
        self.nsub = int(info.N1.iloc[2])
        if self.nsub == 10 or self.nsub == 12:
            self.neutron = "yes"
        tag = info.C1.iloc[4]
        try:
            pattern = re.search("(?P<Z>[0-9\s]*?)-(?P<SYM>[\w\s]*?)-(?P<AM>[0-9\s]*)", tag)
            self.tag = self.hmat = pattern.group("SYM").lower().strip() + pattern.group("AM").lower().strip()
        except:
            self.tag = self.hmat = tag.strip()
        if self.liso > 0:
            self.tag += str(self.liso)
        self.lab = info.C2.iloc[4]
        self.eval = info.L1.iloc[4]
        self.author = info.L2.iloc[4]
        self.dist = info.L1.iloc[5]
        self.rdate = info.L2.iloc[5]
        if not tape.query("MF==3 & MT==18").empty:
            self.totalFission = "yes"
#        self.pureAbsorber = "no" # ???
#        self.nis = 1 # ???
#        self.zai = 1 # ???
#        self.gamma = "no"
#        self.AWP0 = "yes"
#        self.dbrcnuclide = info.L1.iloc[4]
#        self.scattering = info.L1.iloc[4]
        if 2 in tape.MF.values:
            self.file2 = "yes"
#            self.resolved = info.L1.iloc[4]
#            self.resonance = info.L1.iloc[4]
#            self.unresolved = info.L1.iloc[4]
            ehRes = tape.query("MF==2 & MT==151").C2.iloc[2]
            if ehRes == emax: ehRes = 0.
            self.ehRes = ehRes
            self.thnmax = - self.ehRes if self.ehRes != 0 else 2.0E6
        if 3 in tape.MF.values:
            self.file3 = "yes"
        if 4 in tape.MF.values:
            self.file4 = "yes"
        if 5 in tape.MF.values:
            self.file5 = "yes"
        if 6 in tape.MF.values:
            self.file6 = "yes"
        if 7 in tape.MF.values:
            self.file7 = "yes"
        if 8 in tape.MF.values:
            self.file8 = "yes"
        if 9 in tape.MF.values or 10 in tape.MF.values:
            self.file9_10 = "yes"
        if 12 in tape.MF.values:
            self.file12 = "yes"
        if list(filter(lambda x:x>30, tape.MF.values)):
            self.covariances = "yes"
        if self.nsub == 12:
            raise NotImplementedError

    def to_frame(self):
        return pd.DataFrame.from_dict(self.__dict__, orient="index").T


class PyNjoy:

    def __init__(self, **kwargs):
#        self.sab = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sab.csv")
#        self.wdir = os.path.normpath(wdir)
#        self.evaluationName = os.path.basename(self.wdir)
        self.processes = None
        self.capsys = False
        self.NjoyExec = "njoy2016"
        self.wdir = os.getcwd()
        self.evaluationName = os.path.abspath("lib")
        self.iwt = 4
        self.legendre = 1
        self.legendregg = 6
        self.scatteringLaw = None
        self.eFiss = None
        self.branchingNG = None
        self.branchingN2N = None
        self.gstr = 0
        self.oldlib = None
        self.sgref = 1.0E10
        self.yields = None
        self.iburn = -1
        self.ip1opt = 0
        self.ifprod = 0
        self.jp1 = 0
        self.iverw = 4
        # General options
        self.iprint = 1  # by default set verbosity to max
        self.temps = [293.6]
        # RECONR/BROADR default options
        self.err = 0.005
        # PURR/UNRESR default options
        self.sig0 = [1e10] # default is infinite dilution only
        # THERMR default options
        self.iin = 1
        self.icoh = 1
        self.iform = 0
        self.natom = 1
        self.mtref = 221
        # ACER default options
        self.suffixes = [".03"]
        self.newfor = 1
        self.iopp = 1
        self.hk = ''
        # GROUPR default options
        self.ign = 2 # neutron group structure for GROUPR
        self.igg = 0
        self.iwt = 6
        self.lord = 1
        self.nstr = None
        self.gstr = None
        # ERRORR default options
        self.igne = 2 # neutron group structure for ERRORR
        self.iwte = 6
        self.irelco = 1
        self.irespr = 1
        # SKIP modules
        self.no_thermr = False
        self.no_broadr = False
        self.no_purr = False
        self.no_reconr = False
        self.unresr= False
        self.run_acer = False
        self.run_groupr = False
        self.run_errorr = False
        self.__dict__.update(kwargs)
        # Run checks
        if len(self.temps)> 10: raise PyNjoyError("cannot have more than 10 temperatures")



    def runNjoy(self, inputfile, cwd):
        """
        Run ``NJOY``.
        Pass input file to process's stdin (it must be encoded first).
        ..Important::
            In ``Python 3`` you need to convert string to bytes with a
            ``encode()`` function
        """
        from subprocess import Popen, PIPE
        if self.capsys:
            stderr = stdout = PIPE
        else:
            stderr = stdout = None
        process = Popen(self.NjoyExec,
                        shell=True,
                        cwd=cwd,
                        stdin=PIPE,
                        stdout=stdout,
                        stderr=stderr)
        inp = open(inputfile).read().encode()
        stdoutdata, stderrdata = process.communicate(inp)
        if process.returncode not in [0, 24]:
            raise PyNjoyError("NJOY exit status {}, cannot run njoy executable".format(process.returncode))

    def runBash(self, inputfile, cwd):
        from subprocess import Popen, PIPE
        if self.capsys:
            stderr = stdout = PIPE
        else:
            stderr = stdout = None
        process = Popen("bash {}".format(inputfile),
                        shell=True,
                        cwd=cwd,
                        stdin=PIPE,
                        stdout=stdout,
                        stderr=stderr)
        stdoutdata, stderrdata = process.communicate()
        if process.returncode not in [0, 24]:
            raise PyNjoyError("NJOY exit status {}, cannot run njoy executable".format(process.returncode))

    def pendf(self, **fileOptions):
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run pendf for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["endfFile"])): raise PyNjoyError("evaluation file " + kwargs["endfFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        os.makedirs(mydir, exist_ok=True)
        kwargs["nbDil"] = len(kwargs["sig0"])
        kwargs["textDil"] = " ".join(["%E"%dil for dil in kwargs["sig0"]])
        kwargs["nbTmp"] = len(kwargs["temps"])
        kwargs["textTmp"] = " ".join(["%E"%tmp for tmp in kwargs["temps"]])
        kwargs["htime"] = time.ctime(time.time())
        text_bash = """ln -sf %(endfFile)s tape20\n""" % kwargs
        text_data = """#!/bin/bash
set -e
cat > input_pendf.%(hmat)s << EOF
moder
20 -21
""" % kwargs
        kwargs["tapeEndf"] = -21

        # RECONR
        if kwargs["no_reconr"]:
            text_bash += "ln -sf %(pendfFile)s tape30\n" % kwargs
            text_data += """
moder
30 -22
""" % kwargs
        else:
            text_data += """
--
-- *********************************************************
-- Reconstruct XS from resonance parameters
-- *********************************************************
--
reconr
-21 -22
'pendf tape from %(evaluationName)s'/
%(mat)d 1 0/
%(err)E  0.  %(err)E/
'%(hmat)s from %(evaluationName)s at %(htime)s' /
0/
""" % kwargs
        kwargs["tapePendf"] = -22

        # BROADR (OPTIONAL)
        if not kwargs["no_broadr"]:
            text_data += """
--
-- *********************************************************
-- Perform Doppler-broadening
-- *********************************************************
--
broadr
%(tapeEndf)d %(tapePendf)d -23
%(mat)d %(nbTmp)d 0 0 0./
%(err)E %(thnmax)E %(err)E/
%(textTmp)s/
0/
""" % kwargs
            kwargs["tapePendf"] = -23

        # THERMR (OPTIONAL)
        if not kwargs["no_thermr"]:
            if self.scatteringLaw:
                text_data += """
moder
26 -27
--
-- *********************************************************
-- Add thermal scattering data
-- *********************************************************
--
thermr
-27 -23 -25
%(matLaw)d %(mat)d 32 %(nbTmp)d %(typeLaw)d %(elasOpt)d %(iform)d %(nbAtoms)d %(matsab_inc)d 0/
%(textTmp)s/
0.001 4.0
""" % kwargs
            else:
                text_data += """
--
-- *********************************************************
-- Add thermal scattering data for free-gas
-- *********************************************************
--
thermr
0 %(tapePendf)d -24
0 %(mat)d 20 %(nbTmp)d %(iin)d %(icoh)d %(iform)d %(natom)d 221 %(iprint)d
%(textTmp)s/
0.001 4.0
""" % kwargs
                kwargs["tapePendf"] = -24

        # PURR(OPTIONAL)
        if not kwargs["no_purr"]:
            text_data += """
purr
%(tapeEndf)d %(tapePendf)d -25
%(mat)d %(nbTmp)d %(nbDil)d 20 32/
%(textTmp)s/
%(textDil)s/
0/
""" % kwargs
            kwargs["tapePendf"] = -25

        # UNRESR (OPTIONAL)
        if kwargs["unresr"]:
            text_data += """
unresr
%(tapeEndf)d %(tapePendf)d -26
%(mat)d %(nbTmp)d %(nbDil)d 1/
%(textTmp)s/
%(textDil)s/
0/
""" % kwargs
            kwargs["tapePendf"] = -26


        text_data += """
moder
%(tapePendf)d 29
stop
EOF

""" % kwargs
        text_bash += """
%(NjoyExec)s < input_pendf.%(hmat)s
mv output output_pendf.%(hmat)s
mv tape29 %(filename)s.pendf
rm -f tape*
""" % kwargs
        text_data += text_bash
        inputfile = os.path.join(mydir, "run_pendf_{}.sh".format(kwargs["hmat"]))
        with open(inputfile,'w') as f:
            f.write(text_data)
        self.runBash(inputfile, mydir)
        return os.path.join(mydir, "%(filename)s.pendf" % kwargs)

    def groupr(self, **fileOptions):
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run groupr for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["endfFile"])): raise PyNjoyError("evaluation file " + kwargs["endfFile"] + " not found")
        if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise PyNjoyError("evaluation file " + kwargs["pendfFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        os.makedirs(mydir, exist_ok=True)
        kwargs["nbDil"] = len(kwargs["sig0"])
        kwargs["textDil"] = " ".join(["%E"%dil for dil in kwargs["sig0"]])
        kwargs["nbTmp"] = len(kwargs["temps"])
        kwargs["textTmp"] = " ".join(["%E"%tmp for tmp in kwargs["temps"]])
        if kwargs["ign"] == -99: # 1 GROUP (1E-5 to 2E7)
            kwargs["ign"] = 1
            kwargs["nstr"] = """1 /
1E-5 2E7 /
"""
        elif kwargs["ign"] == -98: # 1 GROUP (resonance integral)
            kwargs["ign"] = 1
            kwargs["nstr"] = """1 /
5E-1 1e5 /
"""
        kwargs["htime"] = time.ctime(time.time())
        text_bash = """
ln -sf %(endfFile)s tape20
ln -sf %(pendfFile)s tape29

%(NjoyExec)s < input_gendf.%(hmat)s
mv output output_gendf.%(hmat)s
mv tape30 %(filename)s.gendf
rm -f tape*
""" % kwargs
        text_data = """#!/bin/bash
set -e
cat > input_gendf.%(hmat)s << EOF
moder
20 -21
moder
29 -25
groupr
-21 -25 0 -26
%(mat)d %(ign)d %(igg)d %(iwt)d %(lord)d %(nbTmp)d %(nbDil)d %(iprint)d /
'%(hmat)s from %(evaluationName)s at %(htime)s' /
%(textTmp)s /
%(textDil)s /
""" % kwargs
        if kwargs["nstr"]:
            text_data += """%(nstr)s""" % kwargs
        for tmp in kwargs["temps"]:
            text_data += """3/
6/
0/
"""
        text_data += """0/
moder
-26 30
stop
EOF

""" % kwargs
        text_data += text_bash
        inputfile = os.path.join(mydir, "run_gendf_{}.sh".format(kwargs["hmat"]))
        with open(inputfile,'w') as f:
            f.write(text_data)
        self.runBash(inputfile, mydir)
        return os.path.join(mydir, "%(filename)s.gendf" % kwargs)

    def errorr(self, **fileOptions):
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run errorr for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["endfFile"])): raise PyNjoyError("evaluation file " + kwargs["endfFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        os.makedirs(mydir, exist_ok=True)
        if kwargs["igne"] == -99: # 1 GROUP (1E-5 to 2E7)
            kwargs["igne"] = 1
            kwargs["nstr"] = """1 /
1E-5 2E7 /
"""
        elif kwargs["igne"] == -98: # 1 GROUP (resonance integral)
            kwargs["igne"] = 1
            kwargs["nstr"] = """1 /
5E-1 1e5 /
"""
        kwargs["htime"] = time.ctime(time.time())
        for tmp,suff in zip(kwargs["temps"], kwargs["suffixes"]):
            kwargs.update({"tmp":tmp, "suff":suff})
            text_data = """
moder
20 -21
moder
29 -25
"""
            if "gendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["gendfFile"])): raise PyNjoyError("evaluation file " + kwargs["gendfFile"] + " not found")
                text_data += """
errorr
-21 0 -25 33 /
"""
            elif "pendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise PyNjoyError("evaluation file " + kwargs["pendfFile"] + " not found")
                text_data += """
errorr
-21 -25 0 33 /
"""
            else:
                raise PyNjoyError("user must define either pendfFile or gendfFile attribute")
            text_data += """%(mat)d %(igne)d %(iwte)d %(iprint)d %(irelco)d /
%(iprint)d %(tmp)E /
0 33 %(irespr)d %(lord)d /
""" % kwargs
            if kwargs["nstr"]:
                text_data += """%(nstr)s
""" % kwargs
            text_data += """
stop
"""
            text_bash = """#!/bin/bash
set -e
cat > input_errorr%(suff)s.%(hmat)s << EOF
""" % kwargs + text_data + """
EOF

ln -sf %(endfFile)s tape20
""" % kwargs
            if "gendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["gendfFile"])): raise PyNjoyError("evaluation file " + kwargs["gendfFile"] + " not found")
                text_bash += """
ln -sf %(gendfFile)s tape29
""" % kwargs
            elif "pendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise PyNjoyError("evaluation file " + kwargs["pendfFile"] + " not found")
                text_bash += """
ln -sf %(pendfFile)s tape29
""" % kwargs

            text_bash += """
%(NjoyExec)s < input_errorr%(suff)s.%(hmat)s
mv output output_errorr%(suff)s.%(hmat)s
mv tape33 %(filename)s%(suff)s.errorr
rm -f tape*
""" % kwargs
            inputfile = os.path.join(mydir, "run_errorr_{}{}.sh".format(kwargs["hmat"], suff))
            with open(inputfile,'w') as f:
                f.write(text_bash)
            self.runBash(inputfile, mydir)

    def acer(self, **fileOptions):
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run acer for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["endfFile"])): raise PyNjoyError("evaluation file " + kwargs["endfFile"] + " not found")
        if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise PyNjoyError("evaluation file " + kwargs["pendfFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        kwargs["htime"] = time.ctime(time.time())
        for tmp,suff in zip(kwargs["temps"], kwargs["suffixes"]):
            kwargs.update({"tmp":tmp, "suff":suff})
            text_data = """#!/bin/bash
set -e
cat > input_acer%(suff)s.%(hmat)s << EOF
moder
20 -21
moder
29 -25
acer
-21 -25 0 38 39
1 %(iprint)d 1 %(suff)s 0/
'%(hk)s'/
%(mat)d  %(tmp)E /
%(newfor)d %(iopp)d/
0.001/
-- acer / Check ACE files
-- 0 38 0 40 41
-- 7 1 1 -1/
-- /
stop
EOF

""" % kwargs
            text_data += """
ln -sf %(endfFile)s tape20
ln -sf %(filename)s.pendf tape29
%(NjoyExec)s < input_acer%(suff)s.%(hmat)s
mv output output_acer%(suff)s.%(hmat)s
mv tape38 %(filename)s%(suff)s.ace
mv tape39 %(filename)s%(suff)s.xsdir
rm -f tape*
""" % kwargs
            inputfile = os.path.join(mydir, "run_acer_{}{}.sh".format(kwargs["hmat"], suff))
            with open(inputfile,'w') as f:
                f.write(text_data)
            self.runBash(inputfile, mydir)
            xsdirfile = os.path.join(mydir, kwargs["filename"] + kwargs["suff"] + ".xsdir")
            if os.path.isfile(xsdirfile):
                with open(xsdirfile) as f:
                    xargs = f.read().split()
                xargs[2] = "%(filename)s%(suff)s.ace" % kwargs
                xargs[3] = "0" % kwargs
                xargs[0] = ".".join([str(kwargs['za'] + 300*kwargs["liso"]), xargs[0].split('.')[1]])
                with open(xsdirfile, 'w') as f:
                    f.write(" ".join(xargs))

    def run_modules(self, **fileOptions):
#        kwargs = dict(self.__dict__, **fileOptions)
        fileOptions["pendfFile"] = self.pendf(**fileOptions)
        if self.run_acer:
            self.acer(**fileOptions)
        if self.run_groupr:
            self.groupr(**fileOptions)
        if self.run_errorr:
            self.errorr(**fileOptions)



class evalLib(pd.DataFrame):

    @classmethod
    def from_file(cls, inputfile):
        """
        Populate dataframe using each row of a given inputfile.
        Example of inputfile:
            ENDF-1  [PENDF-1]
            ENDF-2  [PENDF-2]
            ...
            ENDF-N  [PENDF-2]
        Values in brackets [] are optional.
        """
        lines = open(inputfile).read().splitlines()
        evals = []
        for line in lines:
            groups = line.split()
            E = evalFile(groups[0])
            if len(groups) > 1:
                E.pendfFile = os.path.abspath(os.path.realpath(groups[1]))
            evals.append(E.to_frame())
        frame = pd.concat(evals, sort=False)
        return cls(frame)

    def run_njoy(self, **njoyOptions):
        import multiprocessing as mp
        njoy = PyNjoy(**njoyOptions)
        pool = mp.Pool(processes=njoy.processes)
        outs = [pool.apply_async(njoy.run_modules, kwds = {**row.to_dict()}) for i,row in self.iterrows()]
        outs = list(map(lambda x:x.get(), outs))
        if njoy.run_acer:
            xsdirFiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(njoy.evaluationName)) for f in fn if f.endswith(".xsdir")]
            aceFiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(njoy.evaluationName)) for f in fn if f.endswith(".ace")]
            mydir = os.path.join(njoy.evaluationName, "ace")
            os.makedirs(mydir, exist_ok=True)
            with open(os.path.join(mydir, os.path.basename(njoy.evaluationName)+".xsdir"), 'w') as f:
                for file in sorted(xsdirFiles):
                    f.write(open(file).read() + '\n')
            for file in sorted(aceFiles):
                shutil.move(file, os.path.join(mydir, os.path.basename(file)))

#lib = PyNjoy(iwt=-98, capsys=False)
#file = evalFile("1-H-3g.jeff33")
#file.pendfFile = lib.pendf(**file.__dict__)
#lib.errorr(**file.__dict__)
#sys.exit()

def is_valid_file(parser, arg, r=True, w=False, x=False):
    arg = os.path.abspath(os.path.realpath(os.path.normpath(arg)))
    if not os.path.isfile(arg):
        parser.error("File {} does not exist".format(arg))
    if r and not os.access(arg, os.R_OK):
        parser.error("File {} is not readable".format(arg))
    if w and not os.access(arg, os.W_OK):
        parser.error("File {} is not writable".format(arg))
    if x and not os.access(arg, os.X_OK):
        parser.error("File {} is not executable".format(arg))
    return arg

def is_valid_dir(parser, arg, mkdir=False):
    arg = os.path.abspath(os.path.realpath(os.path.normpath(arg)))
    if os.path.isdir(arg):
        return arg
    if mkdir:
        os.makedirs(arg, exist_ok=True)
    else:
        parser.error("Directory {} does not exist".format(arg))
    return arg

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run SANDY')
    parser.add_argument('-i','--inputfile',
                        type=lambda x: is_valid_file(parser, x),
                        required=True,
                        help="<Required> List of endf files to be processed (one file per line).")
    parser.add_argument('--processes',
                        type=int,
                        default=1,
                        help="Number of worker processes (default=1).")
    parser.add_argument('--capsys',
                        type=bool,
                        default=False,
                        help="Capture NJOY stderr and stdout (default=False).")
    parser.add_argument('--NjoyExec',
                        default='njoy2016',
                        help="NJOY executable (default=njoy2016).")
    parser.add_argument('--evaluationName',
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        default="lib",
                        metavar="lib",
                        help="Name of the evaluation.")
    parser.add_argument('--iprint',
                        type=int,
                        choices=range(2),
                        default=1,
                        help="NJOY verbosity: 0=min, 1=max (default=1).")
    parser.add_argument('--no-reconr',
                        action="store_true",
                        help="<Developer only> Skip RECONR module in the NJOY sequence.")
    parser.add_argument('--no-broadr',
                        action="store_true",
                        help="Skip BROADR module in the NJOY sequence.")
    parser.add_argument('--no-thermr',
                        action="store_true",
                        help="Skip THERMR module in the NJOY sequence.")
    parser.add_argument('--no-purr',
                        action="store_true",
                        help="Skip PURR module in the NJOY sequence.")
#    parser.add_argument('--no-pendf',
#                        action="store_true",
#                        help="Skip all NJOY modules necessary to produce a PENDF file.")
    parser.add_argument('--unresr',
                        action="store_true",
                        help="Replace PURR module with UNRESR in the NJOY sequence.")
    parser.add_argument('--run-groupr',
                        action="store_true",
                        help="Run GROUPR module in the NJOY sequence.")
    parser.add_argument('--run-acer',
                        action="store_true",
                        help="Run ACER module in the NJOY sequence.")
    parser.add_argument('--run-errorr',
                        action="store_true",
                        help="Run ERRORR module in the NJOY sequence.")
    parser.add_argument('--temps',
                        type=float,
                        default = [293.6],
                        nargs='+',
                        help="Temperature values (default=[293.6]).")
    parser.add_argument('--sig0',
                        type=float,
                        default=[1e10],
                        nargs='+',
                        help="Sigma0 values (default=[1e10]).")
    parser.add_argument('--err',
                        type=float,
                        default=0.005,
                        help="Fractional tolerance for RECONR and BROADR (default=0.005).")
    parser.add_argument('--ign',
                        type=int,
                        default=2,
                        help="Neutron group structure option for GROUPR (default=2 : csewg 239-group structure).")
    parser.add_argument('--igne',
                        type=int,
                        default=2,
                        help="Neutron group structure option for ERRORR (default=2 : csewg 239-group structure).")
    parser.add_argument('--iwt',
                        type=int,
                        default=6,
                        help="Weight function option for GROUPR (default=6 : Maxwellian - 1/E - fission - fusion).")
    parser.add_argument('--iwte',
                        type=int,
                        default=6,
                        help="Weight function option for ERRORR (default=6 : Maxwellian - 1/E - fission - fusion).")
    parser.add_argument('--suffixes',
                        type=str,
                        default=[".03"],
                        nargs='+',
                        help="Suffixes for ACE files, as many as temperature values (default=[\".03\"]]).")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="")
    args = parser.parse_args()
    evalLib.from_file(args.inputfile).run_njoy(**vars(args))