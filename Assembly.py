from ruffus import *
import ruffus.cmdline as cmdline
import pandas as pd
import numpy as np
import os
import re
import sys
import pathlib
import tempfile
import subprocess
import shutil
import ut_functions
import Run


def runTrinity(infiles, pref, ispaired, strand, outfile, threads=1,
               mem="50G", syst=""):

    in1, in2, in3 = infiles
    outstem = "%s_trinity" % ".".join(outfile.split(".")[:-1])
    log = "%s.log" % outstem

    if ispaired:
        statement = '''Trinity --seqType fq %(strand)s --max_memory %(mem)s \
                       --CPU %(threads)s --left %(in1)s --right %(in2)s \
                       --output %(outstem)s &>%(log)s''' % locals()
    else:
        statement = '''Trinity --seqType fq %(strand)s --max_memory %(mem)s \
                       --CPU %(threads)s --single %(in3)s \
                       --output %(outstem)s &>%(log)s''' % locals()
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)
    ut_functions.writeTime("Trinity", "end", pref)
    pathlib.Path(outfile).touch()
