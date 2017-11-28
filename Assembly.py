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
sys.path.append("/mnt/data6A/functions")
import ut_functions


def runTrinity(infiles, pref, ispaired, strand, outfile):
    in1, in2, in3 = infiles
    outstem = "%s_trinity" % ".".join(outfile.split(".")[:-1])
    log = "%s.log" % outstem

    if ispaired == 'paired':
        statement = '''Trinity --seqType fq %(strand)s --max_memory 50G \
                       --CPU 5 --left %(in1)s --right %(in2)s \
                       --output %(outstem)s &>%(log)s''' % locals()
    else:
        statement = '''Trinity --seqType fq %(strand)s --max_memory 50G \
                       --CPU 5 --single %(in3)s \
                       --output %(outstem)s &>%(log)s''' % locals()
    ut_functions.writeCommand(statement, pref)
    os.system(statement)
    ut_functions.writeTime("Trinity", "end", pref)
