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

def runSalmon(infiles, outfile, ispaired, pref,
              transpath, salmon_opts="", logfile="salmon.log",
              syst=""):
    print (transpath)
    if os.path.exists(transpath):
        in1, in2, in3 = infiles
        if ispaired:
            statement = '''salmon quant -i %(transpath)s  \
            %(salmon_opts)s --output salmon.dir/%(pref)s \
            -1 %(in1)s -2 %(in2)s 2>>%(logfile)s''' % locals()
        else:
            statement = '''salmon quant -i %(transpath)s  \
            %(salmon_opts)s --output salmon.dir/%(pref)s \
            -r %(in3)s 2>>%(logfile)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
    else:
        # touch the output file if there is no transcriptome
        # to compare to
        try:
            os.mkdir("salmon.dir/%(stem)s" % locals())
            pathlib.Path(outfile).touch()
        except:
            pathlib.Path(outfile).touch()
