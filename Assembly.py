'''
Functions to run different assembly software within a pipeline using
user defined parameters.
'''

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
    '''
    Generates and runs a statement to run the Trinity assembly software:
        https://github.com/trinityrnaseq/trinityrnaseq/wiki
    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]

    pref: str
        filename prefix (used for the log files)
    ispaired: bool
        True if paired end else False
    strand: str
        library strandedness string (e.g. --SS_lib_type SF) - this can be
        generated from the salmon output using Infer.getStrand()
    outfile: str
        path to output file, this is a placeholder to show the task has run
        as Trinity generates many output files
    threads: int
        number of threads with which to run Trinity
    mem: str
        maximum memory to pass to Trinity (formated as e.g 50G)
    syst: str
        system option to run statements on different systems    
    '''
    in1, in2, in3 = infiles
    outstem = "%s_trinity" % ".".join(outfile.split(".")[:-1])
    log = "%s.log" % outstem

    if ispaired:
        statement = '''Trinity --seqType fq %(strand)s --max_memory %(mem)s \
                       --CPU %(threads)s --left %(in1)s --right %(in2)s --full_cleanup \
                       --output %(outstem)s 2>%(log)s''' % locals()
    else:
        statement = '''Trinity --seqType fq %(strand)s --max_memory %(mem)s \
                       --CPU %(threads)s --single %(in3)s --full_cleanup \
                       --output %(outstem)s 2>%(log)s''' % locals()
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)
    ut_functions.writeTime("Trinity", "end", pref)
    pathlib.Path(outfile).touch()
