'''
Functions to run different assembly software within a pipeline using
user defined parameters.
'''

from ruffus import *
import ruffus.cmdline as cmdline
import pathlib
import ut_functions
import Run
import os


def runTrinity(infiles, pref, ispaired, strand, outfile, threads=1,
               mem="50G", syst=""):
    '''
    Generates and runs a statement to run the Trinity assembly software.
    
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


def runSpades(infiles, pref, ispaired, outfile, threads=1,
              mem=4, syst=""):
    '''
    Generates and runs a statement to run the rnaspades assembly software
    
    
    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]

    pref: str
        filename prefix (used for the log files)
    ispaired: bool
        True if paired end else False
    outfile: str
        path to output file, this is a placeholder to show the task has run
    threads: int
        number of threads with which to run spades
    mem: str
        maximum memory to pass to spades (formated as an integer)
    syst: str
        system option to run statements on different systems    
    '''
    in1, in2, in3 = infiles
    ut_functions.writeTime("Spades", "start", pref)

    outstem1 = "%s_rnaspades" % ".".join(outfile.split(".")[:-1])
    outstem2 = "%s_rnaviralspades" % ".".join(outfile.split(".")[:-1])
    log1 = "%s.log" % outstem1
    log2 = "%s.log" % outstem2
    if ispaired:
        f1 = in1.replace(".fq.1", ".1.fq")
        f2 = in2.replace(".fq.2", ".2.fq")
        if not os.path.exists(f1):
            os.symlink(in1.split("/")[-1], f1)
        if not os.path.exists(f2):
            os.symlink(in2.split("/")[-1], f2)
        statement1 = '''rnaspades.py -1 %(f1)s -2 %(f2)s \
                       -o %(outstem1)s -t %(threads)s -m %(mem)s \
                       2>%(log1)s''' % locals()
        statement2 = '''rnaviralspades.py -1 %(f1)s -2 %(f2)s \
                       -o %(outstem2)s -t %(threads)s -m %(mem)s \
                       2>%(log2)s''' % locals()
    else:
        statement1 = '''rnaspades.py -s %(in3)s\
                       -o %(outstem1)s -t %(threads)s -m %(mem)s \
                       2>%(log1)s''' % locals()
        statement2 = '''rnaviralspades.py -s %(in1)s\
                       -o %(outstem2)s -t %(threads)s -m %(mem)s \
                       2>%(log2)s''' % locals()
    ut_functions.writeCommand(statement1, pref)
    ut_functions.writeCommand(statement2, pref)
    Run.systemRun(statement1, syst)
    Run.systemRun(statement2, syst)
    ut_functions.writeTime("Spades", "end", pref)
    pathlib.Path(outfile).touch()
