'''
Various functions to make inferences about datasets downloaded from SRA.
Several functions require the output of the preprocessing steps in
sra_to_filtered_fastq.py.
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
sys.path.append("/mnt/data6A/functions")
import ut_functions
import json
import Run


def isPaired(pref, suffix="_sra.tsv", direc="."):
    '''
    Reads a file in the working directory nammed prefsuffix
    containing either "paired" or "single" and returns
    True for paired and False for single.

    Parameters
    ----------
    pref: str
        prefix of the input file name
    suffix: str
        suffix of the input file name
    direc: str
        path to directory containing the input file    
    '''
    ispaired = open("%s/%s%s" % (direc, pref, suffix)).readline().strip()
    return (True if ispaired == "paired" else False)


def getReadLength(pref, suffix="_readlength.tsv", direc="."):
    '''
    Reads a file containing only read length and returns this as an integer.
    
    Parameters
    ----------
    pref: str
        prefix of the input file name
    suffix: str
        suffix of the input file name
    direc: str
        path to directory containing the input file
    '''
    readlen = int(open("%s/%s%s" % (direc, pref,
                                    suffix)).readline().strip())
    return (readlen)


def getStrand(pref, suffix="_strandedness", intype="salmon",
              prog="hisat", direc="."):
    '''
    Reads a file containing only a dataframe describing library strandedness
    (generated with inferStrandednessSalmon below).
    formatted for one software package (currently only salmon is implemented)
    and returns a string formatted as part of the input statement
    for another package
    (currently implemented: hisat, bowtie, trinity).
    Assumes the input is stored in directory "direc" as 
    inputpackage.dir/prefsuffix.tsv
    e.g. for salmon.dir/SRR12345_strandedness.tsv
    intype = "salmon"
    pref = "SRR12345"
    suffix = "_strandedness.tsv"
    
    Parameters
    ----------
    pref: str
        prefix for input filename
    suffix: str
        suffix for input filename
    intype: str
        input format (salmon)
    prog: str
        output format (hisat, bowtie, trinity)
    direc: str
        directory containing the file to be parsed
    '''
    # path to the directory containing the raw salmon output
    path = "%s/%s.dir/%s%s.tsv" % (direc, intype, pref, suffix)
    
    # if the input file doesn't exist, assume unstranded
    try:
        df = pd.read_csv(path, sep="\t")
    except FileNotFoundError:
        return ""
    if prog == "hisat" or prog == "bowtie":
        cmd = "--rna-strandness"
    elif prog == "trinity":
        cmd = "--SS_lib_type"

    val = df['expected_format'].values[0]
    
    # conversion between salmon libtypes and libtypes used by
    # bowtie, hisat and trinity
    if intype == "salmon":
        if val == "ISR":
                return " %s RF " % cmd
        if val == "ISF":
            return " %s FR " % cmd
        if val == "SR":
            return " %s R " % cmd
        if val == "SF":
            return " %s F " % cmd
    return ""


def inferPairedSRA(pref, outfile, syst=""):
    '''
    Infers if an SRA dataset is paired or single end.
    Uses the NCBI fastq-dump --split-files function, which generates
    one output file for single end and two for paired end.
    Parameters
    ----------
    pref: str
        SRA ID for sample
    outfile: str
        output file - contains only "paired" or "single"
    
    '''
    T = tempfile.NamedTemporaryFile()
    Tnam = T.name

    # if there are multiple runs there will be a "." in the
    # outfile name - take the first run ID only

    if "." in pref:
        pref = pref.split(".")[0]


    # Execute fastq-dump using -X1 - take the first read or read pair only
    statement = """fastq-dump -X 1 --split-files %(pref)s """ % locals()
    ut_functions.writeCommand(statement, pref)
    Run.systemReRun(statement, syst)
    
    
    o = open(outfile, "w")
    # if the data is paired end there will be two files labelled as 
    # _1.fastq and _2.fastq, otherwise there will only be one
    # this one is usually _1.fastq but occasionally _2.fastq
    # Files are deleted once the type is known

    if os.path.exists("%s_2.fastq" % pref) and os.path.exists("%s_1.fastq" % pref):
        o.write("paired")
        os.unlink("%s_2.fastq" % pref)
        os.unlink("%s_1.fastq" % pref)
    
    elif os.path.exists("%s_1.fastq" % pref) or os.path.exists("%s_2.fastq" % pref):
        o.write("single")
        if os.path.exists("%s_1.fastq" % pref):
            os.unlink("%s_1.fastq" % pref)
        if os.path.exists("%s_2.fastq" % pref):
            os.unlink("%s_2.fastq" % pref)
    else:
        raise RuntimeError ("Couldn't downloaded SRA data to determine endedness for %s" % pref)
    o.close()
    

def inferReadLenFastQC(infiles, ispaired, outfile):
    '''
    Uses the output of FastQC software to determine the read length of the
    dataset.  Parses the fastqc_data file to determine this.  Where this is a
    range the top end of the range is used.
    Parameters
    ----------
    infiles: list
        list of output files from FastQC - this should always be a list of
        three: ['read_1' , 'read_2', 'single_end'] but for paired end the final
        file is empty and for single end the first two files are empty.
    ispaired: bool
        True if data is paired end else False
    outfile: str
        path to output file - will contain only an integer of maximum read
        length
    
    '''
    out = open(outfile, "w")
    if ispaired:
        dat1 = open(
            infiles[0].replace(
                ".html", "/fastqc_data.txt")).readlines()
        readlen1 = int(dat1[8].split("\t")[1].strip().split("-")[-1])
        dat2 = open(
            infiles[1].replace(".html",
                               "/fastqc_data.txt")).readlines()
        readlen2 = int(dat2[8].split("\t")[1].strip().split("-")[-1])
        out.write(str(min(readlen1, readlen2)))

    else:
        dat = open(infiles[2].replace(".html",
                                      "/fastqc_data.txt")).readlines()
        readlen = int(dat[8].split("\t")[1].strip().split("-")[-1])
        out.write(str(readlen))
    out.close()


def inferStrandednessSalmon(infile, outfile):
    '''
    Uses the output from Salmon software to determine the strandedness
    of an RNA-seq library.
    Running salmon with --libtype "A" outputs a json file containing
    the inferred library type.
    Parameters
    ----------
    infile: str
        path to the salmon lib_format_counts.json file
    outfile: str
        path to output file
    '''
    if os.path.exists(infile):
        # parse the json file into a dataframe
        j = json.load(open(infile, "r"))
        vals = list(j.values())
        cols = list(j.keys())
        D = pd.DataFrame(vals, index=cols).T
        D['sample'] = infile.split("/")[-1]
        D = D[["sample", "expected_format",
               "compatible_fragment_ratio",
               "num_compatible_fragments",
               "num_assigned_fragments",
               "num_consistent_mappings",
               "num_inconsistent_mappings",
               "MSF", "OSF", "ISF", "MSR",
               "OSR", "ISR", "SF", "SR",
               "MU", "OU", "IU", "U"]]

        D.to_csv(outfile, sep="\t", index=None)
    else:
        return None
