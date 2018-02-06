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


def isPaired(pref, suffix="_sra.tsv"):
    '''
    Reads a file in the working directory nammed prefsuffix
    containing either "paired" or "single" and returns
    True for paired and False for single.
    '''
    ispaired = open("%s%s" % (pref, suffix)).readline().strip()
    return (True if ispaired == "paired" else False)


def getReadLength(pref, suffix="_readlength.tsv"):
    '''
    Reads a file in the working directory named prefsuffix
    containing only read length and returns this as an integer.
    '''
    readlen = int(open("%s%s" % (pref, suffix)).readline().strip())
    return (readlen)


def getStrand(pref, suffix="_strandedness", intype="salmon",
              prog="hisat"):

    path = "%s.dir/%s%s.tsv" % (intype, pref, suffix)
    try:
        df = pd.read_csv(path, sep="\t")
    except:
        return ""
    if prog == "hisat":
        cmd = "--rna-strandness"
    elif prog == "trinity":
        cmd = "--SS_lib_type"

    val = df['expected_format'].values[0]
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
    T = tempfile.NamedTemporaryFile()
    Tnam = T.name

    # if there are multiple runs there will be a "." in the
    # outfile name - take the first run ID only

    if "." in pref:
        pref = pref.split(".")[0]

    # writes a single alignment to file and checks if it has
    # 8 lines - paired or 4 lines - single.
    # -M 4 very short reads are filtered as some single end
    # datasets have a second very short read

    statement = """fastq-dump -X 1 --split-spot --stdout -Z \
    %(pref)s \
    | wc -l > %(Tnam)s""" % locals()

    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)
    os.system(statement)

    count = int(open(Tnam).readline().strip())
    o = open(outfile, "w")
    if count == 8:
        o.write("paired")
    elif count == 4:
        o.write("single")
    else:
        raise AssertionError ("can't infer endedness of %s" % pref)
    o.close()


def inferReadLenFastQC(infiles, ispaired, outfile):
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
    if os.path.exists(infile):
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
