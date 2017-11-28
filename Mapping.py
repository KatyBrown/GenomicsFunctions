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


def mapReadsHisat(infiles, outfile, genomepath, genomename, strand,
                  ispaired, options, mismatches, pref, double=False):
    if double is True:
        hisat =  "%s/%s/hisat/%s" % (genomepath, genomename,
                                     genomename)
    else:
        hisat = "%s/hisat/%s" % (genomepath, genomename)

    #log = outfile.split("_")[-2].split("/")[-1]
    log = "logs.dir/%s_%s_mapping.log" % (pref, genomename)
    #met = outfile.split("_")[-2].split("/")[-1]
    met = "logs.dir/%s_%s_mapping.met" % (pref, genomename)

    m = -1 * (float(mismatches) * 6)
    mstring = '''--ignore-quals --score-min L,0,%s''' % m

    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        statement = '''hisat2 -x %(hisat)s %(strand)s -1 %(in1)s -2 %(in2)s \
        %(mstring)s --met-file %(met)s %(options)s 2>%(log)s|\
        samtools view -b > %(outfile)s''' % locals()
    else:
        in1 = infiles[2]
        statement = '''hisat2 -x %(hisat)s %(strand)s -U %(in1)s \
        %(mstring)s --met-file %(met)s %(options)s 2>%(log)s|\
        samtools view -b > %(outfile)s''' % locals()
    ut_functions.writeCommand(statement, pref)
    os.system(statement)


def mapReadsBowtie(infiles, outfile, genomepath, genomename, strand,
                   ispaired, options, mismatches, pref, double=False):
    if double is True:
        bowtie =  "%s/%s/bowtie/%s" % (genomepath, genomename)
    else:
        bowtie = "%s/bowtie/%s" % (genomepath, genomename)

    log = "logs.dir/%s_%s_mapping.log" % (pref, genomename)
    met = outfile.split("_")[-2].split("/")[-1]
    met = "logs.dir/%s_%s_mapping.met" % (met, genomename)

    m = -1 * (float(mismatches) * 6)
    mstring = '''--ignore-quals --score-min L,0,%s''' % m

    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        statement = '''bowtie2 -x %(bowtie)s -1 %(in1)s -2 %(in2)s \
        --met-file %(met)s %(options)s %(mstring)s 2>%(log)s|\
        samtools view -b > %(outfile)s''' % locals()
    else:
        in1 = infiles[2]
        statement = '''bowtie2 -x %(bowtie)s -U %(in1)s \
        --met-file %(met)s %(options)s %(mstring)s 2>%(log)s|\
        samtools view -b > %(outfile)s''' % locals()
    ut_functions.writeCommand(statement, pref)
    os.system(statement)


def filterMappedReads(infile, outfiles, ispaired, genomename, pref,
                q=10):
    out1, out2, out3 = outfiles
    tempout = infile.replace(".bam", ".temp")
    tempout2 = infile.replace(".bam", ".temp2")

    log = "logs.dir/%s_%s_filtering.log" % (pref, genomename)

    if ispaired:
        statement = '''samtools view -F4 -q%(q)s -f1 -f2 -b \
                       -U %(tempout)s  %(infile)s \
                       > %(tempout2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        os.system(statement)
        statement = '''samtools sort -n %(tempout)s |\
                       samtools fastq -i - -1 %(out1)s -2 %(out2)s \
                       &>%(log)s;\
                       rm -rf %(tempout)s %(tempout2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        os.system(statement)
        pathlib.Path(out3).touch()

    else:
        tempout = infile.replace(".bam", ".temp")
        statement = '''samtools view -F4 -q%(q)s -b\
                       -U %(tempout)s %(infile)s \
                       > %(tempout2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        os.system(statement)
        statement = '''samtools sort -n %(tempout)s |
                       samtools fastq -i - -0 %(out3)s &>%(log)s;
                       rm -rf %(tempout)s %(tempout2)s''' % locals()
        os.system(statement)
        pathlib.Path(out1).touch()
        pathlib.Path(out2).touch()


def cleanBam(infile, outfile, ispaired, pref,
             q=10):
    if ispaired == "paired":
        statement = '''samtools view -b -F4 -q%(q)s -f1 -f2 %(infile)s >\
                       %(outfile)s''' % locals()
    else:
        statement = '''samtools view -b -F4 -q%(q)s %(infile)s >\
                       %(outfile)s''' % locals()

    ut_functions.writeCommand(statement, pref)
    os.system(statement)
