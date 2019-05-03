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
import Run
import gzip
import time


def getEBIAddress(run):
    prefix = run[0:6]
    address = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s" % prefix
    if len(run) == 9:
        address = "%s/%s/%s" % (address, run, run)
    elif len(run) == 10:
        suffix = "00%s" % run[-1]
        address = "%s/%s/%s/%s"  % (address, suffix, run, run)
    else:
        raise ValueError ("run ID not recognised")
    return (address)
    
    
def getSRA(pref, ispaired, log, sra_opts,
           outfiles, nreads_download=1000000000,
           syst="", sra='ncbi', asperadir=""):
    '''
    Downloads data "pref" from SRA.
    Pref can be either a single SRA ID or a list seperated by "."
    ispaired is True for paired, False for single ended.
    log is a log file
    sra_opts is any additional string of options to pass to SRA
    outfiles is three outfiles - paired end one, paired end two
    and single end.
    nreads_download is the maximum number of reads to download.
    '''
    out1, out2, out3 = outfiles
    x = 0
    while True and x < 25:
        if ispaired:
            # merge all the runs into one big fasta file
            if "." in pref:
                runs = pref.split(".")
                statements = []
                l = open(log, "w")
                l.close()
                # only take a maximum of nspots reads or pairs
                nspots = "-X %i" % (nreads_download / len(runs))
                for run in runs:
                    if sra == "ncbi":
                        statement = '''fastq-dump %(nspots)s %(sra_opts)s \
                        --split-files -v \
                        --outdir fastqs.dir &>>%(log)s\
                        %(run)s ''' % locals()
                        statements.append(statement)
                    elif sra == "ebi":
                        nrows = nreads_download * 4
                        address = getEBIAddress(run)
                        
                        statement = """curl -s %(address)s_1.fastq.gz | \
                                       gunzip - | \
                                       head -n %(nrows)s > \
                                       fastqs.dir/%(run)s_1.fastq;
                                       curl -s %(address)s_2.fastq.gz | \
                                       gunzip - | \
                                       head -n %(nrows)s > \
                                       fastqs.dir/%(run)s_2.fastq""" % locals()
                        statements.append(statement)
                    elif sra == 'aspera_ebi':
                        nrows = nreads_download * 4
                        address = getEBIAddress(run)
                        asp = address.replace("ftp://ftp.sra.ebi.ac.uk", "")
                        tempname1 = "run_%s_1.fastq.gz" % run
                        tempname2 = "run_%s_2.fastq.gz" % run
                        statement = """
                        ascp -QT -l 300m -P33001 -i %(asperadir)s/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:%(asp)s_1.fastq.gz fastqs.dir;
                        mv fastqs.dir/%(run)s_1.fastq.gz fastqs.dir/%(tempname1)s;\
                        zcat fastqs.dir/%(tempname1)s | head -%(nrows)s > fastqs.dir/%(run)s_1.fastq;
                        ascp -QT -l 300m -P33001 -i %(asperadir)s/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:%(asp)s_2.fastq.gz fastqs.dir;
                        mv fastqs.dir/%(run)s_2.fastq.gz fastqs.dir/%(tempname2)s;\
                        zcat fastqs.dir/%(tempname2)s | head -%(nrows)s > fastqs.dir/%(run)s_2.fastq;
                        rm -rf fastqs.dir/%(tempname1)s;
                        rm -rf fastqs.dir/%(tempname2)s""" % locals()
                        statements.append(statement)

                statement = "; ".join(statements)
                catlist1 = " ".join(["fastqs.dir/%s_1.fastq" % run
                                     for run in runs])
                catlist2 = " ".join(["fastqs.dir/%s_2.fastq" % run
                                     for run in runs])
                # merge the fastq files then delete the originals
                o1 = out1.replace(".gz", "")
                o2 = out2.replace(".gz", "")
                statement += "; cat %(catlist1)s > %(o1)s; gzip -f %(o1)s" % locals()
                statement += "; cat %(catlist2)s > %(o2)s; gzip -f %(o2)s" % locals()
                statement += "; rm -rf  %(catlist1)s" % locals()
                statement += "; rm -rf %(catlist2)s" % locals()
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)

            else:
                run = pref
                nspots = "-X %i" % (nreads_download)
                if sra == "ncbi":
                    statement = '''fastq-dump %(nspots)s %(sra_opts)s \
                    --split-files --gzip -v \
                    --outdir fastqs.dir &>>%(log)s\
                    %(pref)s ''' % locals()
                elif sra == "ebi":
                    nrows = nreads_download * 4
                    address = getEBIAddress(run)
                    statement = """curl -s %(address)s_1.fastq.gz | \
                                   gunzip - | \
                                   head -n %(nrows)s \
                                   > fastqs.dir/%(run)s_1.fastq;
                                   curl -s %(address)s_2.fastq.gz | \
                                   gunzip - | \
                                   head -n %(nrows)s \
                                   > fastqs.dir/%(run)s_2.fastq""" % locals()
  
                elif sra == "aspera_ebi":
                    nrows = nreads_download * 4
                    address = getEBIAddress(run)
                    asp = address.replace("ftp://ftp.sra.ebi.ac.uk", "")
                    tempname1 = "run_%s_1.fastq.gz" % run
                    tempname2 = "run_%s_2.fastq.gz" % run
                    statement = """
                    ascp -QT -l 300m -P33001 -i %(asperadir)s/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:%(asp)s_1.fastq.gz fastqs.dir;
                    mv fastqs.dir/%(run)s_1.fastq.gz fastqs.dir/%(tempname1)s;\
                    zcat fastqs.dir/%(tempname1)s | head -%(nrows)s > fastqs.dir/%(run)s_1.fastq;
                    ascp -QT -l 300m -P33001 -i %(asperadir)s/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:%(asp)s_2.fastq.gz fastqs.dir;
                    mv fastqs.dir/%(run)s_2.fastq.gz fastqs.dir/%(tempname2)s;\
                    zcat fastqs.dir/%(tempname2)s | head -%(nrows)s > fastqs.dir/%(run)s_2.fastq;
                    rm -rf fastqs.dir/%(tempname1)s;
                    rm -rf fastqs.dir/%(tempname2)s""" % locals()
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)
                statement = 'gzip -f fastqs.dir/%s_1.fastq' % run
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)
                statement = 'gzip -f fastqs.dir/%s_2.fastq' % run
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst) 
            pathlib.Path(out3).touch()
    
        else:
            if "." in pref:
                runs = pref.split(".")
                statements = []
                l = open(log, "w")
                l.close()
                nspots = "-X %i" % (nreads_download / len(runs))
                for run in runs:
                    if sra == "ncbi":
                        statement = '''fastq-dump %(nspots)s %(sra_opts)s \
                         -v \
                        --outdir fastqs.dir &>>%(log)s\
                        %(run)s ''' % locals()
                        statements.append(statement)
                    elif sra == "ebi":
                        nrows = nreads_download * 4
                        address = getEBIAddress(run)
                        statement = """curl -s %(address)s.fastq.gz | \
                                       gunzip - | \
                                       head -n %(nrows)s \
                                       > fastqs.dir/%(run)s.fastq""" % locals()
                        statements.append(statement)
                    elif sra == "aspera_ebi":
                        nrows = nreads_download * 4
                        address = getEBIAddress(run)
                        asp = address.replace("ftp://ftp.sra.ebi.ac.uk", "")
                        tempname = "run_%s.fastq.gz" % run
                        statement = """
                        ascp -QT -l 300m -P33001 -i %(asperadir)s/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:%(asp)s.fastq.gz fastqs.dir;
                        mv fastqs.dir/%(run)s.fastq.gz fastqs.dir/%(tempname)s;\
                        zcat fastqs.dir/%(tempname)s | head -%(nrows)s > fastqs.dir/%(run)s.fastq;
                        rm -rf fastqs.dir/%(tempname)s""" % locals()
                        statements.append(statement)
                statement = "; ".join(statements)
                catlist = " ".join(["fastqs.dir/%s.fastq" % run
                                    for run in runs])
                o3 = out3.replace(".gz", "")
                statement += "; cat %(catlist)s > %(o3)s; \
                gzip -f %(o3)s" % locals()
                statement += "; rm -rf  %(catlist)s" % locals()
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)
            else:
                nspots = "-X %i" % (nreads_download)
                run = pref
                if sra == "ncbi":
                    statement = '''fastq-dump %(nspots)s %(sra_opts)s \
                    --gzip -v \
                    --outdir fastqs.dir %(pref)s &>%(log)s''' % locals()
                elif sra == "ebi":
                    nrows = nreads_download * 4
                    address = getEBIAddress(run)
                    statement = """curl -s %(address)s.fastq.gz | \
                                   gunzip - | \
                                   head -n %(nrows)s \
                                   > fastqs.dir/%(run)s.fastq""" % locals()
                elif sra == "aspera_ebi":
                    nrows = nreads_download * 4
                    address = getEBIAddress(run)
                    asp = address.replace("ftp://ftp.sra.ebi.ac.uk", "")
                    tempname = "run_%s.fastq.gz" % run
                    statement = """
                    ascp -QT -l 300m -P33001 -i %(asperadir)s/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:%(asp)s.fastq.gz fastqs.dir;
                    mv fastqs.dir/%(run)s.fastq.gz fastqs.dir/%(tempname)s;\
                    zcat fastqs.dir/%(tempname)s | head -%(nrows)s > fastqs.dir/%(run)s.fastq;
                    rm -rf fastqs.dir/%(tempname)s""" % locals()
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)
                statement = "gzip -f fastqs.dir/%s.fastq" % pref
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)            
            pathlib.Path(out1).touch()
            pathlib.Path(out2).touch()
        if ispaired:
            o1 = gzip.open(out1).readlines(10)
            o2 = gzip.open(out2).readlines(10)
            if len(o1) != 0 and len(o2) != 0:
                break
            else:
                time.sleep(20)
        else:
            o = gzip.open(out3).readlines(10)
            if len(o) != 0:
                break
            else:
                time.sleep(20)
        x += 1



def runFastQC(infiles, outfiles, pref, ispaired, log, syst=""):

    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        out1 = outfiles[0].replace(".html", ".zip")
        out2 = outfiles[1].replace(".html", ".zip")
        statement = '''fastqc %(in1)s %(in2)s \
        -o fastqc.dir &>%(log)s -t 4; unzip %(out1)s -d fastqc.dir;\
        unzip %(out2)s -d fastqc.dir''' % locals()
        pathlib.Path(outfiles[2]).touch()
    else:
        in1 = infiles[2]
        out1 = outfiles[2].replace(".html", ".zip")
        statement = '''fastqc %(in1)s -o fastqc.dir -t 4 &>%(log)s;\
        unzip %(out1)s -d fastqc.dir''' % locals()
        pathlib.Path(outfiles[0]).touch()
        pathlib.Path(outfiles[1]).touch()

    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)


def trimReadsTrimGalore(infiles, outfiles, ispaired, pref, log, opts,
                        syst=""):
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        statement = '''trim_galore --paired %(in1)s %(in2)s \
        -o trimmed.dir %(opts)s &>%(log)s''' % locals()
        pathlib.Path(outfiles[2]).touch()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        stem1 = in1.split("/")[-1]
        shutil.move("trimmed.dir/%s_trimming_report.txt" % stem1,
                    "logs.dir/%s_trimming_report.txt" % stem1)
        stem2 = in2.split("/")[-1]
        shutil.move("trimmed.dir/%s_trimming_report.txt" % stem2,
                    "logs.dir/%s_trimming_report.txt" % stem2)
    else:
        in1 = infiles[2]
        statement = '''trim_galore %(opts)s %(in1)s \
        -o trimmed.dir &>%(log)s''' % locals()
        pathlib.Path(outfiles[0]).touch()
        pathlib.Path(outfiles[1]).touch()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        stem3 = in1.split("/")[-1]
        shutil.move("trimmed.dir/%s_trimming_report.txt" % stem3,
                    "logs.dir/%s_trimming_report.txt" % stem3)    


def renameReads(infiles, outfiles, ispaired, pref, syst=""):
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        if in1.endswith(".gz"):
            out1 = outfiles[0].replace(".gz", "")
            out2 = outfiles[1].replace(".gz", "")
            cat = 'zcat'
        else:
            out1 = outfiles[0]
            out2 = outfiles[1]
            cat = 'cat'
        statement = '''%(cat)s %(in1)s \
                      | \
                      awk '{ if (NR%%4==1) { print $1"_"$2"/1" }
                      else { print } }'  > %(out1)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        statement = '''%(cat)s %(in2)s \
                      | \
                      awk '{ if (NR%%4==1) { print $1"_"$2"/2" } \
                      else { print } }'  > %(out2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        if in1.endswith(".gz"):
            statement = "gzip %(out1)s; gzip %(out2)s" % locals()
            Run.systemRun(statement, syst)          
        pathlib.Path(outfiles[2]).touch()
    else:
        in1 = infiles[2]
        if in1.endswith(".gz"):
            out1 = outfiles[2].replace(".gz", "")
            cat = 'zcat'
        else:
            out1 = outfiles[2]
            cat = 'cat'
        statement = '''%(cat)s %(in1)s \
                      | \
                      awk '{ if (NR%%4==1) { print $1"_"$2"/1" } \
                      else { print } }'  > %(out1)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        if in1.endswith(".gz"):
            statement = "gzip %(out1)s" % locals()
            Run.systemRun(statement, syst)
        pathlib.Path(outfiles[0]).touch()    
        pathlib.Path(outfiles[1]).touch()    


def dustReadsSGA(infiles, outfiles, pref, ispaired, opts, dt,
                 log, syst=""):
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        out1 = outfiles[0].replace(".gz", "")
        out2 = outfiles[1].replace(".gz", "")
        statement = '''sga preprocess -m 20 --pe-mode 1 --dust \
        --dust-threshold %(dt)s %(opts)s %(in1)s %(in2)s 2>%(log)s |\
        paste  - - - - - - - - | tee >(cut -f 1-4 | tr '\\t' '\\n' \
         > %(out1)s)|\
        cut -f 5-8 | tr '\\t' '\\n' > %(out2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        pathlib.Path(outfiles[2]).touch()
        Run.systemRun(["bash", "-c", statement], syst)
        os.system("gzip %s" % out1)
        os.system("gzip %s" % out2)

    else:
        in1 = infiles[2]
        out1 = outfiles[2].replace(".gz", "")
        statement = '''sga preprocess -m 20 --dust  \
        --dust-threshold %(dt)s %(opts)s %(in1)s 2>%(log)s\
        > %(out1)s''' % locals()
        pathlib.Path(outfiles[0]).touch()
        pathlib.Path(outfiles[1]).touch()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        statement = "gzip %s" % out1
        Run.systemRun(statement, syst)

    ut_functions.writeTime("Dust", "end", pref)
    
    
def subsetReads(infiles, outfiles, pref, ispaired, subsetlen,
                syst=""):
    in1, in2, in3 = infiles
    out1, out2, out3 = [o.replace(".gz", "") for o in outfiles]
    if ispaired:
        statement = '''zcat %(in1)s  | head -%(subsetlen)s > %(out1)s;\
        zcat %(in2)s  | head -%(subsetlen)s > %(out2)s; gzip %(out1)s;\
        gzip %(out2)s;\
        touch %(out3)s.gz''' % locals()
    else:
        statement = '''zcat %(in3)s | head -%(subsetlen)s  > %(out3)s;\
        gzip %(out3)s;\
        touch %(out1)s.gz;\
        touch %(out2)s.gz;''' % locals()
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)


def sampleFastq(infiles, ispaired, nreads, pref, outfile, syst=""):
    in1, in2, in3 = infiles
    temp = "%s_temp.out" % pref
    if ispaired:
        nreads = nreads * 2
        statement = """zcat %(in1)s %(in2)s \
        | awk 'NR%%4==2' | shuf -n %(nreads)s > %(temp)s""" % locals()
    else:
        statement = """zcat %(in3)s \
        | awk 'NR%%4==2' | shuf -n %(nreads)s > %(temp)s""" % locals()
    Run.systemRun(statement)
    lines = [line.strip() for line in open(temp).readlines()]
    os.unlink(temp)
    i = 0
    out = open(outfile, "w")
    for line in lines:
        out.write(">%s\n%s\n" % (i, line))
        i += 1
    out.close()


def bowtieToFastq(infile, outfile, pref, syst=""):
    statement = r"""awk -F "\t" '{printf("@%%s/1\n%%s\n+\n%%s\n", $1, $5, $6)}' %s > %s""" % (infile, outfile)
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)
