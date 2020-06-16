'''
Functions to perform various preprocessing tasks on sequencing data
prior to mapping.

In order to use the same pipeline for single and paired end data, for
software which requires different settings for single and paired
end data a list of three input files and three output files is used, with
blank placeholder files used for the other data type

e.g.
infiles = ['SRR1234_1.fastq.gz', 'SRR1234_1.fastq.gz', 'SRRR1234.fastq.gz']
outfiles = ['SRR1234_1.fastqc', 'SRR1234_2.fastqc', 'SRRR1234.fastqc']

for paired end SRRR1234.fastq.gz is a blank file and the pipeline will generate
a blank file SRRR1234.fastqc

for single end SRR1234_1.fastq.gz and SRR1234_2.fastq.gz are blank and the
pipeline generates blank files SRR1234_1.fastqc and SRR1234_2.fastqc

To run the functions outside of the pipeline any string can replace the blank
filenames e.g.

infiles = ['x', 'x', 'SRRR1234.fastq.gz']
outfiles = ['x', 'x', 'SRRR1234.fastqc']

would run on the single end data and generate a blank file "x"


'''

from ruffus import *
import ruffus.cmdline as cmdline
import os
import sys
import pathlib
import shutil
import ut_functions
import Run
import gzip
import time


def getEBIAddress(run):
    '''
    Formats an SRA run ID into the EBI server address for the same ID
    
    e.g.
    SRR2904018 >
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR290/008/SRR2904018/SRR2904018.fastq.gz
    
    Parameters
    ----------
    run: str
        SRA run ID]
    
    Returns
    -------
    str
        path to SRA ID on EBI server
    '''
    prefix = run[0:6]
    address = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s" % prefix
    if len(run) == 9:
        address = "%s/%s/%s" % (address, run, run)
    elif len(run) == 10:
        suffix = "00%s" % run[-1]
        address = "%s/%s/%s/%s"  % (address, suffix, run, run)
    else:
        # The SRA ID should be 9-10 characters long
        raise ValueError ("run ID not recognised")
    return (address)
    
    
def getSRA(pref, ispaired, log, sra_opts,
           outfiles, nreads_download=1000000000,
           syst="", sra='ncbi', asperadir=""):
    '''
    Downloads dataset "pref" from SRA.
    
    Generates a statement which downloads all the replicates for a sample,
    takes the first n reads from each and concatenates these into
    one output file, then removes the original files and zips the concatenated
    file.
    
    Example statement (paired end, aspera,
    two replicates - SRR2075829 and SRR2904396)
        # download mate 1 for SRR2075829
        ascp -QT -l 300m -P33001 -i /home/kab84/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR207/009/SRR2075829/SRR2075829_1.fastq.gz fastqs.dir;
        # move it into a temporary location
        mv fastqs.dir/SRR2075829_1.fastq.gz fastqs.dir/run_SRR2075829_1.fastq.gz;
        # take the first 50,000,000 reads
        zcat fastqs.dir/run_SRR2075829_1.fastq.gz | head -200000000 > fastqs.dir/SRR2075829_1.fastq;
        # download mate 2 for SRR2075829
        ascp -QT -l 300m -P33001 -i /home/kab84/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR207/009/SRR2075829/SRR2075829_2.fastq.gz fastqs.dir;
        mv fastqs.dir/SRR2075829_2.fastq.gz fastqs.dir/run_SRR2075829_2.fastq.gz;
        zcat fastqs.dir/run_SRR2075829_1.fastq.gz | head -200000000 > fastqs.dir/SRR2075829_2.fastq;
        rm -rf fastqs.dir/run_SRR2075829_1.fastq.gz;
        rm -rf fastqs.dir/run_SRR2075829_2.fastq.gz; 
        # download mate 1 for SRR2075829
        ascp -QT -l 300m -P33001 -i /home/kab84/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR290/006/SRR2904396/SRR2904396_1.fastq.gz fastqs.dir;
        mv fastqs.dir/SRR2904396_1.fastq.gz fastqs.dir/run_SRR2904396_1.fastq.gz;
        zcat fastqs.dir/run_SRR2904396_1.fastq.gz | head -200000000 > fastqs.dir/SRR2904396_1.fastq;
        # download mate 2 for SRR2075829
        ascp -QT -l 300m -P33001 -i /home/kab84/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR290/006/SRR2904396/SRR2904396_2.fastq.gz fastqs.dir;
        mv fastqs.dir/SRR2904396_2.fastq.gz fastqs.dir/run_SRR2904396_2.fastq.gz;
        zcat fastqs.dir/run_SRR2904396_1.fastq.gz | head -200000000 > fastqs.dir/SRR2904396_2.fastq;
        # delete temporary files
        rm -rf fastqs.dir/run_SRR2904396_1.fastq.gz;
        rm -rf fastqs.dir/run_SRR2904396_2.fastq.gz;
        # concatenate SRR2904396 mate 1 and SRR2075829 mate 1 and zip
        cat fastqs.dir/SRR2075829_1.fastq fastqs.dir/SRR2904396_1.fastq > fastqs.dir/SRR2075829.SRR2904396_1.fastq;
        gzip -f fastqs.dir/SRR2075829.SRR2904396_1.fastq;
        # concatenate SRR2904396 mate 2 and SRR2075829 mate 2 and zip
        cat fastqs.dir/SRR2075829_2.fastq fastqs.dir/SRR2904396_2.fastq > fastqs.dir/SRR2075829.SRR2904396_2.fastq;
        gzip -f fastqs.dir/SRR2075829.SRR2904396_2.fastq;
        # remove temporary files
        rm -rf  fastqs.dir/SRR2075829_1.fastq fastqs.dir/SRR2904396_1.fastq;
        rm -rf fastqs.dir/SRR2075829_2.fastq fastqs.dir/SRR2904396_2.fastq;
        
    Parameters
    ----------
    pref: str
        Can be either a single SRA ID or several seperated by "."
    ispaired: bool
        True for paired, False for single ended.
    log: str
        path to a log file
    sra_opts: str
        any additional string of options to pass to the download software
        (fastq-dump or aspera)
    outfiles: list
        list of paths to three outfiles - paired end one, paired end two
        and single end.
    nreads_download: int
        maximum number of reads to download.
    syst: str
        system option to run statements on different systems
    sra: str
        method to use to download the files
        Currently implemented:
        sra - NCBI fastq-dump function
        https://github.com/ncbi/sra-tools
        ebi - EBI ftp server with curl
        aspera_ebi - EBI ftp server with aspera - preferred
    asperadir: str
        https://downloads.asperasoft.com/connect2/
        as described here: https://www.ebi.ac.uk/ega/about/ftp-aspera
        if sra is aspera_ebi, path to the aspera connect directory
    '''
    out1, out2, out3 = outfiles
    x = 0
    # There are very often connectivity problems, so download is attempted
    # up to 25 times
    while True and x < 25:
        # generate statements for paired end data
        if ispaired:
            # merge all the runs into one big fasta file
            if "." in pref:
                # split for multiple replicates with same sample
                runs = pref.split(".")
                statements = []
                # initiate the log file
                l = open(log, "w")
                l.close()

                for run in runs:
                    if sra == "ncbi":
                        # runs ncbi fastq-dump
                        statement = '''fastq-dump %(nreads_download)s %(sra_opts)s \
                        --split-files \
                        --outdir fastqs.dir &>>%(log)s\
                        %(run)s ''' % locals()
                        statements.append(statement)
                    elif sra == "ebi":
                        # fastq - 4 lines per read
                        nrows = nreads_download * 4
                        address = getEBIAddress(run)
                        # downloads from EBI ftp site using curl
                        # -s - silent
                        # unzips the file and takes the first nrows rows
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
                        # fastq - 4 lines per read
                        nrows = nreads_download * 4
                        address = getEBIAddress(run)
                        # converts EBI ftp addres to aspera address
                        asp = address.replace("ftp://ftp.sra.ebi.ac.uk", "")
                        tempname1 = "run_%s_1.fastq.gz" % run
                        tempname2 = "run_%s_2.fastq.gz" % run
                        # download the file using aspera
                        # move to a temporary file and take the first nrows
                        # rows and put these in the output file
                        # remove the temporary files
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

                # concatenate the statements generated above
                statement = "; ".join(statements)
                # make a list containing all the fastq files which will be output
                # and need to be concatenated
                catlist1 = " ".join(["fastqs.dir/%s_1.fastq" % run
                                     for run in runs])
                catlist2 = " ".join(["fastqs.dir/%s_2.fastq" % run
                                     for run in runs])
    
                # filenames for unzipped concatenated output files
                o1 = out1.replace(".gz", "")
                o2 = out2.replace(".gz", "")
                # add to the statement to concatenate the files and zip the output
                statement += "; cat %(catlist1)s > %(o1)s; gzip -f %(o1)s" % locals()
                statement += "; cat %(catlist2)s > %(o2)s; gzip -f %(o2)s" % locals()
                # remove the files once they have been concatenated
                statement += "; rm -rf  %(catlist1)s" % locals()
                statement += "; rm -rf %(catlist2)s" % locals()
                ut_functions.writeCommand(statement, pref)
                # run the statement
                Run.systemRun(statement, syst)

            else:
                # if there is only one run for this sample there is no need
                # to concatenate etc.
                run = pref
                nspots = "-X %i" % (nreads_download)
                if sra == "ncbi":
                    # runs ncbi fastq-dump
                    statement = '''fastq-dump %(nspots)s %(sra_opts)s \
                    --split-files --gzip \
                    --outdir fastqs.dir &>>%(log)s\
                    %(pref)s ''' % locals()
                elif sra == "ebi":
                    # downloads from EBI ftp site using curl
                    # -s - silent
                    # unzips the file and takes the first nrows rows
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
                    # download the file using aspera
                    # move to a temporary file and take the first nrows
                    # rows and put these in the output file
                    # remove the temporary files
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
                    zcat fastqs.dir/%(tempname2)s | head -%(nrows)s  > fastqs.dir/%(run)s_2.fastq;
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
            # as above but for single ended data
            if "." in pref:
                # if are multiple runs for a sample
                runs = pref.split(".")
                statements = []
                l = open(log, "w")
                l.close()
                nspots = "-X %i" % (nreads_download / len(runs))
                for run in runs:
                    if sra == "ncbi":
                        # runs ncbi fastq-dump
                        statement = '''fastq-dump %(nspots)s %(sra_opts)s \
                        --outdir fastqs.dir &>>%(log)s\
                        %(run)s ''' % locals()
                        statements.append(statement)
                    elif sra == "ebi":
                        # downloads from EBI ftp site using curl
                        # -s - silent
                        # unzips the file and takes the first nrows rows
                        nrows = nreads_download * 4
                        address = getEBIAddress(run)
                        statement = """curl -s %(address)s.fastq.gz | \
                                       gunzip - | \
                                       head -n %(nrows)s \
                                       > fastqs.dir/%(run)s.fastq""" % locals()
                        statements.append(statement)
                    elif sra == "aspera_ebi":
                        # download the file using aspera
                        # move to a temporary file and take the first nrows
                        # rows and put these in the output file
                        # remove the temporary files
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
                # concatenate the statements generated above
                statement = "; ".join(statements)
                # make a list containing all the fastq files which will be output
                # and need to be concatenated
                catlist = " ".join(["fastqs.dir/%s.fastq" % run
                                    for run in runs])
                # filenames for unzipped concatenated output files
                o3 = out3.replace(".gz", "")
                # add to the statement to concatenate the files and zip the output
                statement += "; cat %(catlist)s > %(o3)s; \
                gzip -f %(o3)s" % locals()
                # remove the files once they have been concatenated
                statement += "; rm -rf  %(catlist)s" % locals()
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)
            else:
                # if there is only one run for this sample there is no need
                # to concatenate etc.
                nspots = "-X %i" % (nreads_download)
                run = pref
                if sra == "ncbi":
                    # runs ncbi fastq-dump
                    statement = '''fastq-dump %(nspots)s %(sra_opts)s \
                    --gzip \
                    --outdir fastqs.dir %(pref)s &>%(log)s''' % locals()
                elif sra == "ebi":
                    # downloads from EBI ftp site using curl
                    # -s - silent
                    # unzips the file and takes the first nrows rows
                    nrows = nreads_download * 4
                    address = getEBIAddress(run)
                    statement = """curl -s %(address)s.fastq.gz | \
                                   gunzip - | \
                                   head -n %(nrows)s \
                                   > fastqs.dir/%(run)s.fastq""" % locals()
                elif sra == "aspera_ebi":
                    # download the file using aspera
                    # move to a temporary file and take the first nrows
                    # rows and put these in the output file
                    # remove the temporary files                    
                    nrows = nreads_download * 4
                    address = getEBIAddress(run)
                    asp = address.replace("ftp://ftp.sra.ebi.ac.uk", "")
                    tempname = "run_%s.fastq.gz" % run
                    statement = """
                    ascp -QT -l 300m -P33001 -i %(asperadir)s/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:%(asp)s.fastq.gz fastqs.dir;
                    mv fastqs.dir/%(run)s.fastq.gz fastqs.dir/%(tempname)s;\
                    zcat fastqs.dir/%(tempname)s | head -%(nrows)s  > fastqs.dir/%(run)s.fastq;
                    rm -rf fastqs.dir/%(tempname)s""" % locals()
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)
                statement = "gzip -f fastqs.dir/%s.fastq" % pref
                ut_functions.writeCommand(statement, pref)
                Run.systemRun(statement, syst)            
            pathlib.Path(out1).touch()
            pathlib.Path(out2).touch()
        if ispaired:
            # check the output exists (the download didn't fail)
            # if it does move on to the next step, else keep trying
            o1 = gzip.open(out1).readlines(10)
            o2 = gzip.open(out2).readlines(10)
            if len(o1) != 0 and len(o2) != 0:
                break
            else:
                time.sleep(20)
        else:
            # check the output exists (the download didn't fail)
            # if it does move on to the next step, else keep trying
            o = gzip.open(out3).readlines(10)
            if len(o) != 0:
                break
            else:
                time.sleep(20)
        x += 1



def runFastQC(infiles, outfiles, pref, ispaired, log, threads=4,
              syst=""):
    '''
    Generates a statement and runs FastQC
    
    Runs FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    Generates a basic quality report for the data.

    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfiles: list
        list of three output files [mate1, mate2, singleend]
    pref: str
        filename prefix (used for the command log file)
    ispaired: bool
        True if paired end else False
    threads: int
        Number of threads with which to run FastQC
    log: str
        path to output FastQC log file
    syst: str
        system option to run statements on different systems 
    '''
    job_threads = threads
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        out1 = outfiles[0].replace(".html", ".zip")
        out2 = outfiles[1].replace(".html", ".zip")
        # statement to run fastqc and unzip the output
        statement = '''fastqc %(in1)s %(in2)s \
        -o fastqc.dir -threads %(threads)i &>%(log)s;
        unzip %(out1)s -d fastqc.dir;
        unzip %(out2)s -d fastqc.dir''' % locals()
        pathlib.Path(outfiles[2]).touch()
    else:
        in1 = infiles[2]
        out1 = outfiles[2].replace(".html", ".zip")
        # statement to run fastqc and unzip the output
        statement = '''fastqc %(in1)s -o fastqc.dir \
        -threads %(threads)i &>%(log)s;\
        unzip %(out1)s -d fastqc.dir''' % locals()
        pathlib.Path(outfiles[0]).touch()
        pathlib.Path(outfiles[1]).touch()

    ut_functions.writeCommand(statement, pref)
    # run the fastqc statement
    Run.systemRun(statement, syst)


def trimReadsTrimGalore(infiles, outfiles, ispaired, pref, log, opts,
                        syst=""):
    '''
    Generates and runs a statement to trim reads using Trim Galore
    
    (https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfiles: list
        list of three output files [mate1, mate2, singleend]
    pref: str
        filename prefix (used for the command log file)
    ispaired: bool
        True if paired end else False
    log: str
        path to output FastQC log file
    opts: str
        string of additional options to incorporate into the statement used
        to run trim galore
    syst: str
        system option to run statements on different systems     
    '''
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        # statement to run trim galore with paired end settings
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
        # statement to run trim galore with single end settings
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
    '''
    Appends /1 and /2 to read names in a fastq file.

    Generates and runs a statement which reads a fastq file and renames
    the reads according to the convention required by various software
    and removes spaces from read names
    Paired end:
         mate one reads end with /1
         mate two reads end with /2
    Single end:
        reads end with /1
        
    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfiles: list
        list of three output files [mate1, mate2, singleend]
    pref: str
        filename prefix (used for the command log file)
    ispaired: bool
        True if paired end else False
    syst: str
        system option to run statements on different systems  
        
    '''
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        # use cat if the input file is not zipped else zcat
        if in1.endswith(".gz"):
            out1 = outfiles[0].replace(".gz", "")
            out2 = outfiles[1].replace(".gz", "")
            cat = 'zcat'
        else:
            out1 = outfiles[0]
            out2 = outfiles[1]
            cat = 'cat'
        # append /1 to mate one read names and replace spaces and tabs with "_"
        # NR%%4 - every 4th line is a read name
        statement = '''%(cat)s %(in1)s \
                      | \
                      awk '{ if (NR%%4==1) { print $1"_"$2"/1" }
                      else { print } }'  > %(out1)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        # append /2 to mate two read names and replace spaces and tabs with "_"
        # NR%%4 - every 4th line is a read name
        statement = '''%(cat)s %(in2)s \
                      | \
                      awk '{ if (NR%%4==1) { print $1"_"$2"/2" } \
                      else { print } }'  > %(out2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        # if the input file was zipped, zip the output file
        if in1.endswith(".gz"):
            statement = "gzip %(out1)s; gzip %(out2)s" % locals()
            Run.systemRun(statement, syst)          
        pathlib.Path(outfiles[2]).touch()
    else:
        # use cat if the input file is not zipped else zcat
        in1 = infiles[2]
        if in1.endswith(".gz"):
            out1 = outfiles[2].replace(".gz", "")
            cat = 'zcat'
        else:
            out1 = outfiles[2]
            cat = 'cat'
        # append /1 to read names and replace spaces and tabs with "_"
        # NR%%4 - every 4th line is a read name
        statement = '''%(cat)s %(in1)s \
                      | \
                      awk '{ if (NR%%4==1) { print $1"_"$2"/1" } \
                      else { print } }'  > %(out1)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        # if the input file was zipped, zip the output file
        if in1.endswith(".gz"):
            statement = "gzip %(out1)s" % locals()
            Run.systemRun(statement, syst)
        pathlib.Path(outfiles[0]).touch()    
        pathlib.Path(outfiles[1]).touch()    


def dustReadsSGA(infiles, outfiles, pref, ispaired, opts, dt,
                 log, syst=""):
    '''
    Generates and runs a statement to run the SGA preprocess "dust" tool
    
    https://vcru.wisc.edu/simonlab/bioinformatics/programs/sga/preprocess.txt
    Removes low complexity reads.

    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfiles: list
        list of three output files [mate1, mate2, singleend]
    pref: str
        filename prefix (used for the command log file)
    ispaired: bool
        True if paired end else False
    dt: int
        Dust threshold - see https://vcru.wisc.edu/simonlab/bioinformatics/programs/sga/preprocess.txt
    log: str
        path to output FastQC log file
    opts: str
        string of additional options to incorporate into the statement used
        to run SGA
    syst: str
        system option to run statements on different systems
    '''
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        out1 = outfiles[0].replace(".gz", "")
        out2 = outfiles[1].replace(".gz", "")
        # generate the statement to run sga preprocess for paired end data
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
        # generate the statement to run sga preprocess for single end data
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
    '''
    Takes the first n reads from a FASTQ file
    
    Generates and runs a statement which takes the first subsetlen reads
    from a FASTQ file.

    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfiles: list
        list of three output files [mate1, mate2, singleend]
    pref: str
        filename prefix (used for the command log file)
    ispaired: bool
        True if paired end else Fals
    subsetlen: int
        How many reads to select
    syst: str
        system option to run statements on different systems

    '''
    in1, in2, in3 = infiles
    # fastq file - 4 lines per read
    subsetlen = subsetlen * 4
    out1, out2, out3 = [o.replace(".gz", "") for o in outfiles]
    if ispaired:
        # paired end - read both files and output the first subsetlen reads
        statement = '''zcat %(in1)s  | head -%(subsetlen)s > %(out1)s;\
        zcat %(in2)s  | head -%(subsetlen)s > %(out2)s; gzip %(out1)s;\
        gzip %(out2)s;\
        touch %(out3)s.gz''' % locals()
    else:
        # single end - read the file and output the first subsetlen reads
        statement = '''zcat %(in3)s | head -%(subsetlen)s  > %(out3)s;\
        gzip %(out3)s;\
        touch %(out1)s.gz;\
        touch %(out2)s.gz;''' % locals()
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)


def sampleFastq(infiles, ispaired, nreads, pref, outfile, syst=""):
    '''
    Takes a random sample of size nreads from a fastq file.
    
    Saves result into a fasta file.
    For paired end reads the sample is taken from both input files combined -
    the number of reads from each mate is not fixed - and twice as many reads
    are selected.

    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    ispaired: bool
        True if paired end else False
    nreads: int
        The number of reads to sample
    pref: str
        filename prefix (used for the command log file)
    outfile: str
        output fasta file
    syst: str
        system option to run statements on different systems
    '''
    in1, in2, in3 = infiles
    temp = "%s_temp.out" % pref
    if ispaired:
        nreads = nreads * 2
        # concatenate the two files, shuffle and select nreads * 2 reads
        # NR%%4==2 - sequence only
        statement = """zcat %(in1)s %(in2)s \
        | awk 'NR%%4==2' | shuf -n %(nreads)s > %(temp)s""" % locals()
    else:
        # stream the file, shuffle and select nreads reads
        # NR%%4==2 - sequence only
        statement = """zcat %(in3)s \
        | awk 'NR%%4==2' | shuf -n %(nreads)s > %(temp)s""" % locals()
    Run.systemRun(statement)
    lines = [line.strip() for line in open(temp).readlines()]
    os.unlink(temp)
    i = 0
    # convert the output to fasta
    out = open(outfile, "w")
    for line in lines:
        out.write(">%s\n%s\n" % (i, line))
        i += 1
    out.close()
