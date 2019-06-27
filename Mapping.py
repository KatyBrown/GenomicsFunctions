'''
Functions to run different mapping software within a pipeline using
user defined parameters.
'''

from ruffus import *
import ruffus.cmdline as cmdline
import sys
import pathlib
import ut_functions
import Run


def mapReadsHisat(infiles, outfile, genomepath, genomename, strand,
                  ispaired, options, mismatches, pref, unmapped=[],
                  threads=1, double=False, maxmem="3.9G",
                  syst=""):
    '''
    Generates and runs a statement to map reads using hisat2.
    
    https://ccb.jhu.edu/software/hisat2/index.shtml
    The genome index needs to be generate prior to running the pipeline.
    
    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfile: str
        path to output file, this file is for the bam output
    genomepath: str
        path to the directory containing the indexed genome you want to map
        the reads to
    genomename: str
        prefix of the files in the indexed reference genome
    strand: str
        if a stranded library prep is used, the string passed to hisat2 with
        strand information e.g. --rna-strandness ISF.  This can be
        generated from the salmon output using Infer.getStrand()
    ispaired: bool
        True if paired end else False
    options: str
        additional options to pass directly to hisat2
    mismatches: float
        proportion of mismatches to allow when mapping
    pref: str
        filename prefix (used for the log files)
    unmapped: list
        list of three output files for the unmapped reads as [mate1, mate2, singleend]
    threads: int
        number of threads with which to run hisat2
    syst: str
        system option to run statements on different systems    
    '''
    if double:
        hisat = "%s/%s/hisat/%s" % (genomepath, genomename, genomename)
    else:
        hisat = "%s/hisat/%s" % (genomepath, genomename)
    job_threads = threads
    job_memory = maxmem

    # hisat makes a .log and a .met file when logging
    log = "logs.dir/%s_%s_mapping.log" % (pref, genomename)
    met = "logs.dir/%s_%s_mapping.met" % (pref, genomename)

    # In ignore-quals mode, by default, each mismatch give a penalty of 6
    # so for 8 mismatches this would be a minimum score of 48.
    # `--score-min L,0,-0.48` allows 8 mismatches in a 100bp read, 6 in a 75bp read and 4 in a 50bp read and so on.
    m = -1 * (float(mismatches) * 6)
    mstring = '''--ignore-quals --score-min L,0,%s''' % m
    threads = " -p %s " % threads
    
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        if len(unmapped) != 0:
            # if unmapped output files are specified then keep a fastq file of
            # the unmapped reads
            unmapped_path = unmapped[0].replace(".fastq.1.gz", ".fastq.%.gz")
            unmapped_str = " --un-conc-gz %(unmapped_path)s" % locals()
        else:
            # otherwise just discard them
            unmapped_str = ""
        # build the hisat2 statement
        statement = '''hisat2 -x %(hisat)s %(strand)s -1 %(in1)s -2 %(in2)s \
        %(mstring)s --met-file %(met)s %(options)s %(unmapped_str)s %(threads)s 2>%(log)s|\
        samtools view -F4 -b > %(outfile)s''' % locals()
    else:
        in1 = infiles[2]
        if len(unmapped) != 0:
            # if unmapped output files are specified then keep a fastq file of
            # the unmapped reads            
            unmapped_path = unmapped[2]
            unmapped_str = " --un-gz %(unmapped_path)s" % locals()
        else:
            # otherwise just discard them
            unmapped_str = ""
        # build the hisat2 statement
        statement = '''hisat2 -x %(hisat)s %(strand)s -U %(in1)s \
        %(mstring)s --met-file %(met)s %(options)s %(unmapped_str)s %(threads)s 2>%(log)s|\
        samtools view -F4 -b > %(outfile)s''' % locals()
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)
    for u in unmapped:
        pathlib.Path(u).touch()


def mapReadsBowtie(infiles, outfile, genomepath, genomename, strand,
                   ispaired, options, mismatches, pref, unmapped=[],
                   threads=1, maxmem="1.9G", double=False,
                   syst=""):
    '''
    Generates and runs a statement to map reads using bowtie2.
    
    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    The genome index needs to be generate prior to running the pipeline.
    Bowtie2 and Hisat2 mostly have identical parameters but bowtie2 doesn't
    need a "strand" parameter

    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfile: str
        path to output file, this file is for the bam output
    genomepath: str
        path to the directory containing the indexed genome you want to map
        the reads to
    genomename: str
        prefix of the files in the indexed reference genome

    ispaired: bool
        True if paired end else False
    options: str
        additional options to pass directly to hisat2
    mismatches: float
        proportion of mismatches to allow when mapping
    pref: str
        filename prefix (used for the log files)
    unmapped: list
        list of three output files for the unmapped reads as [mate1, mate2,
        singleend]
    threads: int
        number of threads with which to run hisat2
    syst: str
        system option to run statements on different systems    
    '''
    if double is True:
        bowtie = "%s/%s/bowtie/%s" % (genomepath, genomename)
    elif double is "exact":
        bowtie = genomepath
    else:
        bowtie = "%s/bowtie/%s" % (genomepath, genomename)

    job_threads = threads
    job_memory = maxmem
    
    # bowtie makes a .log and a .met file when logging
    log = "logs.dir/%s_%s_mapping.log" % (pref, genomename)
    met = "logs.dir/%s_%s_mapping.met" % (pref, genomename)

    # In ignore-quals mode, by default, each mismatch give a penalty of 6
    # so for 8 mismatches this would be a minimum score of 48.
    # `--score-min L,0,-0.48` allows 8 mismatches in a 100bp read, 6 in a 75bp read and 4 in a 50bp read and so on.
    m = -1 * (float(mismatches) * 6)
    mstring = '''--ignore-quals --score-min L,0,%s''' % m
    threads = " -p %s " % threads
    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        if len(unmapped) != 0:
            # if unmapped output files are specified then keep a fastq file of
            # the unmapped reads
            unmapped_path = unmapped[0].replace(".fastq.1.gz", ".fastq.%.gz")
            unmapped_str = " --un-conc-gz %(unmapped_path)s" % locals()
        else:
            # otherwise just discard them
            unmapped_str = ""
        # build the bowtie2 statement
        statement = '''bowtie2 -x %(bowtie)s -1 %(in1)s -2 %(in2)s \
        --met-file %(met)s %(options)s %(unmapped_str)s %(threads)s %(mstring)s 2>%(log)s|\
        samtools view -F4 -b > %(outfile)s''' % locals()
    else:
        in1 = infiles[2]
        if len(unmapped) != 0:
            # if unmapped output files are specified then keep a fastq file of
            # the unmapped reads
            unmapped_path = unmapped[2]
            unmapped_str = " --un-gz %(unmapped_path)s" % locals()
        else:
            # otherwise just discard them
            unmapped_str = ""
        # build the bowtie2 statement
        statement = '''bowtie2 -x %(bowtie)s -U %(in1)s \
        --met-file %(met)s %(options)s %(unmapped_str)s %(threads)s %(mstring)s 2>%(log)s|\
        samtools view -F4 -b > %(outfile)s''' % locals()
        
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)
    for u in unmapped:
        pathlib.Path(u).touch()


def mapReadsBowtie1(infiles, outfile, genomepath, genomename,
                    ispaired, options, mismatches, pref, threads, double=False,
                    maxmem="1.9G",
                    syst=""):
    '''
    Generates and runs a statement to map reads using bowtie1
    
    http://bowtie-bio.sourceforge.net/index.shtml
    
    The genome index needs to be generated prior to running the pipeline.
    
    Keeping unmapped reads is not yet implemented.
    
    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfile: str
        path to output file, this file is for the bam output
    genomepath: str
        path to the directory containing the indexed genome you want to map
        the reads to
    genomename: str
        prefix of the files in the indexed reference genome
    ispaired: bool
        True if paired end else False
    options: str
        additional options to pass directly to hisat2
    mismatches: int
        number of mismatches to allow when mapping
    pref: str
        filename prefix (used for the log files)
    threads: int
        number of threads to run bowtie with
    syst: str
        system option to run statements on different systems    
    '''
    if double is True:
        bowtie = "%s/%s/bowtie1/%s"
    elif double is "exact":
        bowtie = genomepath
    else:
        bowtie = "%s/bowtie1/%s" % (genomepath, genomename)
    job_threads = threads
    job_memory = maxmem

    log = "logs.dir/%s_%s_mapping.log" % (pref, genomename)
    met = "logs.dir/%s_%s_mapping.met" % (pref, genomename)

    if ispaired:
        in1 = infiles[0]
        in2 = infiles[1]
        statement = '''bowtie -S -q %(bowtie)s \
        -1 %(in1)s -2 %(in2)s \
         %(options)s -v %(mismatches)s  2>%(log)s|\
        samtools view -F4 -b > %(outfile)s''' % locals()
    else:
        in1 = infiles[2]
        statement = '''bowtie -S -q %(bowtie)s %(in1)s \
        %(options)s -v %(mismatches)s 2>%(log)s|\
        samtools view -F4 -b > %(outfile)s''' % locals()
    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)


def filterMappedReads(infile, outfiles, ispaired, genomename, pref,
                      q=10, maxmem='4G', syst=""):
    '''
    Extracts unmapped reads from a bam file and converts to FASTQ

    Filters a bam file to remove mapped, properly paired reads with 
    quality scores above q then converts the remainder back into a fastq file
    (the unmapped reads).

    Parameters
    ----------
    infiles: str
        path to the bam file to be filtered
    outfiles: list
        list of three output fastq files [mate1, mate2, unpaired]
    ispaired: bool
        True if paired end else False
    genomename: str
        Name of the genome which was mapped to, only used to name the log file
    pref: str
        filename prefix (used for the log files)
    q: int
        minimum quality score to class a read as "mapped"
    syst: str
        system option to run statements on different systems   
    '''
    job_memory = maxmem
    out1, out2, out3 = outfiles
    tempout = infile.replace(".bam", ".temp")
    tempout2 = infile.replace(".bam", ".temp2")

    log = "logs.dir/%s_%s_filtering.log" % (pref, genomename)

    if ispaired:
        # Filter the bam file using 
        # -F4 (mapped) -q X (minimum quality) -f1 -f2
        # (properly paired) and send the reads passing the filter
        # to tempout2 and the reads failing the filter to tempout
        statement = '''samtools view -F4 -q %(q)s -f1 -f2 -b \
                       -U %(tempout)s  %(infile)s \
                       > %(tempout2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        # sort the unmapped reads and convert to fastq
        statement = '''samtools sort -n %(tempout)s |\
                       samtools fastq -i - -1 %(out1)s -2 %(out2)s \
                       &>%(log)s;\
                       rm -rf %(tempout)s %(tempout2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        pathlib.Path(out3).touch()

    else:
        tempout = infile.replace(".bam", ".temp")
        # Filter the bam file using 
        # -F4 (mapped) -q X (minimum quality) and
        # send the reads passing the filter
        # to tempout2 and the reads failing the filter to tempout
        statement = '''samtools view -F4 -q%(q)s -b\
                       -U %(tempout)s %(infile)s \
                       > %(tempout2)s''' % locals()
        ut_functions.writeCommand(statement, pref)
        Run.systemRun(statement, syst)
        # sort the unmapped reads and convert to fastq
        statement = '''samtools sort -n %(tempout)s | \
                       samtools fastq -i - -0 %(out3)s &>%(log)s;
                       rm -rf %(tempout)s %(tempout2)s''' % locals()
        Run.systemRun(statement, syst)
        pathlib.Path(out1).touch()
        pathlib.Path(out2).touch()


def cleanBam(infile, outfile, ispaired, pref,
             q=10, syst=""):
    '''
    Filters a bam file.
    
    Filters a bam file to keep only mapped, properly paired reads with a 
    quality score >= q
    
    Parameters
    ----------
    infile: str
        path to input bam file
    outfile: str
        path to output bam file
    ispaired: bool
        True if paired end else False
    pref: str
        filename prefix (used for the log files)
    q: int
        minimum acceptable quality score
    '''
    if ispaired == "paired":
        statement = '''samtools view -b -F4 -q%(q)s -f1 -f2 %(infile)s >\
                       %(outfile)s''' % locals()
    else:
        statement = '''samtools view -b -F4 -q%(q)s %(infile)s >\
                       %(outfile)s''' % locals()

    ut_functions.writeCommand(statement, pref)
    Run.systemRun(statement, syst)
