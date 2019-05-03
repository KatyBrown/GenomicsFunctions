'''
Functions to run different quantification software within a pipeline using
user defined parameters.
'''


from ruffus import *
import ruffus.cmdline as cmdline
import os
import pathlib
import ut_functions
import Run

def runSalmon(infiles, outfile, ispaired, pref,
              transpath, salmon_opts="", logfile="salmon.log",
              syst=""):
    '''
    Generates a statement and runs Salmon
    
    https://salmon.readthedocs.io/en/latest/salmon.html
    
    Quantifies reads across a transcriptome.
    
    The transcriptome needs to be indexed prior to processing.
    
    Parameters
    ----------
    infiles: list
        list of three input files [mate1, mate2, singleend]
    outfile: str
        path to output file, this is a dummy file as salmon generates multiple
        output files
    ispaired: bool
        True if paired end else False        
    pref: str
        filename prefix (used for the log files)
    transpath: str
        path to the directory containing the indexed transcriptome
        you want to map the reads to
    salmon_opts: str
        additional options to pass directly to salmon
    logfile: str
        path to log file
    syst: str
        system option to run statements on different systems       
    '''
    # if there is a transcriptome (when running as part of the
    # sra_to_filtered_fastq pipeline, in some cases there is no
    # transcriptome to quantify over).
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