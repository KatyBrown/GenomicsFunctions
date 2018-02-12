import pandas as pd
import numpy as np
import os
import re
import sys
import pathlib
import tempfile
import Run

def checkSums(infile1, infile2, syst=""):
    statement1 = '''md5sum %s''' % infile1
    statement2 = '''md5sum %s''' % infile2
    cs1 = Run.systemPopen(statement1, syst=syst)[0].split(" ")[0]
    cs2 = Run.systemPopen(statement2, syst=syst)[0].split(" ")[0]
    return(cs1 == cs2)
