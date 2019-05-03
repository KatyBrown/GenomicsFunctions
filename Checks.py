import Run

def checkSums(infile1, infile2, syst=""):
    '''
    Checks that two files are identical using checksums
    
    Parameters
    ----------
    infile1: str
        path to file 1
    infile2: str
        path to file 2
    syst: str
        system option to run statements on different systems
    
    Returns
    -------
    bool
        True if files are identical else False
    '''
    statement1 = '''md5sum %s''' % infile1
    statement2 = '''md5sum %s''' % infile2
    cs1 = Run.systemPopen(statement1, syst=syst)[0].split(" ")[0]
    cs2 = Run.systemPopen(statement2, syst=syst)[0].split(" ")[0]
    return(cs1 == cs2)
