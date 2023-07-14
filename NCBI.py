import ut_functions
import Bio.Entrez as Entrez
import numpy as np
import math
import xmltodict
import pandas as pd
import os
import shutil
import time


def getApiKey():
    params = ut_functions.readIni("pipeline.ini")
    return (params['api_key'])

def fixStatement(statement, syst=""):
    # statement = statement.replace("esearch ", "esearch -api %s " % getApiKey())
    # statement = statement.replace("efetch ", "efetch -api %s " % getApiKey())
    return (statement)


def getRecord(acc,
              db='nuccore',
              email,
              api,
              silent=True):
    '''
    Retrieve taxonomic lineage for an NCBI taxonomy ID using their API,
    plus the common name.
    '''
#    # NCBI API requires an email address
    Entrez.email = email
    Entrez.api = api
    # Retrieve all data for this taxonomy ID from NCBI Taxonomy
    x = 0
    while x < 10:
        try:
            search = Entrez.efetch(id=acc,
                                   db=db,
                                   retmode="xml")
            if not silent:
                print ("attempt %s" % x)
            break
        except Entrez.HTTPError:
            time.sleep(10)
            x += 1
    if x != 10:
        try:
            record = Entrez.read(search, validate=False)
        except Entrez.Parser.CorruptedXMLError:
            record = dict()
        return(record)
    else:
        print ("Failed for %s" % acc)
        return (dict())


def getRecordSRA(acc,
              email,
              api):
    '''
    Retrieve taxonomic lineage for an NCBI taxonomy ID using their API,
    plus the common name.
    '''
    # NCBI API requires an email address
    Entrez.email = email
    Entrez.api = api
    x = 0
    while x < 5:
        try:
            # Retrieve all data for this taxonomy ID from NCBI Taxonomy
            search = Entrez.efetch(id=acc,
                                   db='sra',
                                   retmode="xml")
            print ("hi %s" % x)
            break
        except Entrez.HTTPError:
            time.sleep(10)
            x += 1
    if x != 10:
        record = xmltodict.parse(search.read())
        record = record['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
        return(record)
    else:
        print ("Failed for %s" % acc)
        return (dict())


def splitaccs(accs, splitsize):
    nchunks = math.ceil(len(accs) / splitsize)
    accarray = np.array(list(accs))
    chunks = np.array_split(accarray, nchunks)
    return(chunks)

def getRecords(accs, chunksize=50, db='nuccore', silent=False):
    chunks = splitaccs(accs, chunksize)
    records = []
    recordD = dict()
    for i, chunk in enumerate(chunks):
        div = 10**math.floor(math.log10(len(chunks)))
        if not silent:
            if i % int(div) == 0:
                print ("Searched %i / %i accession blocks" % (i, len(chunks)))
        records += getRecord(",".join(chunk), db, silent=silent)
    for record in records:
        if 'GBSeq_accession-version' in record:
            recordD[record['GBSeq_accession-version']] = record
        elif 'TaxId' in record:
            recordD[record['TaxId']] = record
        else:
            print (record)
            raise RuntimeError ("ID not found")
    return (recordD)

def getRecordsSRA(accs):
    records = []
    recordD = dict()
    div = 10**math.floor(math.log10(len(accs)))
    for i, acc in enumerate(accs):
        if i % int(div) == 0:
            print ("Searched %i / %i accessions" % (i, len(accs)))
        records.append(getRecordSRA(acc))
    for record in records:

        trueaccs = []
        rr = dict()
        for level in ['experiment', 'study', 'sample']:
            acc = record[level.upper()]['@accession']
            if acc in accs:
                trueaccs.append(acc)
            rr[level] = acc
        runs = record['RUN_SET']
        if not isinstance(runs, list):
            runs = [runs]
        runL = []
        for run in runs:

            acc = run['RUN']
            if not isinstance(acc, list):
                acc = [acc]
            for a in acc:
                a2 = a['@accession']
                if a2 in accs:
                    trueaccs.append(a2)
                runL.append(a)
        rr['runs'] = runL
        rr['full'] = record
        assert len(trueaccs) != 0
        for acc in trueaccs:
            recordD[acc] = rr
    return (recordD)


def cleanNCBITaxonomy(record):
    tD = dict()
    tD['Current'] = (record['TaxId'], record['ScientificName'], record['Rank'])

    if record['ParentTaxId'] == '1':
        tD[record['Rank']] = (record['ScientificName'], record['TaxId'])
    tD[record['Rank']] = (record['ScientificName'], record['TaxId'])
    segs = record['LineageEx']
    
    for seg in segs:
        ID = seg['TaxId']
        nam = seg['ScientificName']
        rank = seg['Rank']
        tD[rank] = (nam, ID)
    return (tD)


def getStandardTaxRanks():
    ranks = ['superkingdom', 'kingdom', 'phylum', 'subphylum', 'class',
             'order', 'family', 'subfamily', 'genus', 'species']
    return (ranks)


def getTaxCols(ranks):
    cols = ['TaxID', 'NCBI_ScientificName', 'NCBI_Rank']
    for rank in ranks:
        cols.append("%s_name" % rank)
        cols.append("%s_TaxID" % rank)
    return (cols)


def makeTaxTab(taxD, ranks=getStandardTaxRanks()):
    rows = []
    for acc, rec in taxD.items():
        assert acc == rec['Current'][0]
        row = [int(acc), rec['Current'][1], rec['Current'][2]]
        for rank in ranks:
            if rank in rec:
                row.append(rec[rank][0])
                row.append(int(rec[rank][1]))
            else:
                row.append("X")
                row.append(0)
        rows.append(row)
    tab = pd.DataFrame(rows)
    cols = getTaxCols(ranks)
    tab.columns = cols
    return (tab)

def makeBlankTaxTab(ranks=getStandardTaxRanks()):
    cols = getTaxCols(ranks)
    df = pd.DataFrame(columns=['Accession'] + cols)
    return (df)

    
def getNCBISeqs(db, acclist, chunksize, outfile):
    div = 10**math.floor(math.log10(len(acclist)))
    chunks = splitaccs(acclist, chunksize)
    API = getApiKey()
    tempdir = "temp_seqs"
    statement = "export NCBI_API_KEY=%s; mkdir -p %s" % (API, tempdir)
    os.system(statement)
    out_comb = open(outfile, "w")
    for i, chunk in enumerate(chunks):
        if i % int(div) == 0:
            print ("Searched %i / %i accessions" % (i, len(acclist)))
        out = open("%s/%i_list.txt" % (tempdir, i), "w")
        for acc in chunk:
            out.write("%s\n" % (acc))
        out.close()
        statement_c = "efetch -input %s/%i_list.txt \
-format fasta -db %s > %s/%s_seqs.fasta; " % (tempdir, i, db, tempdir, i)
        os.system(statement_c)
        lines = open("%s/%s_seqs.fasta" % (tempdir, i)).readlines()
        for line in lines:
            out_comb.write(line)
    out_comb.close()
    shutil.rmtree(tempdir)
    
