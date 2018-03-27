import os
import pandas as pd
import numpy as np
import glob
import xmltodict
import ut_functions
import Run

def findTaxonID(species_name, syst=""):
    '''
    Return the taxonomy ID for a species name
    '''
    statement = 'esearch -db taxonomy -query "%s" \
                 | efetch' % species_name
    L = Run.systemPopen(statement, syst)
    if len(L) == 0:
        return "000"
    else:
        return L[0]


def findTaxonName(taxid, syst=""):
    '''
    Return the scientific name for a taxonomy ID
    '''
    statement = '''efetch -db taxonomy -id %s -format xml | \
                    xtract -pattern ScientificName \
                    -element ScientificName''' % taxid
    L = Run.systemPopen(statement, syst)
    return L[0]


def findFamilyGenus(taxonid, nodespath, syst=""):
    '''
    Return the family and genus for a taxonomy ID,
    based on the nodes.dmp file downloaded from NCBI
    and stored in nodespath.
    '''
    tab = pd.read_csv(nodespath, sep="|", header=None)
    tab[2] = tab[2].str.replace("\t", "")
    x = tab[tab[0] == int(taxonid)]
    types = set(x[2].values)
    newvals = set(x[1])
    j = 0
    k = 0
    while "family" not in types and j < 15:
        x2 = tab[tab[0].isin(newvals)]
        types = set(x2[2].values)
        newvals = set(x2[1])
        j += 1
    if j == 15:
        fam = "NA"
    else:
        fam = (x2[0].values[0])
    types = set(x[2].values)
    newvals = set(x[1])
    while "genus" not in types and k < 15:
        x2 = tab[tab[0].isin(newvals)]
        types = set(x2[2].values)
        newvals = set(x2[1])
        k += 1
    if k == 15:
        genus = "NA"
    else:
        genus = (x2[0].values[0])
    if taxonid == "2065263":
        genus = "6726"
        fam = "6725"
    return (fam, genus, findTaxonName(fam, syst=syst),
            findTaxonName(genus, syst=syst))

    
def getMetadataSRA(ID, genomesdir, nodespath, outfile,
                   log="on", syst=""):
    genomes = glob.glob("%s/*/" % genomesdir)
    genome_names = [line.split("/")[-2] for line in genomes]

    genometab = [line.strip().split("\t")
                 for line in open("%s/table.tsv" % genomesdir)]
    family = []
    genus = []
    species = []
    speciesnames = []
    for g in genometab:
        speciesnames.append(g[0])
        species.append(g[1])
        genus.append(g[3])
        family.append(g[5])
    fdict = dict(zip(family, speciesnames))
    gdict = dict(zip(genus, speciesnames))

    x = 0
    while True:
        # Retrieve the full XML record for the ID
        # this query often seems to fail for no reason - keep trying
        # until it doesn't or for 10 tries
        try:
            statement = 'efetch -db sra -id %s -format xml' % ID
            if log == "on":
                 ut_functions.writeCommand(statement, ID)
            y = "\n".join(Run.systemPopen(statement, syst))
            xdict = xmltodict.parse(y)
            break
        except:
            if x == 10:
                print ("Failed to connect to NCBI: attempt %i" % x)
                break
            else:
                x += 1
                continue
    # Parse this into a dictionary
    d = xdict['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
    experiment = d['EXPERIMENT']
    expid = experiment['@accession']
    sample = d['SAMPLE']
    libtype = experiment['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY']
    pe = experiment['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_LAYOUT']
    taxonid = sample['SAMPLE_NAME']['TAXON_ID']
    sciname = sample['SAMPLE_NAME']['SCIENTIFIC_NAME']
    biosample = sample['@accession']
    title = experiment['TITLE']
    spots = d['Pool']['Member']['@spots']

    runs = d['RUN_SET']['RUN']
    if '@accession' in runs:
        run = d['RUN_SET']['RUN']['@accession']
    else:
        run = []
        runsize = 0
        for res in runs:
            if '@accession' in res:
                run.append(res['@accession'])
        run = ".".join(run)
    if '@size' in runs:
        runsize = int(d['RUN_SET']['RUN']['@size'])
    else:
        runsize = 0
        if '@size' in res:
            runsize += int(res['@size'])
    r = ((expid, biosample, libtype, taxonid, sciname,
          title, run, runsize, spots))
    r = [r]
    results = pd.DataFrame(r)
    results.columns = ['Experiment', 'Biosample', 'Type', 'TaxonID',
                       'SciName', 'Title', 'Runs', 'Runsize', 'Spots']
    results['Runsize_h'] = (
        results['Runsize'] / 1e9).round(2).astype(str) + "G"
    results['SciName'] = results['SciName'].str.replace(" ", "_")

    familyID, genusID, family, genus = findFamilyGenus(taxonid,
                                                       nodespath)


    genusID = str(genusID)
    familyID = str(familyID)

    results['FamilyID'] = familyID
    results['GenusID'] = genusID
    results['Family'] = family
    results['Genus'] = genus
    statement = '''esearch -db assembly -query "txid%s[Organism]" | \
                 efetch -format uilist''' % taxonid
    if log == "on":
        ut_functions.writeCommand(statement, ID)
    x = Run.systemPopen(statement, syst)
    print (x)
    if len(x) >= 1:
        reftype = "H"
    else:
        statement = '''esearch -db assembly -query "txid%s[Organism]" \
                      | efetch -format uilist''' % genusID
        if log == "on":
            ut_functions.writeCommand(statement, ID)
        x = Run.systemPopen(statement, syst)
        if len(x) >= 1 and not genus.ID.endswith("dae"):
            reftype = "G"
        else:
            if genusID.endswith("dae"):
                familyID = genusID
            statement = '''esearch -db assembly \
                          -query "txid%s[Organism]" \
                          | efetch -format uilist''' % familyID
            if log == "on":
                ut_functions.writeCommand(statement, ID)
            x = Run.systemPopen(statement, syst)
            if len(x) >= 1:
                reftype = "F"
            else:
                reftype = "N"
                reference = "X"
    results['Reference_Type'] = reftype

    if reftype == "H":
        if sciname == "Drosophila pseudoobscura":
            sciname = "Drosophila pseudoobscura pseudoobscura"
        if sciname == "Apis":
            sciname = "Apis mellifera"
        if sciname == "Heliconius melpomene":
            sciname = "Heliconius melpomene melpomene"
        if sciname == "Heliconius cydno":
            sciname = "Heliconius cydno hermogenes"

        assert taxonid in species, "Host %s reference genome exists but not found" % sciname
        reference = sciname
    elif reftype == "G":
        assert genusID in gdict, "Genus %s reference genome exists but not found"  % genus
        reference = gdict[genusID]
        
    elif reftype == "F":
        assert familyID in fdict, "Family %s reference genome exists but not found" % family
        reference = fdict[familyID]
    results['Reference'] = reference
    results.to_csv(outfile, sep="\t", index=None)
