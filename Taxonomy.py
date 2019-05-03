import os
import pandas as pd
import numpy as np
import glob
import xmltodict
import ut_functions
import Run
import time

def findTaxonID(species_name, syst=""):
    '''
    Return the taxonomy ID for a species name
    '''
    i = 0
    while True:
        statement = 'esearch -db taxonomy -query "%s" \
        | efetch' % species_name
        L = Run.systemPopen(statement, syst)
        if len(L) == 0 and i > 50:
            return "000"
        elif len(L) > 0:
            return L[0]
        i += 1

def findTaxonName(taxid, syst=""):
    '''
    Return the scientific name for a taxonomy ID
    '''
    i = 0
    while True:
        statement = '''efetch -db taxonomy -id %s -format xml | \
        xtract -pattern ScientificName \
        -element ScientificName''' % taxid
        L = Run.systemPopen(statement, syst)
        if len(L) == 0 and i > 50:
            return "000"
        elif len(L) != 0:
            return L[0]
        i += 1


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
        fam = "000"
    else:
        try:
            fam = (x2[0].values[0])
        except:
            fam = "000"
    types = set(x[2].values)
    newvals = set(x[1])
    while "genus" not in types and k < 15:
        x2 = tab[tab[0].isin(newvals)]
        types = set(x2[2].values)
        newvals = set(x2[1])
        k += 1
    if k == 15:
        genus = "000"
    else:
        try:
            genus = (x2[0].values[0])
        except:
            genus = "000"
    if taxonid == "2065263":
        genus = "6726"
        fam = "6725"
    if genus == "000":
        genusname = "000"
    else:
        genusname = findTaxonName(genus)
    if fam == "000":
        familyname = "000"
    else:
        familyname = findTaxonName(fam)
    return (fam, genus, familyname,
            genusname)


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
    species.append("000")
    fdict = dict(zip(family, speciesnames))
    gdict = dict(zip(genus, speciesnames))

    x = 0
    while True:
        # Retrieve the full XML record for the ID
        # this query often seems to fail for no reason - keep trying
        # until it doesn't or for 10 tries
        s1 = 'efetch -db sra -format xml -id %s' % ID
        s2 = 'esearch -db sra -query %s | efetch -format xml' % ID
        try:
            if x < 5:
                statement = s1
            else:
                statement = s2
            if log == "on":
                 ut_functions.writeCommand(statement, ID)
            y = "\n".join(Run.systemPopen(statement, syst)).replace("</CENTER_NAME-->", "")
            xdict = xmltodict.parse(y)
            break
        except:
            if x == 100:
                print ("Failed to connect to NCBI: attempt %i" % x)
                raise SystemExit
            else:
                time.sleep(10)
                x += 1
                continue
    # Parse this into a dictionary
    d = xdict['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']
    try:
        experiment = d['EXPERIMENT']
    except:
        d = d[0]
        experiment = d['EXPERIMENT']
    expid = experiment['@accession']
    sample = d['SAMPLE']
    libtype = experiment['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY']
    pe = experiment['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_LAYOUT']
    taxonid = sample['SAMPLE_NAME']['TAXON_ID']
    if 'SCIENTIFIC_NAME' in sample['SAMPLE_NAME']:
        sciname = sample['SAMPLE_NAME']['SCIENTIFIC_NAME']
    else:
        taxid = int(sample['SAMPLE_NAME']['TAXON_ID'])
        sciname = findTaxonName(taxid)
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
    i = 0
    while True:
        statement = '''esearch -db assembly -query "txid%s[Organism]" | \
        efetch -format uilist''' % taxonid
        if log == "on":
            ut_functions.writeCommand(statement, ID)
        x = Run.systemPopen(statement, syst)
        if len(x) >= 1:
            reftype = "H"
            break
        elif i > 50:
            reftype = "X"
            break
        i += 1
    if reftype == "X":
        i = 0
        while True:
            statement = '''esearch -db assembly -query "txid%s[Organism]" \
            | efetch -format uilist''' % genusID
            if log == "on":
                ut_functions.writeCommand(statement, ID)
            x = Run.systemPopen(statement, syst)
            if len(x) >= 1 and not genus.endswith("dae"):
                reftype = "G"
                break
            elif len(x) > 1 and genus.endswith("dae"):
                reftype = "X"
                family = genus
                familyID = genusID
                break
            elif i > 50:
                reftype = "X"
                break
            i += 1
        if reftype == "X":
            i = 0
            while True:
                statement = '''esearch -db assembly \
                -query "txid%s[Organism]" \
                | efetch -format uilist''' % familyID
                if log == "on":
                    ut_functions.writeCommand(statement, ID)
                x = Run.systemPopen(statement, syst)
                if len(x) >= 1:
                    reftype = "F"
                    break
                elif i > 50:
                    reftype = "N"
                    reference = "X"
                    break
                i += 1
    results['Reference_Type'] = reftype

    if reftype == "H":
        if sciname == "Drosophila pseudoobscura":
            sciname = "Drosophila pseudoobscura pseudoobscura"
            taxonid = "46245"
        if sciname == "Apis":
            sciname = "Apis mellifera"
            taxonid = "7460"
        if sciname == "Heliconius melpomene":
            sciname = "Heliconius melpomene melpomene"
            taxonid = "171917"
        if sciname == "Heliconius cydno":
            sciname = "Heliconius cydno hermogenes"
            taxonid = "1501340"
        if sciname == "Hawaiian Drosophila":
            sciname = 'Drosophila grimshawi'
            taxonid = '7222'
        if sciname == "Ceratosolen solmsi":
            sciname = "Ceratosolaen solmsi marchali"
            taxonid = '326594'
        if sciname == "Glossina palpalis":
            sciname = "Glossina palpalis gambiensis"
            taxonid = '67801'
        if sciname == "Odonata" or sciname == "Lepidoptera" or sciname == "Insecta" or sciname == "Curculionidae":
            sciname = "X"
            taxonid = "000"
        if sciname == "Galleria":
            sciname = "Galleria_mellonella"
            taxonid = "7137"
        if sciname == "Penaeus":
            sciname = "Penaeus_monodon"
            taxonid = "6687"
        if sciname == "Ixodes scapularis + Ehrlichia sp.":
            sciname = "Ixodes_scapularis"
            taxonid = "6945"
        if sciname == "L. boulardi":
            sciname = "Leptopilina_boulardi"
            taxonid = "63433"
        
        assert taxonid in species, "Host %s reference genome exists but not found" % sciname
        reference = sciname
    elif reftype == "G":
        assert genusID in gdict, "Genus %s reference genome exists but not found"  % genus
        reference = gdict[genusID]
        
    elif reftype == "F":
        if familyID == "6913":
            familyID = "450948"
        assert familyID in fdict, "Family %s reference genome exists but not found" % family
        reference = fdict[familyID]
    results['Reference'] = reference
    results.to_csv(outfile, sep="\t", index=None)
