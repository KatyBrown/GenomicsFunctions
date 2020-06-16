import pandas as pd
import xmltodict
import ut_functions
import Run
import time
import NCBI

def findTaxonID(species_name, syst=""):
    '''
    Returns the taxonomy ID for a species name

    Uses the NCBI command line toolkit
    https://www.ncbi.nlm.nih.gov/books/NBK179288/

    Parameters
    ----------
    species_name: str
        scientific name of the species
    syst: str
        system option to run statements on different systems
    
    Returns
    -------
    str
        NCBI taxonomy ID for species
    '''

    i = 0
    species_name = species_name.replace("_", " ")
    while True:
        statement = 'esearch -db taxonomy -query "%s" \
        | efetch' % species_name
  #      statement = NCBI.fixStatement(statement)
        L = Run.systemPopen(statement, syst)
        if len(L) == 0 and i > 50:
            return "000"
        elif len(L) > 0:
            return L[0]
        i += 1


def findTaxonName(taxid, syst=""):
    '''
    Return the scientific name for a taxonomy ID.
    
    Uses the NCBI command line toolkit
    https://www.ncbi.nlm.nih.gov/books/NBK179288/
    
    Parameters
    ----------
    taxid: str
        NCBI taxonomy ID of the species
    syst: str
        system option to run statements on different systems
        
    Returns
    -------
    str
        scientific name of taxon
    '''
    i = 0
    while True:
        statement = '''efetch -db taxonomy -id %s -format xml | \
        xtract -pattern ScientificName \
        -element ScientificName''' % taxid
        L = Run.systemPopen(statement, syst)
    #    statement = NCBI.fixStatement(statement)
        if len(L) == 0 and i > 50:
            return "000"
        elif len(L) != 0:
            return L[0]
        i += 1


def makeDictsClassifications(nodespath, namespath):
    tab = pd.read_csv(nodespath, sep="|", header=None)
    tab[2] = tab[2].str.replace("\t", "")
    taxa_tab = pd.read_csv(namespath, sep="\t", header=None)
    taxa_tab[2] = taxa_tab[2].str.replace("\t", "")
    taxa_names_rev = dict(zip(taxa_tab[2], taxa_tab[0]))

    taxa_tab = taxa_tab[taxa_tab[6] == 'scientific name']
    taxa_names = dict(zip(taxa_tab[0], taxa_tab[2]))
    links = dict(zip(tab[0], tab[1]))
    classes = dict(zip(tab[0], tab[2]))
    return (classes, links, taxa_names, taxa_names_rev)


def findClassification(classification, taxonid,
                       classes,
                       links,
                       taxa_names, syst=""):
    allclassifications = ['class',
                          'cohort',
                          'family',
                          'forma',
                          'genus',
                          'infraclass',
                          'infraorder',
                          'kingdom',
                          'order',
                          'parvorder',
                          'phylum',
                          'species',
                          'species group',
                          'species subgroup',
                          'subclass',
                          'subfamily',
                          'subgenus',
                          'subkingdom',
                          'suborder',
                          'subspecies',
                          'subtribe',
                          'superclass',
                          'superfamily',
                          'superorder',
                          'tribe']

    types = set([classes[taxonid]])
    j = 0
    current = taxonid
    full = dict()
    fullid = dict()
    full[classes[taxonid]] = taxa_names[current]
    fullid[classes[taxonid]] = current
    while classification not in types and j < 15:
        nex = links[current]
        nex_class = classes[nex]
        types.add(nex_class)
        current = nex
        full[nex_class] = taxa_names[nex]
        fullid[nex_class] = nex
        j += 1
    for c in allclassifications:
        if c not in full:
            full[c] = '000'
            fullid[c] = 0
    return (full, fullid)


def findFamilyGenus(taxonid, nodespath, syst=""):
    '''
    Return the family and genus for a taxonomy ID.
    
    Based on the nodes.dmp file downloaded from NCBI
    and stored in nodespath.
    
    From the NCBI website:
    nodes.dmp
    This file represents taxonomy nodes. The description for each node includes 
    the following fields:
    
    	tax_id					-- node id in GenBank taxonomy database
     	parent tax_id				-- parent node id in GenBank taxonomy database
     	rank					-- rank of this node (superkingdom, kingdom, ...) 
     	embl code				-- locus-name prefix; not unique
     	division id				-- see division.dmp file
     	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
     	genetic code id				-- see gencode.dmp file
     	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
     	mitochondrial genetic code id		-- see gencode.dmp file
     	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
     	GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
     	hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
     	comments				-- free-text comments and citation

    Parameters
    ----------
    taxonid: str
        NCBI taxonomy ID
    nodespath: str
        path to the NCBI taxonomy nodes.dmp file
    syst: str
        system option to run statements on different systems
    
    Returns
    -------
    tuple
        family taxonomy ID, genus taxonomy ID, family name, genus name
    '''
    # read the nodes.dmp file
    tab = pd.read_csv(nodespath, sep="|", header=None)
    # remove tabs from columns
    tab[2] = tab[2].str.replace("\t", "")
    # take the row of the table corresponding to this taxon id
    if taxonid == "":
        taxonid = "000"
    x = tab[tab[0] == int(taxonid)]
    # rank of current taxon ID
    types = set(x[2].values)
    #  parent of current taxon ID
    newvals = set(x[1])
    j = 0
    k = 0
    # the hierarchy above a taxon ID has various numbers of levels depending
    # on the taxon, this while loop climbs to the "family" level
    while "family" not in types and j < 15:
        # row corresponding to current taxon ID
        x2 = tab[tab[0].isin(newvals)]
        # rank of current taxon ID
        types = set(x2[2].values)
        # parent of current taxon ID
        newvals = set(x2[1])
        j += 1
    #  Sometimes family is not listed, then assign the value '000'
    if j == 15:
        fam = "000"
    else:
        # family will be the current taxon ID
        try:
            fam = (x2[0].values[0])
        except:
            fam = "000"
    # as above but for the genus level
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
    
    # Allows for one exception added to NCBI since the database was downloaded
    if taxonid == "2065263":
        genus = "6726"
        fam = "6725"
    
    # each taxon ID has a corresponding name
    # if the genus is not known the genus name can't be checked
    if genus == "000":
        genusname = "000"
    else:
        # look up the name of this taxon ID on NCBI
        genusname = findTaxonName(genus)
    if fam == "000":
        # if the family is not known the family name can't be checked
        familyname = "000"
    else:
        # look up the name of this taxon ID on NCBI
        familyname = findTaxonName(fam)
    return (fam, genus, familyname,
            genusname)


def getMetadataSRA(ID, genomesdir, nodespath, outfile,
                   log="on", syst="", check=True):
    '''
    Retrieves metadata for an SRA ID.
    
    Takes an SRA ID and retrieves the full XML record for this ID from NCBI
    taxonomy.  The following fields are extracted:
        Experiment: SRA experiment ID
        Biosample: SRA sample ID
        Type: Study type (e.g. RNA-Seq)
        TaxonID: NCBI taxonomy ID of host
        SciName: Scientific name of host
        Title: Experiment title
        Runs: SRA run IDs
        Runsize: Total size of all runs in bytes
        Spots: Number of reads
    The family and genus of the host species (and their NCBI IDs)
    are identified using the findFamilyGenus function and added to
    this table with the column names "FamilyID", "Family",
    "GenusID" and "Genus".
    
    The genomesdir directory needs to contain a file table.tsv containing
    the following columns for each reference genome which is available to map
    reads to:
        Species: scientific name of species
        TaxonID: NCBI taxonomy ID for species
        Genus: scientific name of genus
        GenusID: NCBI taxonomy ID for genus
        Family: scientific name of family
        FamilyID: NCBI taxonomy ID for family
    This table is used to check which of the available reference genomes should
    be mapped to - column Reference in the output table -
    this will ideally be a reference for the same species
    (Reference_Type H in the output table), or the same genus if the same
    species in not available (G in the output table) or the
    same family if the same genus is not available (F in the output table).
    A query to NCBI Assembly is performed to check that a closer reference
    genome is not available online.

    Parameters
    ----------
    ID: str
        SRA experiment, biosample or run ID
    genomesdir: str
        path to directory containing available reference genomes and the
        table.tsv file
    nodespath: str
        path to the NCBI taxonomy nodes.dmp file
    outfile: str
        path to output table
    log: bool
        can be "on" or "off" depending if you want to keep a log file (turn off
        for very large sets of queries)
    syst: str
        system option to run statements on different systems
    '''
    # read the reference genome table.tsv
    genometab = [line.strip().split("\t")
                 for line in open("%s/table.tsv" % genomesdir)]
    # make a list of families and genera with available reference genomes
    family = []
    genus = []
    species = []
    speciesnames = []
    for g in genometab:
        speciesnames.append(g[0])
        species.append(g[1])
        genus.append(g[3])
        family.append(g[5])
    print (species)
    # allows for SRA datasets with no species ID
    species.append("000")
    # store the family and the genus of each reference in a dictionary
    fdict = dict(zip(family, speciesnames))
    gdict = dict(zip(genus, speciesnames))

    x = 0
    while True:
        # Retrieve the full XML record for the ID
        # this query often seems to fail for no reason - keep trying
        # until it doesn't or for 10 tries
        s1 = 'efetch -db sra -format xml -id %s' % ID
     #   s1 = NCBI.fixStatement(s1)
        s2 = 'esearch -db sra -query %s | efetch -format xml' % ID
     #   s2 = NCBI.fixStatement(s2)
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
    # This is dependent on the layout of the XML file provided by NCBI
    d = xdict['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']

    try:
        experiment = d['EXPERIMENT']
    except:
        d = d[0]
        experiment = d['EXPERIMENT']

    expid = experiment['@accession']
    sample = d['SAMPLE']
    # Type
    libtype = experiment['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY']
    pe = experiment['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_LAYOUT']
    # TaxonID
    taxonid = sample['SAMPLE_NAME']['TAXON_ID']

    if 'SCIENTIFIC_NAME' in sample['SAMPLE_NAME']:
        sciname= "_".join(sample['SAMPLE_NAME']['SCIENTIFIC_NAME'].replace(" ", "_").split("_")[:2])
        sciname = sciname.replace(" ", "_")
        taxonid = findTaxonID(sciname)
    else:
        taxid = int(sample['SAMPLE_NAME']['TAXON_ID'])
        sciname = "_".join(findTaxonName(taxid).replace(" ", "_").split("_")[:2])
        sciname = sciname.replace(" ", "_")
        taxonid = findTaxonID(sciname)


    # Biosample
    biosample = sample['@accession']
    # Title
    title = experiment['TITLE']
    # Spots
    try:
        spots = d['Pool']['Member']['@spots']
    except:
        spots = 0

    # There might be multiple runs for a single experiment or sample ID
    runs = d['RUN_SET']['RUN']
    # if there is only one run it will be a dictionary with a key "@accession"
    if '@accession' in runs:
        run = d['RUN_SET']['RUN']['@accession']
    # otherwise it is a list of dictionaries - loop through all of them
    else:
        run = []
        runsize = 0
        for res in runs:
            if '@accession' in res:
                run.append(res['@accession'])
        # combine run IDs into a single string delimited by "."
        run = ".".join(run)
    # as for the run ID, run size can be a dictionary or a list of dictionaries
    if '@size' in runs:
        runsize = int(d['RUN_SET']['RUN']['@size'])
    else:
        runsize = 0
        # add together the size if there are multiple runs
        if '@size' in res:
            runsize += int(res['@size'])
    
    # make a table row with these columns
    r = ((expid, biosample, libtype, taxonid, sciname,
          title, run, runsize, spots))
    r = [r]
    # convert to pandas dataframe
    results = pd.DataFrame(r)
    results.columns = ['Experiment', 'Biosample', 'Type', 'TaxonID',
                       'SciName', 'Title', 'Runs', 'Runsize', 'Spots']
    # make the human readable run size - convert to GB and add "G"
    results['Runsize_h'] = (
        results['Runsize'] / 1e9).round(2).astype(str) + "G"
    results['SciName'] = results['SciName'].str.replace(" ", "_")

    # Look up the family and genus ID for the host species
    familyID, genusID, family, genus = findFamilyGenus(taxonid,
                                                       nodespath)


    genusID = str(genusID)
    familyID = str(familyID)

    # add these to the table row
    results['FamilyID'] = familyID
    results['GenusID'] = genusID
    results['Family'] = family
    results['Genus'] = genus

    i = 0
    while True:
        statement = '''esearch -db assembly -query "txid%s[Organism]" | \
        efetch -format uilist''' % taxonid
       #  statement = NCBI.fixStatement(statement)
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
         #   statement = NCBI.fixStatement(statement)
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
            #    statement = NCBI.fixStatement(statement)
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



    statement = '''esearch -db assembly -query "txid%s[Organism]" | \
                 efetch -format uilist''' % taxonid
  #  statement = NCBI.fixStatement(statement)
    if log == "on":
        ut_functions.writeCommand(statement, ID)
    x = Run.systemPopen(statement, syst)
    if len(x) >= 1:
        # if this is true there is a reference for this species
        reftype = "H"
    else:
        # try genus level if species level does not exist
        statement = '''esearch -db assembly -query "txid%s[Organism]" \
                      | efetch -format uilist''' % genusID
       # statement = NCBI.fixStatement(statement)
        if log == "on":
            ut_functions.writeCommand(statement, ID)
        x = Run.systemPopen(statement, syst)
        # occasionally family is listed as genus - genus should never end with
        # dae and family always should
        if len(x) >= 1 and not genus.endswith("dae"):
            # if this is true there is a genus level reference
            reftype = "G"
        else:
            # a "dae" ending is a family name
            if genus.endswith("dae"):
                family = genus
                familyID = genusID
            # look for a family level reference
            statement = '''esearch -db assembly \
                          -query "txid%s[Organism]" \
                          | efetch -format uilist''' % familyID
         #   statement = NCBI.fixStatement(statement)
            if log == "on":
                ut_functions.writeCommand(statement, ID)
            x = Run.systemPopen(statement, syst)
            if len(x) >= 1:
                # if this is true there is a family level reference
                reftype = "F"
            else:
                # there is no reference for any member of this family
                reftype = "N"
                reference = "X"

    results['Reference_Type'] = reftype

    # The NCBI databases are messy - these are specific edge cases
    # where the appropriate reference exists but can't be determined
    # automatically
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
        if sciname == "Pleurobrachia_":
            sciname = "Pleurobrachia_bachei"
            taxonid = "34499"
        if "metagenome" in sciname or "mixed" in sciname:
            sciname = "X"
            taxonid = "000"
        taxonid = str(taxonid)
        genusID = str(genusID)
        familyID = str(familyID)
        print (taxonid)
        # check that the reference genome is available locally
        if check:
            assert taxonid in species, "Host %s reference genome exists but not found" % sciname
        reference = sciname
    elif reftype == "G":
        if genus == "Penaeus":
            genusID = "133894"
        # check that the reference genome is available locally
        if check:
            if genus != "Limacina":
                assert genusID in gdict, "Genus %s reference genome exists but not found"  % genus
        if genusID in gdict:
            reference = gdict[genusID]
        else:
            reference = "000"
        
        
    elif reftype == "F":

        if familyID == "6913":
            familyID = "450948"

        # check that the reference genome is available locally
        if check:
            if family != "Limacinidae":
                assert familyID in fdict, "Family %s reference genome exists but not found" % family
        if familyID in fdict:
            reference = fdict[familyID]
        else:
            reference = "000"
    # Store the appropriate reference in the Reference column
    results['Reference'] = reference
    # Save the dataframe
    results.to_csv(outfile, sep="\t", index=None)
