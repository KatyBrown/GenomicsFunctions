#!/usr/bin/env python3
import ut_functions
import numpy as np
import NCBI
import pandas as pd
import glob
import re
import Bio.Seq
import ORFs


def showAllProts(recordD):
    prots = set()
    for acc, record in recordD.items():
        for feat in record['GBSeq_feature-table']:
            quals = feat['GBFeature_quals']
            qD = dict()
            for qual in quals:
                if 'GBQualifier_value' in qual:
                    qD[qual['GBQualifier_name']] = qual['GBQualifier_value']
            if 'product' in qD:
                prots.add(qD['product'])
    print (prots)


def findSplit(recordD, n):
    rand = np.random.choice(list(recordD.keys()), n)
    for key, record in recordD.items():
        nam = record['GBSeq_definition']
        print (nam)


def getTaxa(recordD, DD, typ="aa"):
    if typ == "nt":
        recordD = {x: recordD[DD['nuc_nam_D'][x]] for x in DD['protD']}

    tax_ids = dict()
    for record_nam, record in recordD.items():
        txid = None
        for feature in record['GBSeq_feature-table']:
            quals = feature['GBFeature_quals']
            for qual in quals:
                if qual['GBQualifier_name'] == 'db_xref':
                    if qual['GBQualifier_value'].startswith('taxon:'):
                        txid = qual['GBQualifier_value'].split(":")[1]
        assert txid is not None
        tax_ids[record_nam] = txid
    records_tax = NCBI.getRecords(list(tax_ids.values()), 50, 'taxonomy')
    taxD = dict()
    for acc, record in records_tax.items():
        taxD[acc] = NCBI.cleanNCBITaxonomy(record)
    taxtab = NCBI.makeTaxTab(taxD)
    DD['tax_ids'] = tax_ids
    DD['taxtab'] = taxtab
    return (DD)


def RecordsToProts(recordD, rdrp_strings, outf, names,
                   DD, suffix, typ='aa'):
    out = open(outf, "w")
    if typ == "nt":
        outnt = open(outf.replace(".fasta", "_nt.fasta"), "w")
    for acc, record in recordD.items():
        x = 0
        pp = dict()
        for feat in record['GBSeq_feature-table']:
            if not isinstance(feat, list):
                feat = [feat]
            qD = dict()
            for f in feat:
                if 'GBFeature_quals' in f:
                    quals = f['GBFeature_quals']
                    for qual in quals:
                        if 'GBQualifier_value' in qual:
                            qD[qual['GBQualifier_name']] = qual['GBQualifier_value']
            if 'product' in qD:
                y = 0
                for string in rdrp_strings:
                    if string in qD['product']:
                        y += 1
                if y != 0:
                    prod = record['GBSeq_organism']
                    virus_nam = ut_functions.clean_string(prod)
                    if virus_nam[0] == virus_nam[0].lower():
                        virus_nam = virus_nam[0].upper() + virus_nam[1:]
                    if 'protein_id' in qD:
                        acc_prot = qD['protein_id']
                        prot_seq = qD['translation']
                        nuc_seq = record['GBSeq_sequence']
                        virus_id = "%s|%s|%s" % (acc_prot, virus_nam, suffix)
    
                        DD['accD'][acc_prot] = virus_id
                        DD['nucD'][acc_prot] = nuc_seq
                        DD['protD'][acc_prot] = prot_seq
                        DD['namD'][acc_prot] = names[acc]
                        DD['nuc_nam_D'][acc_prot] = acc
                        x += 1

                else:
                    if typ == "nt":
                        if 'translation' in qD and 'protein_id' in qD:
                            pp["%s|%s|%s" % (acc, qD['protein_id'], qD['product'].replace(" ", "_"))] = qD['translation']

        if x == 0:
            if typ == "nt":
                for prod in pp:
                    current = pp[prod]
                    out.write(">%s\n%s\n" % (prod, current))
                if len(pp) == 0:
                    outnt.write(">%s\n%s\n" % (acc, record['GBSeq_sequence']))
            else:
                out.write(">%s\n%s\n" % (acc, record['GBSeq_sequence']))
    out.close()
    return (DD)


def makeMDTable(DD):
    df = pd.DataFrame(DD['accD'].items(), columns=['Accession', 'Virus_ID'])
    df['Virus_Name'] = df['Virus_ID'].str.split("|").str.get(1)
    df['Virus_Source'] = df['Virus_ID'].str.split("|").str.get(2)
    df = df.merge(pd.DataFrame(DD['namD'].items(),
                               columns=['Accession', 'Original_Name']),
                  'right')
    df['Nuccore_Accession'] = [DD['nuc_nam_D'][x] for x in df['Accession']]
    df = df.merge(pd.DataFrame(DD['tax_ids'].items(),
                               columns=['Accession', 'TaxID']), 'outer')
    df['TaxID'] = df['TaxID'].fillna(0)
    df['TaxID'] = df['TaxID'].astype(int)
    df = df.merge(DD['taxtab'], 'left')
    df['TaxID'] = df['TaxID'].astype(int)
    df['Accession_Short'] = df['Accession'].str.split("ntpp").str.get(0)
    df['Original_Taxon'] = [DD['namD'][x] for x in df['Accession']]
    print (len(df))
    return (df)


def writeAll(DD, df, auth):
    donef = glob.glob("../done*txt")
    M = set()
    for d in donef:
        M.add(int(d.split("/")[-1].split(".")[0].replace("done", "")))
    done = set([line.strip() for line in open("../done%i.txt" % max(M)).readlines()])
    df = df[~df['Nuccore_Accession'].isin(done)]
    df = df[~df['Accession'].isin(done)]
    df = df[~df['Accession_Short'].isin(done)]
    df.to_csv("metadata_table.tsv", sep="\t", index=None)
    out = open("clean_gb_fasta_nucleotide_%s.fasta" % auth, "w")
    for nam, seq in DD['nucD'].items():
        if nam not in done:
            out.write(">%s\n%s\n" % (DD['accD'][nam],
                                     seq.upper().replace(
                                         "-", "").replace("X", "")))
    out.close()
    out = open("clean_gb_fasta_%s.fasta" % auth, "w")
    for nam, seq in DD['protD'].items():
        if nam not in done:
            out.write(">%s\n%s\n" % (DD['accD'][nam],
                                     seq.upper().replace(
                                         "-", "").replace("X", "")))
    out.close()
 
    done = done | set(df['Accession_Short']) | set(
        df['Accession']) | set(
            df['Nuccore_Accession'])
    outd = open("../done%i.txt" % (max(M) + 1), "w")
    for d in done:
        outd.write("%s\n" % d)
    outd.close()


def getAccRange(startacc, endacc):
    mm = re.match("[A-Z]+", startacc).span()
    txt = startacc[mm[0]:mm[1]]
    startstr = startacc[mm[1]:]
    startind = int(startacc[mm[1]:])
    endind = int(endacc[mm[1]:])
    accs = []
    for ind in np.arange(startind, endind+1):
        acc = "%s%s" % (txt, str(ind).zfill(len(startstr)))
        accs.append(acc)
    return (accs)


def namesFromOrganism(record, acc):
    org = record['GBSeq_organism']
    org = ut_functions.clean_string(org)
    if org[0] == org[0].lower():
        org = org[0].upper() + org[1:]
    return (org)

    
def BlastToSeqsNt(blast, F, records, names, DD, suffix):
    for ind in blast.index.values:
        row = blast.loc[ind]
        start = row['start_pos_query']
        end = row['end_pos_query']
        if start < end:
            subseq = F[row['query']][start:end]
            trans = Bio.Seq.translate(subseq[2:])
            strand = "p"
            ss = start
            ee = end
        else:
            subseq = F[row['query']][end:start]
            rc = Bio.Seq.reverse_complement(subseq)
            trans = Bio.Seq.translate(rc)
            strand = "m"
            ss = end
            ee = start
        virus_nam = names[row['query']]
        acc = "%s_ntpp%i_%i_%s" % (row['query'], ss, ee, strand)
        virus_ID = "%s|%s|%s" % (acc, virus_nam, suffix)
        DD['accD'][acc] = virus_ID
        DD['namD'][acc] = row['query']
        DD['nuc_nam_D'][acc] = row['query']
        DD['protD'][acc] = trans
        DD['nucD'][acc] = F[row['query']]
    return (DD)


def BlastToSeqsORFs(nt_fasta, blast, oD, sD, names, suffix, DD):
    for ind in blast.index.values:
        row = blast.loc[ind]
        q = row['query']
        qnam, orfn = q.split("|")
        qnam2 = qnam.replace("_", " ")
        orfn = int(orfn)
        orf_details = oD[qnam2][orfn]
        start, end, strandid = orf_details[0], orf_details[1], orf_details[2]
        if strandid == "+":
            strand = "p"
        else:
            strand = "m"
        ss = nt_fasta[qnam2][start:end]
        if strand == "m":
            ss = Bio.Seq.reverse_complement(ss)
        transeq = Bio.Seq.translate(ss)
        assert transeq == sD[qnam2][orfn]
        pos = "ntpp%s_%s_%s" % (start, end, strand)

        virus_name = names[q]
        acc = "noacc_%s" % (virus_name)
        virus_id = "%s|%s|%s" % (acc, virus_name, suffix)
        DD['accD'][acc] = virus_id
        DD['protD'][acc] = transeq
        DD['nucD'][acc] = ss
        DD['nuc_nam_D'][acc] = qnam
        DD['namD'][acc] = virus_name
    return (DD)

def BlastToSeqsAA(blast, F, records, suffix, DD):
    for ind in blast.index.values:
        row = blast.loc[ind]
        query = row['query']
        prot_acc = query.split("|")[1]
        nuc_acc = query.split("|")[0]
        start = row['start_pos_query']
        end = row['end_pos_query']
        full_seq = F[query]
        subseq = full_seq[start:end]
        assert start < end
        strand = "p"
        record = records[nuc_acc] 
        virus_nam = namesFromOrganism(record, nuc_acc)
        acc = "%s_aapp%i_%i_%s" % (prot_acc, start, end, strand)
        virus_ID = "%s|%s|%s" % (acc, virus_nam, suffix)
        DD['protD'][acc] = subseq
        DD['nucD'][acc] = record['GBSeq_sequence']
        DD['namD'][acc] = record['GBSeq_organism']
        DD['accD'][acc] = virus_ID
        DD['nuc_nam_D'][acc] = nuc_acc
    return (DD)

def ntDictToORFs(F, outf, minsize):
    sD, oD, tD = ORFs.getORFs(F, minsize)
    out = open(outf, "w")
    for key, orfs in sD.items():
        for o, seq in enumerate(orfs):
            out.write(">%s|%i\n%s\n" % (key.replace(" ", "_"), o, seq))
    out.close()
    return (oD, sD)

def splitProtBLAST(blast):
    tab = pd.read_csv(blast, sep="\t", header=None)
    tab.columns = ut_functions.getBLASTcolumns()
    tab['query_nt'] = tab['query'].str.split("|").str.get(0)
    tab = tab.sort_values('bit_score', ascending=False)
    best = tab.groupby('query_nt').first()
    best['query_nt'] = best.index.values
    best.index = np.arange(len(best))
    return(best)