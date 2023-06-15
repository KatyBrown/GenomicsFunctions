#!/usr/bin/env python3
import ut_functions
import numpy as np
import NCBI
import pandas as pd
import glob
import re
import Bio.Seq
import ORFs
import os


def filterDD(accs, DD):
    DD2 = dict()
    for key in DD:
        if key != 'taxtab':
            DD2[key] = dict()
    
        for acc in accs:
            DD2[key][acc] = DD[key][acc]

    if 'taxtab' in DD:
        DD2['taxtab'] = DD['taxtab']
    return (DD2)


def getORFPosHMMER(hmmer_tab, orfF, orfD, ntF):
    posD = dict()
    for query in hmmer_tab['query_name']:
        nt, ind = query.split("|")
        ind = int(ind)
        orf_seq = orfF[nt][ind]
        nt_seq = ntF[nt]
        start, end, strand, details = orfD[nt][ind]
        nt_subseq = nt_seq[start:end]
        if strand == "-":
            nt_subseq = Bio.Seq.reverse_complement(nt_subseq)
        trans = Bio.Seq.translate(nt_subseq)
        assert trans == orf_seq
        pos = "%s:%s-%s%s_pp" % (nt, start, end, strand)
        posD[nt] = pos, orf_seq
    return (posD)

def posDToProts(posD, recordD, suffix, DD):
    for query, (pos, oseq) in posD.items():
        record = recordD[query]
        nt_seq = record['GBSeq_sequence'].upper()
        virus_nam = ut_functions.clean_string(record['GBSeq_organism'])
        if virus_nam[0] == virus_nam[0].lower():
            virus_nam = virus_nam[0].upper() + virus_nam[1:]
        virus_id = "%s|%s|%s" % (pos, virus_nam, suffix)
        desc = record['GBSeq_definition']
        DD['accD'][pos] = virus_id
        DD['protD'][pos] = oseq
        DD['namD'][pos] = virus_nam
        DD['nuc_nam_D'][pos] = query
        DD['descD'][pos] = desc
        DD['nucD'][pos] = nt_seq
    return (DD)        


def showAllProts(recordD):
    prots = set()
    for acc, record in recordD.items():
        for feat in record['GBSeq_feature-table']:
            if 'GBFeature_quals' in feat:
                quals = feat['GBFeature_quals']
                qD = dict()
                for qual in quals:
                    if 'GBQualifier_value' in qual:
                        qD[qual['GBQualifier_name']] = qual['GBQualifier_value']
                if 'product' in qD:
                    prots.add(qD['product'])
    return (prots)


def findSplit(recordD, n):
    rand = np.random.choice(list(recordD.keys()), n)
    for key, record in recordD.items():
        nam = record['GBSeq_definition']
        print (nam)


def getTaxa(recordD, DD, typ="aa", silent=False):
    if typ == "nt":
        recordD = {x: recordD[DD['nuc_nam_D'][x]] for x in DD['protD']}

    tax_ids = dict()
    for record_nam, record in recordD.items():
        txid = None
        for feature in record['GBSeq_feature-table']:
            if 'GBFeature_quals' in feature:
                quals = feature['GBFeature_quals']
                for qual in quals:
                    if qual['GBQualifier_name'] == 'db_xref':
                        if qual['GBQualifier_value'].startswith('taxon:'):
                            txid = qual['GBQualifier_value'].split(":")[1]
        assert txid is not None
        tax_ids[record_nam] = txid
    records_tax = NCBI.getRecords(list(tax_ids.values()), 50, 'taxonomy', silent=silent)
    taxD = dict()
    for acc, record in records_tax.items():
        taxD[acc] = NCBI.cleanNCBITaxonomy(record)
    taxtab = NCBI.makeTaxTab(taxD)
    DD['tax_ids'] = tax_ids
    DD['taxtab'] = taxtab
    return (DD)

def makeBlankDict():
    DD = dict()
    DD['accD'] = dict()
    DD['nucD'] = dict()
    DD['protD'] = dict()
    DD['namD'] = dict()
    DD['nuc_nam_D'] = dict()
    DD['descD'] = dict()
    return (DD)


def addNtSeqs(recordD, DD):
    for protacc, ntacc in DD['nuc_nam_D'].items():
        record = recordD[ntacc]
        record_seq = record['GBSeq_sequence']
        DD['nucD'][protacc] = record_seq
    return (DD)


def ProtRecordsToProts(recordD, DD, suffix):  
    for acc, record in recordD.items():
        nt_acc = record['GBSeq_source-db'].split(" ")[-1]
        prot_seq = record['GBSeq_sequence'].upper()
        virus_nam = ut_functions.clean_string(record['GBSeq_organism'])
        if virus_nam[0] == virus_nam[0].lower():
            virus_nam = virus_nam[0].upper() + virus_nam[1:]
        virus_id = "%s|%s|%s" % (acc, virus_nam, suffix)
        desc = record['GBSeq_definition']
        DD['accD'][acc] = virus_id
        DD['protD'][acc] = prot_seq
        DD['namD'][acc] = virus_nam
        DD['nuc_nam_D'][acc] = nt_acc
        DD['descD'][acc] = desc
    return (DD)


def RecordsToProts(recordD, rdrp_strings, rdrp_not, outf, names,
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
                    if string in qD['product'].lower():
                        w = 0
                        b = 0
                        for string2 in rdrp_not:
                            if string2 in qD['product'].lower():
                                w += 1
                        if w == 0 and b == 0:
                            y += 1
                if y != 0 or len(rdrp_strings) == 0:
                    prod = record['GBSeq_organism']
                    virus_nam = ut_functions.clean_string(prod)
                    if virus_nam[0] == virus_nam[0].lower():
                        virus_nam = virus_nam[0].upper() + virus_nam[1:]
                    if 'protein_id' in qD and 'translation' in qD:
                        acc_prot = qD['protein_id']
                        prot_seq = qD['translation']
                        nuc_seq = record['GBSeq_sequence']
                        virus_id = "%s|%s|%s" % (acc_prot, virus_nam, suffix)
    
                        DD['accD'][acc_prot] = virus_id
                        DD['nucD'][acc_prot] = nuc_seq
                        DD['protD'][acc_prot] = prot_seq
                        DD['namD'][acc_prot] = names[acc]
                        DD['nuc_nam_D'][acc_prot] = acc
                        DD['descD'][acc_prot] = qD['product']

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
    #df['Original_Taxon'] = [DD['namD'][x] for x in df['Accession']]
    #print (len(df))
    return (df)


def writeAll(DD, df, auth):

    df.to_csv("data/metadata_table.tsv", sep="\t", index=None)
    out = open("data/clean_gb_fasta_nucleotide_%s.fasta" % auth, "w")
    for nam, seq in DD['nucD'].items():
            out.write(">%s\n%s\n" % (DD['accD'][nam],
                                     seq.upper().replace(
                                         "-", "").replace("X", "")))
    out.close()
    out = open("data/clean_gb_fasta_%s.fasta" % auth, "w")
    for nam, seq in DD['protD'].items():
            out.write(">%s\n%s\n" % (DD['accD'][nam],
                                     seq.upper().replace(
                                         "-", "").replace("X", "")))
    out.close()
    try:
        os.mkdir("data/subsets")
    except:
        pass
    for clasi in set(df['Classification']):
        subtab = df[df['Classification'] == clasi]
        out = open("data/subsets/%s.fasta" % clasi, "w")
        for nam in set(subtab['Accession']):
            out.write(">%s\n%s\n" % (DD['accD'][nam],
                                     seq.upper().replace(
                                         "-", "").replace("X", "")))
        out.close()


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