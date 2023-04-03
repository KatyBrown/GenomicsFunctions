#!/usr/bin/env python3
import orfipy_core
import Bio.Seq
import pandas as pd

def getORFs(sD, minnt, stopstop=True, p3=True, p5=True, strand='b'):

    seqD = dict()
    orfD = dict()
    tabD = dict()
    for i, (acc, seq) in enumerate(sD.items()):
        seq = seq.upper()
        # Check the number of characters in the sequence -
        # sometimes they are amino acid sequences
        # and sometimes nucleotide
        resi = set(list(seq))
        # Less than 10 characters should always be nt
        # assert len(resi) < 15, "%s This might be an amino acid sequence" % seq
        # Find ORFs which are >minnt nt between two stop codons including partial at either end
        # on both strands
        # Seems to combine all tables if a table isn't specified
        orfs = orfipy_core.orfs(seq,
                                minlen=minnt,
                                between_stops=stopstop,
                                partial3=True,
                                partial5=True,
                                strand=strand)
        # If there are any ORFs meeting these criteria
        if len(orfs) != 0:
            trans = []
            tabs = []
            for o in orfs:
                # Get the corresponding sequence segment
                oseq = seq[o[0]:o[1]]
                oseq = oseq.replace("X", "N")
                # RC if - strand
                if o[2] == "-":
                    oseq = Bio.Seq.reverse_complement(oseq)
                # translate
                transeq = Bio.Seq.translate(oseq)
                table = 1
                # Check for stop codons - if there are any, try other
                # tables until you find one without
                if "*" in transeq:
                    while "*" in transeq:
                        transeq = Bio.Seq.translate(oseq, table=table)
                        table += 1

                trans.append(transeq)
                tabs.append(table)
            seqD[acc] = trans
            orfD[acc] = orfs
            tabD[acc] = tabs
    return (seqD, orfD, tabD)


def makeORFTable(seqD, orfD, tabD):
    tab = pd.DataFrame(columns=['ID', 'Sequence',
                                'ORF_Type', 'ORF_Len', 'ORF_Frame',
                                'Start_Codon', 'Stop_Codon',
                                'Start', 'End', 'Table'])
    for key in seqD:
        seqs = seqD[key]
        md = orfD[key]
        tabs = tabD[key]
        Z = zip(seqs, md, tabs)
        z = pd.DataFrame(Z)
        z['ID'] = key
        z['Sequence'] = z[0]
        z['Strand'] = [x[2] for x in z[1]]
        for section in ['type', 'len', 'frame']:
            z["ORF_%s" % section.title()] = [
                x[3].split("%s=" % section)[1].split(";")[0] for x in z[1]]
        for section in ['Start', 'Stop']:
            z["%s_Codon" % section] = [
                x[3].split("%s:" % section)[1].split(";")[0] for x in z[1]]
        z['Start'] = [x[0] for x in z[1]]
        z['End'] = [x[1] for x in z[1]]
        z['Table'] = z[2]
        z = z.drop(0, axis=1)
        z = z.drop(1, axis=1)
        z = z.drop(2, axis=1)
        tab = pd.concat([tab, z])
    return (tab)
