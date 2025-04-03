#!/usr/bin/env python3
import orfipy_core
import Bio.Seq
import pandas as pd
import re


def getORFs(sD, minnt, stopstop=True, p3=True, p5=True, strand='b',
            starts=['TTG','CTG','ATG'],
            stops=['TAA','TAG','TGA']):

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
        assert len(resi) < 15, "%s This might be an amino acid sequence" % seq
        # Find ORFs which are >minnt nt between two stop codons including partial at either end
        # on both strands
        # Seems to combine all tables if a table isn't specified
        orfs = orfipy_core.orfs(seq,
                                minlen=minnt,
                                between_stops=stopstop,
                                partial3=True,
                                partial5=True,
                                strand=strand,
                                starts=starts,
                                stops=stops)
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
    x = 0
    for key in seqD:
        if x % 1000 == 0:
            print (x)
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
        x += 1
    return (tab)


def getLongestORFs(orfD, seqD):
    longestD = dict()
    for orf, resL in orfD.items():
        longest = ""
        longest_len = 0
        for i, res in enumerate(resL):
            start, end, strand, details = res
            length = end - start
            seq = seqD[orf][start:end]
            if strand == '-':
                seq = Bio.Seq.reverse_complement(seq)
            assert len(seq) == length
            tseq = Bio.Seq.translate(seq)
            sections = re.split("XXX", tseq)
            for j, sec in enumerate(sections):
                if len(sec)> longest_len:
                    longest = i, j
                    longest_len = len(sec)
        lstart, lend, lstrand, ldetails = resL[longest[0]]
        lseq = seqD[orf][lstart:lend]
        if lstrand == '-':
            lseq = Bio.Seq.reverse_complement(lseq)
        ltseq = Bio.Seq.translate(lseq)
        lsections = re.split("XXX", ltseq)
        lpos = re.search(lsections[longest[1]], ltseq)
        if lpos:
            ts, te = lpos.span()
            if lstrand == '+':
                true_start = (ts * 3) + lstart
                true_end = (te * 3) + lstart
                true_seq = seqD[orf][true_start:true_end]
            else:
                true_start = lend - (te * 3)
                true_end = lend - (ts * 3)
                true_seq = seqD[orf][true_start:true_end]
                true_seq = Bio.Seq.reverse_complement(true_seq)
            t_true_seq = Bio.Seq.translate(true_seq)
            assert t_true_seq == lsections[longest[1]]
            longestD[orf] = longest, longest_len, true_start, true_end, lsections[longest[1]]
        else:
            print("ORF not found in %s" % orf)
    return (longestD)


def getLongORFs(orfD, seqD, minlen):
    longestD = dict()
    for orf, resL in orfD.items():
        longestD.setdefault(orf, [])
        longest = []
        longest_len = []
        for i, res in enumerate(resL):
            start, end, strand, details = res
            length = end - start
            seq = seqD[orf][start:end]
            if strand == '-':
                seq = Bio.Seq.reverse_complement(seq)
            assert len(seq) == length
            tseq = Bio.Seq.translate(seq)
            sections = re.split("XXX", tseq)
            for j, sec in enumerate(sections):
                if len(sec)> minlen:
                    longest.append([i, j])
                    longest_len.append(len(sec))
        for lest, lest_len in zip(longest, longest_len):
            lstart, lend, lstrand, ldetails = resL[lest[0]]
            lseq = seqD[orf][lstart:lend]
            if lstrand == '-':
                lseq = Bio.Seq.reverse_complement(lseq)
            ltseq = Bio.Seq.translate(lseq)
            lsections = re.split("XXX", ltseq)
            lpos = re.search(lsections[lest[1]], ltseq)
            if lpos:
                ts, te = lpos.span()
                if lstrand == '+':
                    true_start = (ts * 3) + lstart
                    true_end = (te * 3) + lstart
                    true_seq = seqD[orf][true_start:true_end]
                else:
                    true_start = lend - (te * 3)
                    true_end = lend - (ts * 3)
                    true_seq = seqD[orf][true_start:true_end]
                    true_seq = Bio.Seq.reverse_complement(true_seq)
                t_true_seq = Bio.Seq.translate(true_seq)
                assert t_true_seq == lsections[lest[1]]
                longestD[orf].append([lest, lest_len, true_start, true_end,
                                     lsections[lest[1]]])
            else:
                print("ORF not found in %s" % orf)
    return (longestD)