#!/usr/bin/env python3
import orfipy_core
import Bio.Seq


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
        assert len(resi) < 10, "%s This might be an amino acid sequence" % seq
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

