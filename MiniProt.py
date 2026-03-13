import Bio.Seq
import re
import numpy as np
import pandas as pd
import ut_functions
import os


def integerSeg(seg, paramD):
    # Type 1
    typi = 1
    # If the segment is an integer, there are this number of codons
    # in which the amino acids match
    nmatches = int(seg)
    # Iterate throug them
    for n in np.arange(nmatches):
        # Take the next three remaining nucleotides in the nt seq
        codon_nt = paramD['rem_nt'][0:3]
        paramD['codons_nt'].append(codon_nt)
        paramD['rem_nt'] = paramD['rem_nt'][3:]

        # Take the next three remaining nucleotides in the ntaa seq
        codon_prot = paramD['rem_prot'][0:3]
        paramD['codons_prot'].append(codon_prot)
        paramD['rem_prot'] = paramD['rem_prot'][3:]

        # Translate both
        trans_prot = Bio.Seq.translate(codon_prot)
        trans_nt = Bio.Seq.translate(codon_nt)
        paramD['nt_aa'].append(trans_nt)
        paramD['prot_aa'].append(trans_prot)

        # Check they do translate to the same thing
        assert trans_prot == trans_nt, paramD
        # Move to the next position
        paramD['posi_nt'] += 3
    return (paramD)


def codonAASeg(seg, paramD):
    # Type 2
    typi = 2
    # Three lower case letters and one upper case letter - the
    # lower case letters are the codon in the nt sequence
    # the upper case letter is the aa in the protein sequence

    # Get the relevant codons
    codon_nt = paramD['rem_nt'][0:3]
    paramD['codons_nt'].append(codon_nt)
    paramD['rem_nt'] = paramD['rem_nt'][3:]

    codon_prot = paramD['rem_prot'][0:3]
    paramD['codons_prot'].append(codon_prot)
    paramD['rem_prot'] = paramD['rem_prot'][3:]

    trans_nt = Bio.Seq.translate(codon_nt)
    trans_prot = Bio.Seq.translate(codon_prot)

    actual = seg[-1]
    # Check that the protein sequence matches the upper case letter
    assert trans_prot == actual, (trans_prot, actual)
    # Check that the nt sequence matches the lower case letter
    assert codon_nt == seg[0:3].upper()
    paramD['nt_aa'].append(trans_nt)
    paramD['prot_aa'].append(trans_prot)
    # Move to the next position
    paramD['posi_nt'] += 3
    return (paramD)


def gapsProtSeg(seg, paramD):
    # Type 3
    typi = 3
    for indi in np.arange(0, len(seg), 3):
        # Take this fragment
        frag = seg[indi:indi+3]
        # Take from the nt sequence up to the end of the fragment length
        # then pad to three nt with gaps
        codon_nt = paramD['rem_nt'][0:len(frag)] + "-" * (3 - len(frag))
        paramD['codons_nt'].append(codon_nt)
        paramD['rem_nt'] = paramD['rem_nt'][len(frag):]

        # The matching aa sequence here is just three gaps
        codon_prot = "---"
        paramD['codons_prot'].append(codon_prot)

        # Check the fragment sequence are in the identified codon
        assert frag.upper() in codon_nt, (typi, seg, paramD)

        # If there are gaps in the nt codon, don't translate this codon or 
        # add it to the final sequence
        if "-" in codon_nt:
            trans_nt = "-"
            for i in np.arange(len(frag)):
                paramD['skipped'].append(paramD['posi_nt'] + i)
        else:
            trans_nt = Bio.Seq.translate(codon_nt)

        # If the nt sequence was translated, append to the nt sequence
        trans_prot = "-"
        if not (trans_prot == "-" and trans_nt == "-"):
            paramD['nt_aa'].append(trans_nt)
            paramD['prot_aa'].append(trans_prot)

        paramD['posi_nt'] += len(frag)
    return (paramD)


def gapsNtSeg(seg, paramD):
    # Type 4
    # If the segment is just upper case letters, these are in the protein
    # sequence but skipped in the nt sequence.
    typi = 4
    segl = len(seg)
    for i in np.arange(segl):
        # Add a blank codon
        codon_nt = "---"
        paramD['codons_nt'].append(codon_nt)

        # Translate the aant codon and check it's correct
        codon_prot = paramD['rem_prot'][0:3]
        paramD['codons_prot'].append(codon_prot)
        paramD['rem_prot'] = paramD['rem_prot'][3:]
        trans_prot = Bio.Seq.translate(codon_prot)
        assert trans_prot == seg[i], (trans_prot, seg[i], typi, seg, paramD)
        paramD['nt_aa'].append("-")
        paramD['prot_aa'].append(trans_prot)
    return (paramD)


def stopCodonProtSeq(seg, paramD):
    # Type 5
    typi = 5
    actual = "*"
    # The position is a stop codon in the aant seq
    # Translate the next three positions in the nt sequence
    codon_nt = paramD['rem_nt'][0:3]
    paramD['codons_nt'].append(codon_nt)
    paramD['rem_nt'] = paramD['rem_nt'][3:]

    # Translate the next three postions in the aant sequence
    codon_prot = paramD['rem_prot'][0:3]
    paramD['codons_prot'].append(codon_prot)
    paramD['rem_prot'] = paramD['rem_prot'][3:]

    trans_nt = Bio.Seq.translate(codon_nt)
    trans_prot = Bio.Seq.translate(codon_prot)

    # Check that the translated aant sequence is a stop codon
    assert trans_prot == actual, (typi, seg, paramD)

    paramD['nt_aa'].append(trans_nt)
    paramD['prot_aa'].append(trans_prot)
    paramD['posi_nt'] += 3
    return (paramD)


def gapsNtSeg1nt(seg, paramD):
    # Type 6 and 8
    typi = 6
    # Single nucleotide insertion in the nucleotide sequence
    # but unlike in gapNtSeg the other two positions are aligned
    # with the aant sequence
    # Pad the nt sequence with one dash
    codon_nt = paramD['rem_nt'][0:2] + "-"
    paramD['codons_nt'].append(codon_nt)
    paramD['rem_nt'] = paramD['rem_nt'][2:]

    codon_prot = paramD['rem_prot'][0:3]
    paramD['codons_prot'].append(codon_prot)
    paramD['rem_prot'] = paramD['rem_prot'][3:]

    # Don't add to the nt amino acid sequence
    trans_nt = "-"
    # The aant segment can be translated
    trans_prot = Bio.Seq.translate(codon_prot)
    assert trans_prot == seg[-1], (trans_prot, typi, seg, paramD)
    paramD['nt_aa'].append(trans_nt)
    paramD['prot_aa'].append(trans_prot)
    # Track that these two nt positions are not translated
    paramD['skipped'].append(paramD['posi_nt'])
    paramD['skipped'].append(paramD['posi_nt'] + 1)
    paramD['posi_nt'] += 2
    return (paramD)


def gapsNtSeg2nt(seg, paramD):
    # Type 7
    typi = 7
    # 2 nucleotide insertion in the nucleotide sequence
    # but the other position is aligned
    # with the aant sequence
    # Add the one nt and two gaps
    codon_nt = "%s--" % (paramD['rem_nt'][0])
    paramD['codons_nt'].append(codon_nt)
    paramD['rem_nt'] = paramD['rem_nt'][1:]

    codon_prot = paramD['rem_prot'][0:3]
    paramD['codons_prot'].append(codon_prot)
    paramD['rem_prot'] = paramD['rem_prot'][3:]

    trans_prot = Bio.Seq.translate(codon_prot)
    assert seg[0].upper() == codon_nt[0], (codon_nt, typi, seg, paramD)
    assert trans_prot == seg[-1], (typi, seg, paramD)
    #  Don't add to the nt amino acid sequence
    paramD['nt_aa'].append("-")
    # Track that this nt position is not translated
    paramD['skipped'].append(paramD['posi_nt'])
    paramD['prot_aa'].append(trans_prot)
    paramD['posi_nt'] += 1
    return (paramD)


def splice(seg, paramD, typis, pdelim, ndelim):
    segs = seg.split("~")
    leng = int(segs[1][2:-2])
    s = segs[0]
    typi1, paramD, typis = getTypeBasic(s, pdelim, ndelim, paramD, typis,
                                        splice=True)
    for m in np.arange(0, leng):
        codon_nt = paramD['rem_nt'][0]
        paramD['codons_nt'].append(codon_nt)
        paramD['skipped'].append(paramD['posi_nt'])
        paramD['posi_nt'] += 1
        paramD['rem_nt'] = paramD['rem_nt'][1:]
    return (paramD)


def gapsNtNoSeg(seg, paramD):
    # A single nucleotide in the nt seq not aligned to anything in the aa seq
    typi = 9
    codon_nt = "%s--" % paramD['rem_nt'][0]
    paramD['codons_nt'].append(codon_nt)
    paramD['rem_nt'] = paramD['rem_nt'][1:]
    codon_prot = "---"
    paramD['codons_prot'].append(codon_prot)
    paramD['skipped'].append(paramD['posi_nt'])
    paramD['posi_nt'] += 1
    return (paramD)


def gaps2NtStop(seg, paramD):
    # Two nt in the nt seq aligned to a stop in the aa seq
    typi = 13
    codon_nt = "%s%s-" % (paramD['rem_nt'][0], paramD['rem_nt'][1])
    paramD['codons_nt'].append(codon_nt)
    paramD['rem_nt'] = paramD['rem_nt'][2:]

    codon_prot = paramD['rem_prot'][0:3]
    paramD['codons_prot'].append(codon_prot)
    paramD['rem_prot'] = paramD['rem_prot'][3:]
    trans_prot = Bio.Seq.translate(codon_prot)

    assert trans_prot == "*", seg
    paramD['prot_aa'].append(trans_prot)

    paramD['skipped'].append(paramD['posi_nt'])
    paramD['posi_nt'] += 1

    paramD['skipped'].append(paramD['posi_nt'])
    paramD['posi_nt'] += 1
    return (paramD)


def stopProtNoSeg(seg, paramD):
    # Stop codon in the aant sequence aligns to nothing in the nt sequence
    # type 11
    typi = 11
    codon_nt = "---"
    paramD['codons_nt'].append(codon_nt)
    codon_prot = paramD['rem_prot'][0:3]
    paramD['codons_prot'].append(codon_prot)
    paramD['rem_prot'] = paramD['rem_prot'][3:]

    trans_prot = Bio.Seq.translate(codon_prot)

    assert trans_prot == "*", (typi, seg, paramD)
    paramD['nt_aa'].append("-")
    paramD['prot_aa'].append(trans_prot)
    return (paramD)


def insertionSCProt(seg, paramD):
    # string of upper case letters - insertion in prot seq with nothing in
    # nt seq
    typi = 10
    codon_test1 = Bio.Seq.translate(paramD['rem_prot'][0:3])
    codon_test2 = Bio.Seq.translate(paramD['rem_prot'][3:6])
    if codon_test2 == "*":
        sseg = "**" + seg
    elif codon_test1 == "*":
        sseg = "*" + seg
    else:
        sseg = seg
    for pos in sseg:
        codon_nt = "---"
        paramD['codons_nt'].append(codon_nt)
        codon_prot = paramD['rem_prot'][0:3]
        paramD['codons_prot'].append(codon_prot)
        paramD['rem_prot'] = paramD['rem_prot'][3:]
        trans_prot = Bio.Seq.translate(codon_prot)
        assert trans_prot == pos, (codon_prot, pos, trans_prot, typi, seg,
                                   paramD)
        paramD['nt_aa'].append("-")
        paramD['prot_aa'].append(trans_prot)
    return paramD


def cleanIUPAC(seq):
    # Replace IUPAC codes in the nt sequence 
    for nuc in ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']:
        seq = seq.replace(nuc, "N")
    return (seq)


def getTypeBasic(seg, pdelim, ndelim, paramD, typis, splice=False):
    isint = None
    typi = None
    try:
        int(seg)
        isint = True
    except ValueError:
        isint = False

    if isint:
        typi = 1
        paramD = integerSeg(seg,
                            paramD)
        typis.append(typi)
    elif len(seg) == 4 and seg[0:3].lower() == seg[0:3] and seg[3].upper() == seg[3]:
        typi = 2
        paramD = codonAASeg(seg, paramD)
        typis.append(typi)
    elif pdelim == "-" or (pdelim == "!" and seg.lower() == seg):
        typi = 3
        paramD = gapsProtSeg(seg, paramD)
        typis.append(typi)
    elif pdelim == "+":
        typi = 4
        paramD = gapsNtSeg(seg, paramD)
        typis.append(typi)
    elif len(seg) == 3 and ndelim == "*":
        if seg == seg.lower():
            typi = 5
            paramD = stopCodonProtSeq(seg, paramD)
            typis.append(typi)
        elif seg != seg.upper():
            typi = 6
            paramD = gapsNtSeg1nt(seg, paramD)
            typis.append(typi)
        else:
            typi = 14
            paramD = insertionSCProt(seg, paramD)
            typis.append(typi)
    else:
        segl = len(seg[:-1])
        if segl == 1 and seg.upper() != seg and seg.lower() != seg:
            typi = 7
            paramD = gapsNtSeg2nt(seg, paramD)
            typis.append(typi)
        elif segl == 1 and seg.upper != seg and seg.lower() == seg:
            typi = 13
            paramD = gaps2NtStop(seg, paramD)
            typis.append(typi)
        elif segl == 2 and seg.upper() != seg:
            typi = 8
            paramD = gapsNtSeg1nt(seg, paramD)
            typis.append(typi)
        elif segl == 0 and seg.upper() != seg:
            typi = 9
            paramD = gapsNtNoSeg(seg, paramD)
            typis.append(typi)
        elif seg.upper() == seg:
            typi = 10
            paramD = insertionSCProt(seg, paramD)

    return (typi, paramD, typis)


def align_from_cs(cs_string,
                  full_nt_for_prot_seq,
                  full_nt_for_nt_seq,
                  full_prot_seq,
                  prot_start,
                  prot_end, 
                  cont_start,
                  cont_end,
                  strand, ID):
    '''
    In the comments:
    nt seq = the sequence which was nucleotide in the original alignment
    aant seq = the nt version of the sequence which was amino acids in the
    original alignment
    protein seq = the aa version of the sequence which was amino acids in the
    original alignment
    '''
    p = 0
    q = 0
    # Divide the cigar string into codons (or integers where the integer is
    # the number of codons which aligned perfectly)
    segs = re.split("\*|:|\+|\-", cs_string)[3:]
    # Split again but include delimiters
    segs_delim = re.split(r'(\*|:|\+|\-)', cs_string)[6:]

    full_nt_for_nt_seq = cleanIUPAC(full_nt_for_nt_seq)
    full_nt_for_prot_seq = cleanIUPAC(full_nt_for_prot_seq)


    # Crop to only the part of the nt sequence aligned by MiniProt
    cropped_nt_for_nt_seq = full_nt_for_nt_seq[cont_start:cont_end]

    # Reverse translate if needed
    if strand == "-":
        cropped_nt_for_nt_seq = Bio.Seq.reverse_complement(
            cropped_nt_for_nt_seq)

    # Crop the protein sequence based on the miniprot output
    cropped_prot_seq = full_prot_seq[prot_start:prot_end].upper()

    # Crop the aant sequence to match
    cropped_nt_for_prot_seq = full_nt_for_prot_seq[
        (prot_start * 3):(prot_end * 3)]

    paramD = {'codons_nt': [],
              'codons_prot': [],
              'rem_nt': cropped_nt_for_nt_seq,
              'rem_prot': cropped_nt_for_prot_seq,
              'nt_aa': [],
              'prot_aa': [],
              'posi_nt': 0,
              'skipped': [],
              'cs_string': cs_string,
              'paf': ID}
    typis = [""]

    # Iterate through the codons
    for seg in segs:
        if len(seg) != 0:
            #seg = seg.replace("~", "")
            # Check if the codon is an integer
            # Check there are remaining nucleotides
            assert len(paramD['rem_nt']) != 0
            if "~" not in seg:
                if p < len(segs_delim)-1:
                    typi, paramD, typis = getTypeBasic(seg,
                                                    segs_delim[p-1],
                                                    segs_delim[p+1],
                                                    paramD, typis)
                else:
                    typi, paramD, typis = getTypeBasic(seg,
                                                    segs_delim[p-1],
                                                    "-",
                                                    paramD, typis)
            elif "~" in seg:
                typi = 15
                paramD = splice(seg, paramD, typis,
                                segs_delim[p-1], segs_delim[p+1])
                typis.append(typi)
                segs_delim[p+1] = "!"
            else:
                raise RuntimeError("couldn't process segment", (seg,
                                                                paramD))

            px = 0
        else:
            if segs_delim[p-1] == "*":
                if segs[q-1].lower() == segs[q-1] and len(segs[q-1]) == 3:
                    typi = 12
                    typis.append(typi)
                elif typis[-1] != 13:
                        typi = 11
                        paramD = stopProtNoSeg(seg, paramD)
                        typis.append(typi)
                elif len(segs[q-1]) == 0:
                    typi = 14
                    paramD = stopProtNoSeg(seg, paramD)
                    typis.append(typi)                   

            px = 1
        p += 2
        q += 1

    translated_nt = "".join(paramD['nt_aa'])
    translated_prot = "".join(paramD['prot_aa'])
    raw_nt = "".join(paramD['codons_nt'])
    raw_prot = "".join(paramD['codons_prot'])
    assert len(paramD['rem_nt']) == 0
    assert translated_prot.replace("-", "") == cropped_prot_seq, (translated_prot, cropped_prot_seq)
    assert raw_nt.replace("-", "") == cropped_nt_for_nt_seq
    assert paramD['posi_nt'] == len(cropped_nt_for_nt_seq)
    return ([typis, translated_nt, translated_prot, raw_nt, raw_prot,
             paramD['skipped']])


def processPAF(paffile, stem, minlen,
               aant_seq, prot_seq, nt_seq,
               posiD, rD_nt, rD_aa):
    # Read the PAF file
    paf = pd.read_csv(paffile, sep="\t", header=None,
                      names=ut_functions.getMiniProtColumns(),
                      comment="#")

    if len(paf) != 0:
        paf['Score'] = paf['AS'].str.split(":").str.get(2).astype(int)
        paf = paf.sort_values('Score', ascending=False)
        paf.index = np.arange(len(paf))
        # Iterate through the lines
        for ind in paf.index.values:
            row = paf.loc[ind]
            nmn = int(row['Number_of_matching_nucleotides'])
            if nmn >= minlen:
                rD_aa.setdefault(stem, dict())
                rD_nt.setdefault(stem, dict())
                posiD.setdefault(stem, dict())
                cs_string = row['cs']
                cont_start = row['Contig_start_coordinate']
                cont_end = row['Contig_end_coordinate']
                ref_prot_start = row['Protein_start_coordinate']
                ref_prot_end = row['Protein_end_coordinate']
                strand = row['Strand']
                prot_id = row['Protein_sequence_name']

                R = align_from_cs(cs_string,
                                  aant_seq,
                                  nt_seq,
                                  prot_seq,
                                  ref_prot_start,
                                  ref_prot_end,
                                  cont_start,
                                  cont_end,
                                  strand,
                                  paffile)
                c_aa_ali = R[1]
                c_skipped = R[5]
                if strand == "-":
                    nt_cropped = Bio.Seq.reverse_complement(
                        nt_seq[cont_start:cont_end])
                else:
                    nt_cropped = nt_seq[cont_start:cont_end]
                nt_arr = np.array(list(nt_cropped))
                mask = np.ones(len(nt_arr), dtype=bool)
                mask[c_skipped] = False
                c_nt_recreate = cleanIUPAC("".join(nt_arr[mask]))
                c_aa_recreate = Bio.Seq.translate(c_nt_recreate)
                c_orig = "".join(c_aa_ali.replace("-", ""))
                assert c_aa_recreate == c_orig, (c_nt_recreate, c_aa_recreate,
                                                 c_orig,
                                                 cs_string, R[2], paffile)
                posiD[stem][ind] = R + [prot_id, strand, cont_start, cont_end]
                rD_aa[stem][ind] = c_aa_ali.replace("-", "")
                rD_nt[stem][ind] = c_nt_recreate
    return (posiD, rD_aa, rD_nt)


def runMiniProt(outdir, roundi,
                ntfile, protfile,
                stemnt, stemprot,
                settings='-S -k6 -M0 -F6 -n6 -l6 --outs 0.01 --outc 0.01'):
    statement = f'''miniprot {settings} \
                    {ntfile} \
                    {protfile} 2>{outdir}/{roundi}_{stemnt}_{stemprot}.log \
                    > {outdir}/{roundi}_{stemnt}_{stemprot}.paf'''
    os.system(statement)


def miniProtAllAll(ntfasta, protfasta, outdir, 
                   settings='-S -k6 -M0 -F6 -n6 -l6 --outs 0.01 --outc 0.01',
                   roundi=1):
    # Compare every sequence in ntfasta with every sequence in protfasta
    # using miniprot

    # Output files need to be 
    Fnt = ut_functions.FastaToDict(ntfasta, spliton="|")
    Faa = ut_functions.FastaToDict(protfasta, spliton="|")
    for key, val in Fnt.items():
        out = open(f"{outdir}/nt_{key}.fasta", "w")
        out.write(f">{key}\n{val}\n")
        out.close()
    for key, val in Faa.items():
        out = open(f"{outdir}/aa_{key}.fasta", "w")
        out.write(f">{key}\n{val}\n")
        out.close()

    for i, ntkey in enumerate(Fnt):
        if i % 100 == 0:
            print(i, len(Fnt))
        for protkey in Faa:
            if not os.path.exists(f"{outdir}/{ntkey}/{roundi}_{ntkey}_{protkey}.paf"):
                try:
                    os.mkdir(f"{outdir}/{ntkey}")
                except:
                    pass
                runMiniProt(f"{outdir}/{ntkey}",
                            roundi,
                            f"{outdir}/nt_{ntkey}.fasta",
                            f"{outdir}/aa_{protkey}.fasta", ntkey, protkey,
                             settings=settings)


def cleanAllAlloutput(ntfasta, protfasta, aantfasta, outdir,
                      roundi=1):
    # Read the sequences for nt, aant and prot
    Fnt = ut_functions.FastaToDict(ntfasta, spliton="|")
    Fprot = ut_functions.FastaToDict(protfasta, spliton="|")
    Faant = ut_functions.FastaToDict(aantfasta, spliton="|")
    resiD = dict()
    for j, ntkey in enumerate(Fnt):
        if j % 100 == 0:
            print(j)
        nt_seq = Fnt[ntkey]
        i = 0
        for protkey in Fprot:
            # !
            paffile = f'{outdir}/{ntkey}/{roundi}_{ntkey}_{protkey}.paf'
            if os.path.exists(paffile):
                paf = pd.read_csv(paffile, sep="\t", header=None,
                                  names=ut_functions.getMiniProtColumns(),
                                  comment="#")
            else:
                paf = ""
            if len(paf) != 0:
                indi = f"{ntkey}!{protkey}"
                posiD = dict()
                rD_aa = dict()
                rD_nt = dict()
                aant_seq = Faant[protkey]
                prot_seq = Fprot[protkey]
                posiD, rD_aa, rD_nt = processPAF(paffile,
                                                 indi,
                                                 30,
                                                 aant_seq,
                                                 prot_seq,
                                                 nt_seq,
                                                 posiD,
                                                 rD_nt,
                                                 rD_aa)
                resiD.update(posiD)
                i += 1
    return (resiD)


def allallFindLongest(resiD, ntfasta, refset=None, common=False, comb=False):
    nts = set([x.split("!")[0] for x in resiD])
    maxiD = dict()
    commonD = dict()
    combD = dict()
    ntD = ut_functions.FastaToDict(ntfasta)
    for nt in nts:
        maxiD[nt] = dict()
        maxiD[nt]['length'] = 0
        maxiD[nt]['protein'] = 0
        commonD[nt] = dict()
        commonD[nt]['count'] = 0
        combD[nt] = dict()
        combD[nt]['count'] = 0
        combD[nt]['length'] = 0
        combD[nt]['comb'] = 0

    if not refset:
        refset = set([x.split("!")[1] for x in resiD])

    for key, D in resiD.items():
        nt, aa = key.split("!")
        if aa in refset:
            alllens = []
            for ind, arr in D.items():
                seq = arr[1].replace("-", "")
                leni = len(seq)
                alllens.append(leni)
                nrows = len(D)
                if maxiD[nt]['length'] < leni:
                    maxiD[nt]['length'] = leni
                    maxiD[nt]['sequence'] = seq
                    maxiD[nt]['cont_start'] = arr[8]
                    maxiD[nt]['cont_end'] = arr[9]
                    maxiD[nt]['strand'] = arr[7]
                    maxiD[nt]['skip'] = arr[5]
                    maxiD[nt]['protein'] = aa
                    seq_orig = ntD[nt]
                    seq_orig_c = seq_orig[arr[8]:arr[9]]

                    if arr[7] == "-":
                        seq_orig_c = Bio.Seq.reverse_complement(
                            seq_orig_c)
                    nt_arr = np.array(list(seq_orig_c))
                    mask = np.ones(len(nt_arr), dtype=bool)
                    skippi = arr[5]
                    mask[skippi] = False
                    c_nt_recreate = "".join(nt_arr[mask])
                    c_nt_recreate = cleanIUPAC(c_nt_recreate)
                    c_nt_trans = Bio.Seq.translate(c_nt_recreate)
                    assert c_nt_trans == seq, (c_nt_trans, seq)
            if nrows > commonD[nt]['count']:

                for ind, arr in D.items():
                    seq = arr[1].replace("-", "")
                    leni = len(seq)
                    nrows = len(D)
                    commonD[nt][ind] = dict()
                    commonD[nt][ind]['length'] = leni
                    commonD[nt][ind]['sequence'] = seq
                    commonD[nt][ind]['cont_start'] = arr[8]
                    commonD[nt][ind]['cont_end'] = arr[9]
                    commonD[nt][ind]['strand'] = arr[7]
                    commonD[nt][ind]['skip'] = arr[5]
                    commonD[nt][ind]['protein'] = aa
                    seq_orig = ntD[nt]
                    seq_orig_c = seq_orig[arr[8]:arr[9]]

                    if arr[7] == "-":
                        seq_orig_c = Bio.Seq.reverse_complement(
                            seq_orig_c)
                    nt_arr = np.array(list(seq_orig_c))
                    mask = np.ones(len(nt_arr), dtype=bool)
                    skippi = arr[5]
                    mask[skippi] = False
                    c_nt_recreate = "".join(nt_arr[mask])
                    c_nt_recreate = cleanIUPAC(c_nt_recreate)
                    c_nt_trans = Bio.Seq.translate(c_nt_recreate)
                    assert c_nt_trans == seq, (c_nt_trans, seq)
                    commonD[nt]['count'] = nrows
            if comb:
                totlen = np.sum(alllens)
                if totlen > combD[nt]['comb']:
                    for ind, arr in D.items():
                        seq = arr[1].replace("-", "")
                        leni = len(seq)
                        nrows = len(D)
                        combD[nt][ind] = dict()
                        combD[nt][ind]['length'] = leni
                        combD[nt][ind]['sequence'] = seq
                        combD[nt][ind]['cont_start'] = arr[8]
                        combD[nt][ind]['cont_end'] = arr[9]
                        combD[nt][ind]['strand'] = arr[7]
                        combD[nt][ind]['skip'] = arr[5]
                        combD[nt][ind]['protein'] = aa
                        combD[nt]['comb'] = totlen
                        combD[nt]['length'] = totlen / nrows
                        combD[nt]['count'] = nrows
                        seq_orig = ntD[nt]
                        seq_orig_c = seq_orig[arr[8]:arr[9]]

                        if arr[7] == "-":
                            seq_orig_c = Bio.Seq.reverse_complement(
                                seq_orig_c)
                        nt_arr = np.array(list(seq_orig_c))
                        mask = np.ones(len(nt_arr), dtype=bool)
                        skippi = arr[5]
                        mask[skippi] = False
                        c_nt_recreate = "".join(nt_arr[mask])
                        c_nt_recreate = cleanIUPAC(c_nt_recreate)
                        c_nt_trans = Bio.Seq.translate(c_nt_recreate)
                        assert c_nt_trans == seq, (c_nt_trans, seq)
    if common:
        return (commonD)
    elif comb:
        return (combD)
    else:
        return (maxiD)
