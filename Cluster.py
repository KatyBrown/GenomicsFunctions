#!/usr/bin/env python3
import ut_functions
import os
import shutil
import networkx as nx
import sys
sys.path.insert(0, '/home/katy/skbio/scikit-bio')
import skbio
from skbio import alignment
import copy
import sys
import CIAlign.utilityFunctions as utilityFunctions
import CIAlign.similarityMatrix as similarityMatrix
import sklearn
import numpy as np
import sklearn.manifold
import sklearn.cluster
import string
import pathlib
import math

def SWalign(seq1, seq2):
    seq1 = skbio.Protein(seq1)
    seq2 = skbio.Protein(seq2)
    
    BLOSUM = ut_functions.BLOSUM()
    try:
        return alignment.local_pairwise_align_ssw(seq1, seq2, substitution_matrix=BLOSUM)
    except:
        return None


def SWaligntest(seq1, seq2, nam1, nam2, min_score, min_perc,
                rscore=False):
    '''
    Test if two sequences are more similar than min_score and min_perc
    using the Smith Waterman algorithm
    
    If one sequence overlaps the end of the other, min_perc doesn't
    need to be satisfied but if one sequence is internal to the other
    it does.
    '''
    orig_len_1 = len(seq1)
    orig_len_2 = len(seq2)

    # put the longer sequeunce first
    if orig_len_1 > orig_len_2:
        seq1, seq2 = seq2, seq1
        nam1, nam2 = nam2, nam1
    orig_len_1 = len(seq1)
    orig_len_2 = len(seq2)
    
    seq1 = skbio.Protein(seq1)
    seq2 = skbio.Protein(seq2)

    BLOSUM = ut_functions.BLOSUM()

    try:
        ali = skbio.alignment.local_pairwise_align_ssw(seq1, seq2, substitution_matrix=BLOSUM)
        score = ali[1]
        adj_score = score / orig_len_1
    except:
        adj_score = 0
        ali = None

    if adj_score > min_score:
        p1, p2 = ali[2][0], ali[2][1]

        # Does one sequence overlap the end of the other?
        one_left = p1[0] in range(0, 2)
        two_left = p2[0] in range(0, 2)

        one_right = p1[1] in range(orig_len_1 - 2, orig_len_1 + 2)
        two_right = p2[1] in range(orig_len_2 - 2, orig_len_2 + 2)
        
        # the percentage of the shorter sequence which aligns well
        perc = (p1[1] - p1[0]) / orig_len_1
        
        # either the aligned region is at the end of one sequence
        # or the shorter sequence is completely inside the longer sequence with high identity
        if ((one_left and two_right) or (one_right and two_left)) or perc > min_perc:
            if rscore:
                return adj_score
            return True
        else:
            if rscore:
                return 0
            return False
    else:
        if rscore:
            return 0
        return False
    
    
def clusterCyclesSW(infile, outfile_prefix, outpath, final_out, batchsize=200,
                    min_score=6, min_perc=0.95,
                    min_length=50, keep_symbol="", calign_path="."):
    '''
    Cluster a set of ORFS using the Smith Waterman algorithm

    Parameters
    ----------
    infile : TYPE
        DESCRIPTION.
    outfile_prefix : TYPE
        DESCRIPTION.
    outpath : TYPE
        DESCRIPTION.
    final_out : TYPE
        DESCRIPTION.
    batchsize : TYPE, optional
        DESCRIPTION. The default is 200.
    min_score : TYPE, optional
        DESCRIPTION. The default is 6.
    min_perc : TYPE, optional
        DESCRIPTION. The default is 0.99.
    min_length : TYPE, optional
        DESCRIPTION. The default is 50.
    keep_symbol : TYPE, optional
        DESCRIPTION. The default is "".
    calign_path : TYPE, optional
        DESCRIPTION. The default is ".".

    Returns
    -------
    None.

    '''
    # open an empty table to track how many sequences are in each cluster
    count_tab = open("%s/%s_clusters.tsv" % (outpath, outfile_prefix), "w")
    # open the file which will be clustered
    F = ut_functions.FastaToDict(infile)
    currentF = dict()
    # keep only the consensus sequences / only the sequences with the "keep"
    # symbol
    for nam in F:
        if keep_symbol in nam or "cons" in nam:
            currentF[nam] = F[nam]

    # make an output file with only the sequences which will be kept
    orig_out = "%s/%s_1.fasta" % (outpath, outfile_prefix)
    out = open(orig_out, "w")
    for key in currentF:
        out.write(">%s\n%s\n" % (key, currentF[key]))
    out.close()

    # keep track of the length of the previous and the current list
    # of sequences
    pfasta_len = len(currentF) + 1
    cfasta_len = len(currentF)
    
    p = 1
    # while the number of sequences is decreasing
    while pfasta_len > cfasta_len:
        # split the fasta file into batches
        # the batches need to be quite big or nothing will cluster
        nbatches = math.ceil(len(currentF) / batchsize)
        F_L = ut_functions.batchFasta(currentF, nbatches)
        
        n = 1
        Cdict = dict()
        # for each batch of sequences
        for k, batch in enumerate(F_L):
            nams = []
            seqs = []
            for key in batch:
                nams.append(key)
                seqs.append(batch[key])
            
            # sort in ascending sequence length
            z = list(zip(nams, seqs))
            z.sort(key = lambda t: len(t[1]), reverse=False)
            nams = [o[0] for o in z]
            seqs = [o[1] for o in z]
            G = nx.Graph()
            # make every pairwise comparison
            # build a network with connections between sequences which are
            # sufficiently similar
            for i, seq1 in enumerate(seqs):
                nam1 = nams[i]
                G.add_node(nam1)
                for j, seq2 in enumerate(seqs[i+1:]):
                    nam2 = nams[j+i+1]
                    G.add_node(nam2)
                    # test identity using the Smith-Waterman algorithm
                    ali_test = SWaligntest(seq1, seq2,
                                           nam1, nam2,
                                           min_score, min_perc)
                    if ali_test:
                        G.add_edge(nam1, nam2)
    
            # extract the connected groups from the network
            ccs = nx.connected_components(G)
            cdict = dict()
            # for each cluster
            for cc in ccs:
                # put all the sequences from the cluster into an output file
                out = open("%s/%s_%s_%i.fasta" % (outpath, outfile_prefix, n+1, p), "w")
                m = 0
                for val in cc:
                    out.write(">%s\n%s\n" % (val, currentF[val]))
                    m += 1
                out.close()
                # if there is only one sequence in the cluster
                if m == 1:
                    x = ut_functions.FastaToDict("%s/%s_%s_%i.fasta" % (outpath,
                                                                        outfile_prefix,
                                                                     n+1, p))
                    if len(x) != 0:
                        # copy the sequence into the output file
                        out2 = open("%s/%s_%s_%i_consensus.fasta" % (outpath, outfile_prefix,
                                                                  n+1, p), "w")
                        out2.write(">%s_%s_%i_cons\n%s\n" % (outfile_prefix, n+1, p, list(x.values())[0]))
                        out2.close()
                        # keep track of how many sequences this represents
                        count_tab.write("%s\t%s_%s_%i\n" % (list(x.keys())[0], outfile_prefix, n+1, p))
                        cdict.update(ut_functions.FastaToDict("%s/%s_%s_%i_consensus.fasta" % (outpath, outfile_prefix, n+1, p)))
                        n += 1
                elif m > 1:
                        # if there is more than one sequence in the cluster, check
                        # if it should be subdivided by reclustering using a different
                        # algorithm and see if it forms one cluster or more than one

                        # align the current set of sequences
                        statement = '''mafft --localpair %s/%s_%s_%i.fasta \
                        > %s/%s_%s_%i_ali_test.fasta''' % (outpath, outfile_prefix, n+1, p,
                        outpath, outfile_prefix, n+1, p)
                        os.system(statement)

                        arr, nams = utilityFunctions.FastaToArray(
                            "%s/%s_%s_%i_ali_test.fasta" % (outpath, outfile_prefix, n+1, p),
                            outfile_stem=None)
                        
                        # calculate the pairwise similarity between the current set of sequecnes
                        ident = similarityMatrix.calculateSimilarityMatrix(arr, nams, outfile="%s/%s_%s_%i_mat.txt" % (
                            outpath, outfile_prefix, n+1, p))
                        
                        # convert the distance between the sequences into a set of multi-dimensional scaling co-ordinates
                        M = sklearn.manifold.MDS(5, dissimilarity='precomputed')
                        m = M.fit_transform(1-ident)
                        
                        # apply agglomerative clustering to the sequences
                        c = sklearn.cluster.AgglomerativeClustering(n_clusters=None, distance_threshold=0.25).fit_predict(m)
                        
                        # for each subcluster
                        for i in range(max(c) + 1):
                            out = open("%s/%s_%s_%i.fasta" % (outpath, outfile_prefix, n+1, p), "w")
                            m = 0
                            t = np.array(nams)[c == i]
                            # put the current subcluster into a file
                            for val in t:
                                out.write(">%s\n%s\n" % (val, currentF[val]))
                                m += 1
                            out.close()
                            # if there is one sequence in the subcluster, just copy
                            # it to make the consensus
                            if m == 1:
                                x = ut_functions.FastaToDict("%s/%s_%s_%i.fasta" % (outpath,
                                                                                    outfile_prefix,
                                                                                 n+1, p))
                                if len(x) != 0:
                                    out2 = open("%s/%s_%s_%i_consensus.fasta" % (outpath, outfile_prefix,
                                                                              n+1, p), "w")
                                    out2.write(">%s_%s_%i_cons\n%s\n" % (outfile_prefix, n+1, p, list(x.values())[0]))
                                    out2.close()
                                    count_tab.write("%s\t%s_%s_%i\n" % (list(x.keys())[0], outfile_prefix, n+1, p))
                                    cdict.update(ut_functions.FastaToDict("%s/%s_%s_%i_consensus.fasta" % (outpath, outfile_prefix, n+1, p)))
                                    n += 1
                            elif m > 1:
                                # if there is more than one sequence in the subcluster,
                                # align it and generate a consensus
                                statement = '''mafft --localpair %s/%s_%s_%i.fasta \
                                > %s/%s_%s_%i_ali.fasta''' % (outpath, outfile_prefix, n+1, p,
                                outpath, outfile_prefix, n+1, p)
                                os.system(statement)
                                # parse the alignment and generate a consensus
                                statement = '''CIAlign --infile %s/%s_%s_%i_ali.fasta \
                                --remove_min_length %s --outfile_stem %s/%s_%s_%i \
                                --remove_insertions --crop_ends --remove_short \
                                --make_consensus --consensus_type majority_nongap \
                                --consensus_name %s_%s_%i_cons''' % (outpath, outfile_prefix, n+1,
                                p, min_length, outpath, outfile_prefix, n+1, p, outfile_prefix, n+1, p)
                                ut_functions.writeCommand(statement, "%s/%s" % (outpath, outfile_prefix))
                                os.system(statement)
                                ###### !!!!!!
                                # sometimes there is no sequence left after the previous step
                                if os.path.exists("%s/%s_%s_%i_consensus.fasta"  %(outpath, outfile_prefix, n+1, p)):
                                    x = ut_functions.FastaToDict("%s/%s_%s_%i_cleaned.fasta" % (outpath,
                                                                                               outfile_prefix,
                                                                                               n+1, p))
                                    for key in x:
                                        # keep count of the number of sequences represented
                                        count_tab.write("%s\t%s_%s_%i\n" % (key, outfile_prefix, n+1, p))
                                    if os.path.exists("%s/%s_%s_%i_removed.txt" % (outpath, outfile_prefix, n+1, p)):
                                        # keep track of which sequences have been removed because they didn't match well
                                        lines = [line.strip() for line in open("%s/%s_%s_%i_removed.txt" % (outpath, outfile_prefix, n+1, p)).readlines()]
                                        for line in lines:
                                            count_tab.write("%s\tremoved\n" % (line))
                                    cdict.update(ut_functions.FastaToDict("%s/%s_%s_%i_consensus.fasta" % (outpath, outfile_prefix, n+1, p)))
                                    n += 1
                
            Cdict.update(cdict)
        p += 1
        # Cdict contains all the consensus sequences
        out = open("%s/%s_%i.fasta" % (outpath, outfile_prefix, p), "w")
        for val in Cdict:
            out.write(">%s\n%s\n" % (val, Cdict[val]))
        out.close()
        
        pfasta_len = len(currentF)
        cfasta_len = len(Cdict)
        currentF = copy.copy(Cdict)
    count_tab.close()


def clusterCyclesCDHit(infile, outfile_prefix, outpath, final_out, C=0.95, L=0.25,
                  S=0.25, M=50, keep_symbol="", calign_path="."):
    '''
    Run cd-hit on a fasta file and generate a consensus sequence for each cluster
    Keep repeating this until the number of clusters doesn't change with additional runs
    '''
    F = ut_functions.FastaToDict(infile)
    res = dict()
    for nam in F:
        if keep_symbol in nam or "cons" in nam:
            res[nam] = F[nam]
    orig_out = "%s/%s_1.fasta" % (outpath, outfile_prefix)
    out = open(orig_out, "w")
    for key in res:
        out.write(">%s\n%s\n" % (key, res[key]))
    out.close()
    
    j = 2
    currentL = len(F)
    currentF = orig_out
    while True:
        unclusteredL = len(ut_functions.FastaToDict(currentF))
        out2 = "%s/%s_%i.fasta" % (outpath, outfile_prefix, j)
        # run the clustering command
        statement = """cd-hit -i %s -o %s -M 0 -d 1000 \
        -T 1 -c %s -aL %s -aS %s -G 0""" % (currentF, out2, C, L, S)
        ut_functions.writeCommand(statement, "%s/%s" % (outpath, outfile_prefix))
        os.system(statement)
        F = ut_functions.FastaToDict(currentF)
        clusterfile = "%s/%s_%i.fasta.clstr" % (outpath, outfile_prefix, j)
        lines = [line.strip() for line in open(clusterfile).readlines()]
        current = []
        currents = []
        # make a list of sequences in each cluster
        for line in lines:
            if not line.startswith(">"):
                current.append(line.split(">")[1].split("...")[0])
            else:
                currents.append(current)
                current = []
        currents.append(current)
        currents = currents[1:]
        currentL = len(currents)
        if unclusteredL == currentL:
            break
        cdict = dict()
        for i, current in enumerate(currents):
            out = open("%s/%s_%s_%i.fasta" % (outpath, outfile_prefix, i + 1, j), "w")
            for nam in current:
                if keep_symbol in nam or "_cons" in nam:
                    out.write(">%s\n%s\n" % (nam, F[nam]))
            out.close()
            if len(current) == 1:
                x = ut_functions.FastaToDict("%s/%s_%s_%i.fasta" % (outpath,
                                                                    outfile_prefix,
                                                                 i+1, j))
                if len(x) != 0:
                    out2 = open("%s/%s_%s_%i_consensus.fasta" % (outpath, outfile_prefix,
                                                              i+1, j), "w")
                    out2.write(">%s_%s_%i_cons\n%s\n" % (outfile_prefix, i+1, j, list(x.values())[0]))
                    out2.close()
            elif len(current) > 1:
                    statement = '''mafft --localpair %s/%s_%s_%i.fasta \
                    > %s/%s_%s_%i_ali.fasta''' % (outpath, outfile_prefix, i+1, j,
                    outpath, outfile_prefix, i+1, j)
                    os.system(statement)
                    statement = '''CIAlign --infile %s/%s_%s_%i_ali.fasta \
                    --remove_min_length %s --outfile_stem %s/%s_%s_%i \
                    --remove_insertions --crop_ends --remove_short \
                    --make_consensus --consensus_type majority_nongap \
                    --consensus_name %s_%s_%i_cons''' % (outpath, outfile_prefix, i+1,
                    j, M, outpath, outfile_prefix, i+1, j, outfile_prefix, i+1, j)
                    ut_functions.writeCommand(statement, "%s/%s" % (outpath, outfile_prefix))
                    os.system(statement)    
            cdict.update(ut_functions.FastaToDict("%s/%s_%s_%i_consensus.fasta" % (outpath, outfile_prefix,
                                                                                i+1, j)))
        out = open("%s/%s_cleaned_%i.fasta" % (outpath, outfile_prefix, j), "w")
        for nam in cdict:
            out.write(">%s\n%s\n" % (nam, cdict[nam]))
        out.close()
        currentF = "%s/%s_cleaned_%i.fasta" % (outpath, outfile_prefix, j)
        j += 1
    final = "%s/%s_cleaned_%i.fasta" % (outpath, outfile_prefix, j-1)
    shutil.copy(final, final_out)


def expandClusters(infile_fasta, infile_clusters,
                   infile_dir, outfile_dir, calign_path):
    clusters = open(infile_clusters).readlines()
    
    F = ut_functions.FastaToDict(infile_fasta)
    D = dict()
    D2 = dict()
    bigF = dict()
    for line in clusters:
        line = line.strip().split("\t")
        if "_samples_" in line[0] or "~" in line[0]:
            D[line[0]] = line[1]
            D2.setdefault("%s_cons" % line[1], [])
            D2['%s_cons' % line[1]].append(line[0])
            if line[1] != "removed":
                bigF.update(
                        ut_functions.FastaToDict("%s/%s_consensus.fasta" % (infile_dir, line[1])))

        elif "emov" not in line[0]:

            origs = D2[line[0]]
            for orig in origs:
                D[orig] = line[1]
                D2.setdefault("%s_cons" % line[1], [])
                D2["%s_cons" % line[1]].append(orig)
                if line[1] != "removed":
                    bigF.update(
                            ut_functions.FastaToDict("%s/%s_consensus.fasta" % (infile_dir, line[1])))
    
    D_rev = dict()
    for key, val in D.items():
        D_rev.setdefault(val, [])
        D_rev[val].append(key)
    for i, key in enumerate(D_rev):
        if key != "removed":
            out = open("%s/%s_cleaned2.fasta" % (infile_dir, key), "w")
            f = ut_functions.FastaToDict("%s/%s_consensus.fasta" % (infile_dir, key))
            if len(f) != 0:
                cons = list(f.values())[0]
                for val in D_rev[key]:
                    seq = F[val]
                    C = SWalign(seq, cons)
                    if C is not None:
                        coords_cons = C[2][1]    
                        cons_subseq = "".join([x.decode() for x in C[0][1].values])
                        new_subseq = "".join([x.decode() for x in C[0][0].values])
                    else:
                        C = SWalign(cons, seq)
                        coords_cons = C[2][0]
    
                        cons_subseq = "".join([x.decode() for x in C[0][0].values])
                        new_subseq = "".join([x.decode() for x in C[0][1].values])
    
                    cons_array = np.array(list(cons_subseq)) != "-"
                    new_array = np.array(list(new_subseq))[cons_array]
                    
                    new_cleaned = "".join(new_array)
                    buffer_start = coords_cons[0]
                    buffer_end = len(cons) - (buffer_start + sum(cons_array))
                    new_cleaned_buffered = "-" * buffer_start + new_cleaned + "-" * buffer_end
                    
                    out.write(">%s\n%s\n" % (val, new_cleaned_buffered))
            out.close()
            os.system("""CIAlign --infile %s/%s_cleaned2.fasta \
                      --outfile_stem %s/%s \
                      --remove_insertions \
                      --insertion_min_size 1 \
                      --insertion_max_size 20 \
                      --make_consensus \
                      --consensus_type majority_nongap""" % (
                      infile_dir, key, outfile_dir, key))    


def getRanges(boolarray, runlength):
    binarray = boolarray.astype(int)
    diffs = np.diff(binarray)
    range_starts = np.where(diffs == 1)[0]
    range_ends = np.where(diffs == -1)[0]

    if len(range_starts) == 0:
        range_starts = np.array([0])
    if len(range_ends) == 0:
        range_ends = np.array([len(binarray)])
    
    if min(range_starts) > min(range_ends):
            range_starts = np.append([0], range_starts)
    if max(range_starts) > max(range_ends):
            range_ends = np.append(range_ends, [len(binarray)])

    coords = zip(range_starts, range_ends)
    clist = []
    for c in coords:
        if c[1] - c[0] > runlength:
            clist.append(c)
    return (clist)


def parseClusters(infile, outfile, minlength, runlength, overlap, calign_path):
    k = 0
    m = 1
    letters = list(string.ascii_lowercase)

    arr, nams = utilityFunctions.FastaToArray(infile, outfile_stem=None)

    F = ut_functions.FastaToDict(infile)
    binmatrix = arr != "-"
    arrlen = np.shape(arr)[1]
    if arrlen > 1:
        clist = []
        for i, row in enumerate(binmatrix):
            clist += getRanges(row, runlength)
        intervals = set()
        for windowstart in range(arrlen):
            for windowsize in range(2, overlap+1):
                windowend = windowstart + windowsize
                covered = 0
                for coords in clist:
                    if windowstart > coords[0] and windowend < coords[1]:
                        covered += 1
                if covered == 0:
                    intervals = intervals | set(list(range(windowstart, windowend)))
        bool_intervals = np.isin(np.arange(0, arrlen), list(intervals))
        ranges_intervals = getRanges(bool_intervals, 1)
        isend = [0 in x or arrlen in x for x in ranges_intervals]
        ranges_intervals_nonend = np.array(ranges_intervals)[np.invert(isend)]
        s = 0
        o = 0
        if len(ranges_intervals_nonend) > 0:
            f = infile.split("/")[-1], ranges_intervals_nonend
            l = 0
            for r in f[1]:
                if r[0] - s > minlength:
                    outf = outfile.replace("_cleaned.fasta", "_uncleaned_%s.fasta" % letters[l])
                    out = open(outf, "w")
                    for nam in nams:
                        out.write(">%s\n%s\n" % (nam, F[nam][s:r[0]]))
                    out.close()
                    outstem = "%s_%s" % (outfile.replace("_cleaned.fasta", ""), letters[l])
                    pref = outstem.split("/")[-1]
                    statement = '''CIAlign --infile %s \
                    --remove_min_length %s --outfile_stem %s \
                    --remove_short --make_consensus \
                    --consensus_type majority_nongap \
                    --consensus_name %s_cons''' % (outf,
                                                   minlength, outstem,
                                                   pref)
                    os.system(statement)       
                    o += 1
                    l += 1
                s = r[1]
            k += 1
            if arrlen - s > minlength and s != 0:
                outf = outfile.replace("_cleaned.fasta", "_uncleaned_%s.fasta" % letters[l])
                out = open(outf, "w")
                for nam in nams:
                    out.write(">%s\n%s\n" % (nam, F[nam][s:arrlen]))
                out.close()
                outstem = "%s_%s" % (outfile.replace("_cleaned.fasta", ""), letters[l])
                pref = outstem.split("/")[-1]
                statement = '''CIAlign --infile %s \
                --remove_min_length %s --outfile_stem %s \
                --remove_short --make_consensus \
                --consensus_type majority_nongap \
                --consensus_name %s_cons''' % (outf,
                                               minlength, outstem, pref)
                os.system(statement)
                o += 1
            pathlib.Path(outfile).touch()
            
        else:
            outf = outfile.replace("_cleaned.fasta", "_uncleaned.fasta")
            shutil.copy(infile, outf)
            outstem = outfile.replace("_cleaned.fasta", "")
            pref = outstem.split("/")[-1]
            statement = '''CIAlign --infile %s \
                --remove_min_length %s --outfile_stem %s \
                --remove_short --make_consensus \
                --consensus_type majority_nongap \
                --consensus_name %s_cons''' % (outf, minlength,
                                                outstem, pref)
            os.system(statement)
            pathlib.Path(outfile).touch()
            m += 1
