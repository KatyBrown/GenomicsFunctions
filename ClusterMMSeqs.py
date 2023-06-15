import subprocess
import pandas as pd
import os
import ut_functions


def makeDB(infile, dbout, statementout, sout, slog):

    statementdb = ['mmseqs', 'createdb',
                   infile, dbout]
    statementout.write("%s\n" % " ".join(statementdb))
    subprocess.run(statementdb,
                   stdout=sout,
                   stderr=slog,
                   check=True)

    return (dbout)


def cluster(db, clustout, tempout,
            statementout, sout, slog,
            min_sim, min_cov, cov_mode):

    statementclust = ['mmseqs', 'cluster', '--cluster-reassign',
                      '--min-seq-id', '%.2f' % min_sim,
                      '-c', '%.2f' % min_cov,
                      '--cov-mode', "%i" % cov_mode,
                      db, clustout, tempout]
    statementout.write("%s\n" % " ".join(statementclust))
    subprocess.run(statementclust,
                   stdout=sout,
                   stderr=slog,
                   check=True)


def makeTxt(clusters, db, outfile, statementout, sout, slog):
    
    statementtxt = ['mmseqs', 'createtsv',
                    db, clusters, outfile]
    statementout.write("%s\n" % " ".join(statementtxt))
    subprocess.run(statementtxt,
                   stdout=sout,
                   stderr=slog,
                   check=True)    



def makeTablesFastas(lookup, txtfile, orig_fasta, rep_fasta_out, tab_out, fastadir):
    try:
        os.mkdir(fastadir)
    except:
        pass
    F = ut_functions.FastaToDict(orig_fasta)
    inds = pd.read_csv(lookup, sep="\t", header=None, index_col=0)
    indD = dict(zip(inds.index.values, inds[1]))
    tab = pd.read_csv(txtfile, sep="\t", header=None)
    tab.columns = ['Rep_Name', 'Sequence_Ind']
    tab['Cluster_Members'] = [indD[x] for x in tab['Sequence_Ind']]
    clustered = tab.groupby('Rep_Name')
    outr = open(rep_fasta_out, "w")
    outt = open(tab_out, "w")

    outt.write("group_id\tgroup_member\tgroup_rep\tgroup_size\n")
    gn = 1
    for group in set(tab.Rep_Name):
        groupt = clustered.get_group(group)
        group_id = "g%i" % (gn)
        group_size = len(groupt)
        out = open("%s/%s.fasta" % (fastadir, group.replace("|", "_").replace("/", "")[0:100]), "w")
        for member in groupt['Cluster_Members']:
            outt.write("%s\t%s\t%s\t%s\n" % (group_id, member, group,
                                             group_size))
            out.write(">%s_%s\n%s\n" % (group_id, member, F[member]))
        out.close()
        outr.write(">%s|%s\n%s\n" % (group_id, group, F[group]))
        gn += 1
    outr.close()
    outt.close()