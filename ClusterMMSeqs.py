import pathlib
import subprocess
import shutil


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



