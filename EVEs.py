import glob
import os
import HMMs
import pandas as pd
import numpy as np
import pyfaidx
import pybedtools
import Bio.Seq
import ut_functions


def getGenome(genome_id, outdir):
    statement1 = f"""datasets download --no-progressbar \
                     genome accession {genome_id} --include genome"""
    os.system(statement1)
    statement2 = f"""unzip -d {outdir}/{genome_id} ncbi_dataset.zip;
                     unlink ncbi_dataset.zip"""
    os.system(statement2)
    statement3 = f"""fn=$(ls {outdir}/{genome_id}/ncbi_dataset/data/*/*fna \
                          | cut -d "/" -f2-);
                     ln -s $fn {outdir}/{genome_id}.fna"""
    os.system(statement3)
    pyfaidx.Faidx(f"{outdir}/{genome_id}.fna")
    return (f"{outdir}/{genome_id}.fna")


def translateGenome(genome_id, genome_path, outdir):
    statement = f"""transeq -frame 6 -sequence {genome_path} \
                     -outseq {outdir}/{genome_id}_6f.fasta"""
    os.system(statement)
    pyfaidx.Faidx(f"{outdir}/{genome_id}_6f.fasta")
    return (f"{outdir}/{genome_id}_6f.fasta")


def rmGenome(genome_id, outdir):
    statement = f"""rm -rf {outdir}/{genome_id};
                    rm -rf {outdir}/{genome_id}.fna"""
    os.system(statement)


def getORFs(genome_id, trans_path, outdir):
    out = open(f"{outdir}/{genome_id}_segs.fasta", "w")
    with open(f"{outdir}/{genome_id}_6f.fasta") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                header = line.strip()
                seqi = ""
                i = 1
            else:
                for char in line:
                    if char != "*":
                        seqi += char
                    else:
                        h = header.split(" ")[0]
                        if len(seqi) > 50:
                            out.write(f"{h}_{i}\n{seqi}\n")
                            i += 1
                        seqi = ""
    if len(seqi) != 0:
        out.write(f">{h}_{i}\n{seqi}\n")
    out.close()
    pyfaidx.Faidx(f"{outdir}/{genome_id}_segs.fasta")
    return (f"{outdir}/{genome_id}_segs.fasta")


def runHMMER(orfs, genome_id, outdir, hmm, hmmpref):
    HMMs.runHMMER(orfs,
                  [f"{outdir}/hmmer_{genome_id}_{hmmpref}.tab",
                   f"{outdir}/hmmer_{genome_id}_{hmmpref}.out",
                   f"{outdir}/hmmer_{genome_id}_{hmmpref}.dom"],
                  hmm,
                  f"{outdir}/hmmer_{genome_id}_{hmmpref}.log")


def filterHMMER(genome_id, outdir, thresh):
    hmmer = pd.read_csv(f"{outdir}/{genome_id}_hmmer_combined.tsv",
                        sep="\t")
    hmmer = hmmer[hmmer['score'] >= thresh]
    hmmer = hmmer.sort_values('score', ascending=False)
    gd = hmmer.groupby('query_name').first()
    gd['query_name'] = gd.index.values
    gd.index = np.arange(len(gd))
    gd.to_csv(f"{outdir}/{genome_id}_filt.tsv", sep="\t", index=None)
    return (f"{outdir}/{genome_id}_filt.tsv")


def getPositionsTrans(genome_id, outdir, seqs):
    posis = []
    j = 0
    seqi = ""
    with open(f"{outdir}/{genome_id}_6f.fasta") as infile:
        js = []
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                js.append(j)
                if seqi in seqs:
                    if js[-1] < js[0]:
                        posis.append([header, js[0], js[-1]])
                js = []
                header = line.strip()
                seqi = ""
                i = 1
                j = 0
            else:
                for char in line:
                    if char != "*":
                        seqi += char
                        js.append(j)
                    else:
                        js.append(j)
                        if seqi in seqs:
                            if js[-1] < js[0]:
                                print (header, seqi, js)
                                exit
                            posis.append([header, js[0], js[-1]])
                            i += 1
                        seqi = ""
                        js = []
                    j += 1
    return (posis)


def getPositionsGenome(genome_id, outdir, posis, seqs):
    dat_genome = pyfaidx.Fasta(f"{outdir}/{genome_id}.fna")
    rows = []
    for p in posis:
        pp = p[0].split(" ")[0].replace(">", "")
        px = "_".join(pp.split("_")[:-1])
        frame = int(pp[-1])
        snt = p[1] * 3
        ent = p[2] * 3
        sff = None
        eff = None
        if frame >= 4:
            strand = "-"
            L = dat_genome[px].unpadded_len
            for adj in np.arange(-3, 3):
                sf = L - ent + adj
                ef = L - snt + adj
                if sf > 0 and ef > 0:

                    subseq = Bio.Seq.reverse_complement(dat_genome[px][sf:ef].seq)
                    trans = Bio.Seq.translate(subseq)
                    if trans in seqs:
                        sff = sf
                        eff = ef
        else:
            strand = "+"
            for adj in np.arange(-3, 3):
                sf = snt + adj
                ef = ent + adj
                if sf > 0 and ef > 0:
                    subseq = dat_genome[px][sf:ef].seq
                    trans = Bio.Seq.translate(subseq)
                    if trans in seqs:
                        sff = sf
                        eff = ef
        assert sff and eff
        rows.append([px, sff, eff, pp, len(subseq), strand])
    beddf = pd.DataFrame(rows)
    beddf[1] = beddf[1].astype(int)
    beddf[2] = beddf[2].astype(int)
    bed = pybedtools.BedTool.from_dataframe(beddf).sort()
    ff = bed.getfasta(fi=f"{outdir}/{genome_id}.fna", s=True)
    ff.save_seqs(f"{outdir}/{genome_id}_orfs.fasta")
    return f"{outdir}/{genome_id}_orfs.fasta"


def getRegions(genome_id, filt, outdir):
    filt_tab = pd.read_csv(filt, sep="\t")
    accs = set(filt_tab['query_name'])
    dat_segs = pyfaidx.Fasta(f"{outdir}/{genome_id}_segs.fasta")
    seqs = []
    for acc in accs:
        seq = dat_segs[acc][:]
        seqs.append(seq.seq)
    seqs = set(seqs)
    posis = getPositionsTrans(genome_id, outdir, seqs)
    bedout = getPositionsGenome(genome_id, outdir, posis, seqs)
    return (bedout)


def translateRegions(genome_id, outdir, bedout):
    F = ut_functions.FastaToDict(bedout)
    out = open(f"{outdir}/{genome_id}_orfs_trans.fasta", "w")
    for nam, seq in F.items():
        trans = Bio.Seq.translate(seq)
        assert "*" not in trans
        if len(set(list(trans))) > 12:
            out.write(f">{nam}\n{trans}\n")
    out.close()
    return f"{outdir}/{genome_id}_orfs_trans.fasta"


def runAll(genome_id, outdir, hmmlist, hmmpref, thresh):
    genome = getGenome(genome_id, outdir)
    trans = translateGenome(genome_id, genome, outdir)
    orfs = getORFs(genome_id, trans, outdir)
    if not os.path.exists("f{outdir}/{genome_id}_filt.tsv"):
        for hh, hp in zip(hmmlist, hmmpref):
            runHMMER(orfs, genome_id, outdir, hh, hp)
    HMMs.combineHMMEROutput(
        glob.glob(f"{outdir}/hmmer_{genome_id}_*.dom"),
        f"{outdir}/{genome_id}_hmmer_combined.tsv", suffix=".tsv")
    filt = filterHMMER(genome_id, outdir, thresh)
    hmmer = pd.read_csv(filt,
                        sep="\t")    
    bed_reg = getRegions(genome_id, filt, outdir)
    if len(filt) != 0:
        translateRegions(genome_id, outdir, bed_reg)
    rmGenome(genome_id, outdir)
