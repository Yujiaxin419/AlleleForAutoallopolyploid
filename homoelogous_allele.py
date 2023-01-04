import pysam
import sys
import os

def check_faType(fa):
    faType = 'cds'
    with open(fa, 'r') as fin:
        for line in fin:
            if line.startswith('>') == False:
                lineStr = list(line.rstrip().upper())
                for char in lineStr:
                    if char not in ["A", "G", "T", "C", "N"]:
                        return 'pep'
                return faType

def run_balst(faA, faB, faType, threads, prefix):
    if faType == 'cds':
        builddb_cmd = "makeblastdb -in {} -dbtype nucl".format(faA)
        blast_cmd = "blastn -db {} -num_threads {} \
                     -outfmt 6 -out {}.blast -query {} -evalue 1e-5 \
                     -best_hit_score_edge 0.05 -best_hit_overhang 0.25  -max_target_seqs 1".format(faA, threads, prefix, faB)
    else:
        builddb_cmd = "makeblastdb -in {} -dbtype nucl".format(faA)
        blast_cmd = "blastn -db {} -num_threads {} \
                     -outfmt 6 -out {}.blast -query {} -evalue 1e-5 \
                     -best_hit_score_edge 0.05 -best_hit_overhang 0.25  -max_target_seqs 1".format(faA, threads, prefix, faB)
    print("Blast CMD is : {}".format(builddb_cmd))
    os.system(builddb_cmd)
    print("Blast CMD is : {}".format(blast_cmd))
    os.system(blast_cmd)
    print("Blast done!")

def read_faSize(fa):
    f = pysam.FastaFile(fa)
    faLenDic = dict(zip(f.references, f.lengths))
    return faLenDic

def read_blast(bf, faSize, identCutf, covCutf):
    blastDic = {}
    with open(bf, 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split()
            [qry, ref, ident, alen] , bhit = tmpLst[:4], float(tmpLst[-1])
            cov = float(alen)/float(faSize[qry])*100
            if cov < float(covCutf) or float(ident) < float(identCutf):
                continue
            if qry not in blastDic:
                blastDic[qry] = [("",0)]
            if bhit > blastDic[qry][0][1]:
                blastDic[qry] = [(ref, bhit)]
            elif bhit == blastDic[qry][0][1]:
                blastDic[qry].append((ref,bhit))
            else:
                continue
    return blastDic

def reciprocal_blast(blastDicA, blastDicB):
    homologDic = {}
    for qryA in blastDicA:
        if len(blastDicA[qryA]) == 1:
            refA = blastDicA[qryA][0][0]
            if refA not in blastDicB:
                continue
            qryB = refA
            if len(blastDicB[qryB]) == 1:
                refB = blastDicB[qryB][0][0]
                if qryA == refB:
                    homologDic[qryA] = qryB
            else:  # blastDicB[qryB] = [(geneA,100), (geneB, 100)]
                qbLst = blastDicB[qryB]
                if qryA in qbLst:
                    for subqb in qbLst:
                        qb = subqb[0]
                        homologDic[qryA] = qb
        else:
            for raLst in blastDicA[qryA]:
                ra = raLst[0]
                if ra not in blastDicB:
                    continue
                qryB = ra
                if len(blastDicB[qryB]) == 1:
                    refB = blastDicB[qryB][0][0]
                    if qryA == refB:
                        homologDic[qryA] = qryB
                else:  # blastDicB[qryB] = [(geneA,100), (geneB, 100)]
                    qbLst = blastDicB[qryB]
                    if qryA in qbLst:
                        for subqb in qbLst:
                            qb = subqb[0]
                            homologDic[qryA] = qb
    with open('homoeologous.list', 'w') as fout:
        for ha in homologDic:
            fout.write("{}\t{}\n".format(ha, homologDic[ha]))
    return homologDic

def output_allele(homologDic, allele):
    alleleDic = {}
    with open(allele, 'r') as fin:
        header = fin.readline().rstrip()
        length_allele = len(header.split('\t'))
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            ref = tmpLst[0]
            alleleDic[ref] = line.rstrip()
    with open('Allele.homoeologous.table', 'w') as fout:
        fout.write("{}\t{}\n".format(header, header))
        for ha in homologDic:
            hb = homologDic[ha]
            if ha not in alleleDic and hb in alleleDic:
                tmp = []
                for i in range(length_allele):
                    tmp.append("NA")
                lineA = "\t".join(tmp)
                lineB = alleleDic[hb]
                fout.write("{}\t{}\n".format(lineA, lineB))
            elif hb not in alleleDic and ha in alleleDic:
                tmp = []
                for i in range(length_allele):
                    tmp.append("NA")
                lineB = "\t".join(tmp)
                lineA = alleleDic[ha]
                fout.write("{}\t{}\n".format(lineA, lineB))
            elif ha not in alleleDic and hb not in alleleDic:
                continue
            else:
                lineA = alleleDic[ha]
                lineB = alleleDic[hb]
                fout.write("{}\t{}\n".format(lineA, lineB))



def pipe(faA, faB, allele, threads, identCutf, covCutf):
    faType = check_faType(faA)
    prefixA = 'A-B'
    run_balst(faA, faB, faType, threads, prefixA)
    prefixB = 'B-A'
    run_balst(faB, faA, faType, threads, prefixB)
    faASize = read_faSize(faA)
    faBSize = read_faSize(faB)
    blastDicA = read_blast('A-B.blast', faASize, identCutf, covCutf)
    blastDicB = read_blast('B-A.blast', faBSize, identCutf, covCutf)
    homologDicc = reciprocal_blast(blastDicA, blastDicB)
    output_allele(homologDicc, allele)

if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("Usage: {}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sys.argv[0], "<cds/pep of A>", "<cds/pep of B>", "<allele_table>", '<threads>', '<identity cutoff>', '<covrage cutoff>'))
        sys.exit()
    faA, faB, allele, threads, identCutf, covCutf = sys.argv[1:7]
    pipe(faA, faB, allele, threads, identCutf, covCutf)
