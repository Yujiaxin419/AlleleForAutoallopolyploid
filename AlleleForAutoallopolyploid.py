#import numpy as np
import pandas as pd
import os
import time
import re
import pysam
import networkx as nx
import argparse
import collections
import itertools as it
#import backbone as bb
import sys


## 整个过程分为几步：
## 对于骨架搭建
### 1. 相似度：所有基因比对到mono祖先基因组，鉴定最优的比对骨架
### 2. 距离：将全部的基因比对到mono祖先基因组上，存在比对区间有覆盖的认为可能是等位基因
### 3. 共线性：和mono基因组存在共线性的可能是等位基因

## 对于contig回补
### 1. 相似度: 所有的基因进行自身比对，鉴定相似的基因聚类
### 2. 距离： 未比对上的基因比对回自身基因组，比对区间和骨架上的基因区间存在覆盖的认为可能是等位基因
### 3. 共线性： 和骨架上的基因组存在共线性的可能是等位基因

## 最后有时间的话还可以生成一个homolog exchange的表

def time_print(info):
	print("\033[32m%s\033[0m %s"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), info))

def make_conf():
    with open('allele_identify.conf', 'w') as fout:
        fout.write("monoPep:\tpepfile\t# pep file of monoploid reference.\n")
        fout.write("selfPep:\tpepfile\t# pep file of this genome\n")
        fout.write("monoCds:\tcdsfile\t# cds file of monoploid reference.\n")
        fout.write("selfCds:\tcdsfile\t# cds file of this genome. \n")
        fout.write("monoGff3:\tmonoGff3\t# gff3 file of monoploid reference\n")
        fout.write("selfGff3:\tmonoGff3\t# gff3 file of this genome. \n")
        fout.write("monoFasta:\tmonoFasta\t# genome dna sequnce of monoploid reference.\n")
        fout.write("selfFasta:\tselfFasta\t# genome dna sequnce of this genome.\n")
        fout.write("threads:\tthreads\t# threads for running diamond, blastp, minimap2.\n")
        fout.write("monoCutoff:\t70,60\t#identity,coverage.\n")
        fout.write("selfCutoff:\t70,60\t#identity,coverage.\n")
        fout.write("monohap:\tnumber\t# number of genome in the merge mono fasta.\n")
        fout.write("selfhap:\tnumber\t# haplotype of this genome.\n")
        fout.write("subgenome:\tsubgenoemeFile\t# file records what chromosomes contanins in each subgenomes.\n")
    print("config file 'allele_identify.conf' generated!")
    sys.exit()

def read_conf(config):
    with open(config,'r') as fin:
        config_dic = {}
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            if "monoPep" in tmpLst:
                config_dic['monoPep'] = tmpLst[1]
            elif "selfPep" in tmpLst:
                config_dic['selfPep'] = tmpLst[1]
            elif "monoCds" in tmpLst:
                config_dic['monoCds'] = tmpLst[1]
            elif "selfCds" in tmpLst:
                config_dic['selfCds'] = tmpLst[1]
            elif "monoGff3" in tmpLst:
                config_dic['monoGff3'] = tmpLst[1]
            elif "selfGff3" in tmpLst:
                config_dic['selfGff3'] = tmpLst[1]
            elif "monoFasta" in tmpLst:
                config_dic['monoFasta'] = tmpLst[1]
            elif "selfFasta" in tmpLst:
                config_dic['selfFasta'] = tmpLst[1]
            elif "threads" in tmpLst:
                config_dic['threads'] = int(tmpLst[1])
            elif "monoCutoff" in tmpLst:
                config_dic['monoCutoff'] = tmpLst[1].split(',')
            elif "selfCutoff" in tmpLst:
                config_dic['selfCutoff'] = tmpLst[1].split(',')
            elif "monohap" in tmpLst:
                config_dic['monohap'] = int(tmpLst[1])
            elif "selfhap" in tmpLst:
                config_dic['selfhap'] = int(tmpLst[1])
            elif "subgenome" in tmpLst:
                config_dic['subgenome'] = tmpLst[1]
            else:
                print("Warning: Cannot recognize term '{}' in config file. please check!".format(line))
                sys.exit()
    return config_dic

def bestHitBlast_with_mono(monoPep, totalPep, thread):
    cmd_buildDB = "makeblastdb -in {} -dbtype prot".format(monoPep)
    cmd_blastp = "blastp -db {} -num_threads {} -outfmt 6 \
                    -out self-mono.blast -query {} -evalue 1e-5 -best_hit_score_edge 0.05 \
                    -best_hit_overhang 0.25  -max_target_seqs 1".format(monoPep, thread, totalPep)
    os.system(cmd_buildDB)
    time.sleep(30)
    os.system(cmd_blastp)
    
def minimap_anno(cds, fasta, hapNum, thread, output):
    mapping_region = {}
    minimap_cmd = "minimap2 -x splice -t {} -k 12 -a -N {} {} {} > {}.sam".format(thread, hapNum, fasta, cds, output)
    print("Minimap CMD: {}".format(minimap_cmd))
    os.system(minimap_cmd)
    samf = pysam.AlignmentFile('{}.sam'.format(output),'r', check_sq=False)
    for item in samf:
        geneName, refChr, refStart, refEnd = item.query_name, item.reference_name, item.reference_start, item.reference_end
        if geneName not in mapping_region:
            mapping_region[geneName] = [(refChr, refStart, refEnd)]
        else:
            mapping_region[geneName].append((refChr, refStart, refEnd))
    return mapping_region
    # return mapping_region ;mapping_region[geneName] = [(refStart, refEnd)]

def read_gff3(gff3):
    gene_db = collections.OrderedDict()
    with open(gff3, 'r') as fin:
        for line in fin:
            if line.strip() == '' or line[0] == '#':
                continue
            data = line.strip().split()
            if data[2] != 'gene':
                continue
            chrn = data[0]
            if "Name" in data[8]:
                id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
            else:
                id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
            gene_db[id] = [(chrn, int(data[3]), int(data[4]))]
    return gene_db
    # return gene_db; gene_db[id] = [(chrn, int(data[3]), int(data[4]))]

def check_overlap(qryName, refName, gene_db, mapping_db):
    qryRegionLst = mapping_db[qryName]
    refRegionLst = gene_db[refName]
    overlap_flag = False
    for qr in qryRegionLst:
        for rr in refRegionLst:
            if qr[0] != rr[0]:
                continue
            if max(qr[1], rr[1]) <= min(qr[2], rr[2]):
                ovlRatio = float(min(qr[2], rr[2])-max(qr[1], rr[1]))/float(min((qr[2]-qr[1]), (rr[2]-rr[1])))
                overlap_flag = True if ovlRatio > 0.5 else False
    return overlap_flag  
    # return overlap_flag True or FALse

def run_daimond(ref, qry, hapnum, thread, outputFile):
    dbName = os.path.splitext(ref)[0]
    diamond_db_CMD = "diamond makedb --in {} --db {}".format(ref, dbName)
    diamond_CMD = "diamond blastp --db {} -q {} -p {} -k {} -e 1e-5 -o {}".format(dbName, qry, thread, hapnum, outputFile)
    time_print("Diamond CMD: '{}'".format(diamond_db_CMD))
    os.system(diamond_db_CMD)
    time_print("Diamond CMD: '{}'".format(diamond_CMD))
    os.system(diamond_CMD)
    time_print("Diamond done!")

def run_synteny(monopep, selfpep, monogff3, hapgff3, hapnum, thread):
    # for fasta
    faDic = {}
    fMono = pysam.FastaFile(monopep)
    for ref in fMono.reference_name:
        faDic[ref] = fMono.fetch(ref)
    fSelf = pysam.FastaFile(selfpep)
    for ref in fSelf.reference_name:
        faDic[ref] = fSelf.fetch(ref)
    os.mkdir('xyz')
    mergePep = 'xyz/merge.pep'
    outputFile = 'xyz/merge.blast'
    with open(mergePep, 'w') as fout:
        for ref in faDic:
            fout.write('>{}\n{}\n'.format(ref, faDic[ref]))
    # for gff
    gene_db = read_gff3(monogff3)
    hapgene_db = read_gff3(hapgff3)
    gene_db.update(hapgene_db)
    mergeGff3 = 'xyz/merge.gff'
    with open(mergeGff3, 'w') as fout2:
        for g in gene_db:
            gchr, gs, ge = gene_db[g][0]
            fout2.write("{}\t{}\t{}\t{}\n".format(gchr, g, gs, ge))
    # diamond
    run_daimond(mergePep, mergePep, hapnum, thread, outputFile)
    # MCScanx
    mcscanx_cmd = "MCScanX xyz/merge"
    time_print("MCScanX CMD: {}".format(mcscanx_cmd))
    os.system(mcscanx_cmd)
    time_print("MCScanx done!")

def read_faLength(fa):
    f = pysam.FastaFile(fa)
    faLenDic = dict(zip(f.references, f.lengths))
    return faLenDic

def read_balst(csv):
    blast_result = pd.read_table(csv, header=None)
    blast_result.columns = ["qry", "ref", "ident", "al", "mismach", 'gap', "rs", "re", "qs", "qe", "evalue", "bhit"]
    blast_result.drop_duplicates(subset=["qry", "ref"])
    return blast_result
    # reruen pd.datafram

def read_synteny(synteny, genedb):
    syntenyDic = {}
    with open(synteny,'r') as fin:
        for line in fin:
            if line.startswith('## Alignment'):
                block_name = line.split(':')[0][13:]
            elif line.startswith("#") ==False:
                geneA, geneB = line.split('\t')[1:3]
                # 保证得到的共线性区块是不同染色体的
                #if genedb[geneA][0][0] == genedb[geneB][0][0]:
                #    continue
                if block_name not in syntenyDic:
                    syntenyDic[block_name] = [[geneA], [geneB]]
                else:
                    syntenyDic[block_name][0].append(geneA)
                    syntenyDic[block_name][1].append(geneB)
    with open('synBlock.txt', 'w') as fin:
        for block in syntenyDic:
            fin.write("{}\t{}\t{}\n".format(block, ','.join(syntenyDic[block][0]), ','.join(syntenyDic[block][1])))
    geneDic = {}
    for gene in genedb:
        chr, pos = genedb[gene][0][:2]
        if chr not in geneDic:
            geneDic[chr] = [(gene,pos)]
        else:
            geneDic[chr].append((gene,pos))
    for chr in geneDic:
        geneDic[chr] = sorted(geneDic[chr], key=lambda x: x[1])
        geneDic[chr] = list(map(lambda x: x[0], geneDic[chr]))
#    geneLst = list(genedb.keys())
    transSynDict = {}
    for block in syntenyDic:
        gaLst, gbLst = syntenyDic[block]
        chra, chrb = genedb[gaLst[0]][0][0], genedb[gbLst[0]][0][0]
        geneLsta, geneLstb = geneDic[chra], geneDic[chrb]
        gas, gae = geneLsta.index(gaLst[0]), geneLsta.index(gaLst[-1])
        gbs, gbe = geneLstb.index(gbLst[0]), geneLstb.index(gbLst[-1])
        gas, gae = (gae, gas) if gas > gae else (gas, gae)
        gbs, gbe = (gbe, gbs) if gbs > gbe else (gbs, gbe)
        gaLst = geneLsta[gas: gae+1]
        gbLst = geneLstb[gbs: gbe+1]
        for gLst in (gaLst, gbLst):
            for g in gLst:
                if g not in transSynDict:
                    transSynDict[g] = [block]
                else:
                    transSynDict[g].append(block)
    return transSynDict
    # return syntenyDic; syntenyDic[gene] = [block1, block2....]


def read_tandem(tandem):  # 这个地方可能需要改一下
    tandemDic = {}
    with open(tandem,'r') as IN:
        for line in IN:
            tmpLst = line.rstrip().split(',')
            tandemDic[tmpLst[1]] = tmpLst[0]
    return tandemDic
    # return tandemDic; tandemDic[tmpLst[1]] = tmpLst[0]

def check_synteny(geneA, geneB, syntenyDic):
    synFlag = False
    try:
        geneA_syn = syntenyDic[geneA]
    except:
       # print("can't find synteny {}".format(geneA))
        return synFlag
    try:
        geneB_syn = syntenyDic[geneB]
    except:
       # print("can't find synteny {}".format(geneB))
        return synFlag
    if len(set(geneA_syn).union(set(geneB_syn))) < len(geneA_syn) + len(geneB_syn):
        synFlag = True
    return synFlag
    # return True or False

def build_simiGene_network(filtered_svdf):
    homologDic = {}
    #new_svdf = pd.DataFrame(columns=["qry", "ref", "ident", "col"])
    ## 构建相似性网络
    G = nx.from_pandas_edgelist(filtered_svdf,'ref', 'qry')
    ## 遍历子图
    for c in nx.connected_components(G):
        subG = G.subgraph(c)
        nodeSet = list(subG.nodes)
        ### 选择图中度最多的点作为reference重构相似表
        degreLst = [(node, int(subG.degree(node))) for node in nodeSet ]
        degreLst.sort(reverse=True, key=lambda x: x[1])
        ref = degreLst[0][0]
        ### 因为构图的时候是采用贪婪构图，只要任意一条有相似性就会被拉进来，还要考虑一下这样合不合理？
        for node in nodeSet:
            homologDic[node] = ref
    return homologDic

def read_info_db(db):
    infoDB = {}
    with open(db, 'r') as fin:
        fin.readline()
        for line in fin:
            tmpLst = line.strip().split('\t')
            ref, qry, ident, cov, overlap, syn = tmpLst
            infoDB[(ref, qry)] = {'ident': float(ident), 'cov':float(cov), 'overlap': overlap==str(True), 'synteny': syn==str(True) }
    return infoDB


def check_homo_info(ref, qry, infoDB, info, defaultValue):
    if (ref, qry) in infoDB:
        value = infoDB[(ref, qry)][info]
    elif (qry, ref) in infoDB:
        value = infoDB[(qry, ref)][info]
    else:
        value = defaultValue
    return value


def back_bone(config_dic):
    if not os.path.exists("self-mono.blast"):
        time_print("Running blastp")
        bestHitBlast_with_mono(config_dic['monoPep'], config_dic['selfPep'], config_dic['threads'])
    else:
        time_print("Find bestHit result!")
    minimap_output_prefix = 'self_self'
    if not os.path.exists("{}.sam".format(minimap_output_prefix)):
        time_print("Running minimap2")
        mapping_region = minimap_anno(config_dic['selfCds'], config_dic['selfFasta'], config_dic['selfhap'], config_dic['threads'], minimap_output_prefix)
    else:
        time_print("Find mapping result!")
        mapping_region = {}
        samf = pysam.AlignmentFile('{}.sam'.format(minimap_output_prefix),'r')
        for item in samf:
            geneName, refChr, refStart, refEnd = item.query_name, item.reference_name, item.reference_start, item.reference_end
            if geneName not in mapping_region:
                mapping_region[geneName] = [(refChr, refStart, refEnd)]
            else:
                mapping_region[geneName].append((refChr, refStart, refEnd))
    if not os.path.exists('xyz/merge.collinearity'):
        time_print("Running Diamond & MCScanx")
        run_synteny(config_dic['monoPep'], config_dic['selfPep'], config_dic['monoGff3'], config_dic['selfGff3'], config_dic['selfhap'] + config_dic['monohap'], config_dic['threads'])
    else:
        time_print("Find synteny file!")
    if not os.path.exists('self_self.blast'):
        time_print("Running Diamond")
        run_daimond(config_dic['selfPep'], config_dic['selfPep'], config_dic['selfhap'], config_dic['threads'], 'self_self.blast')
    homoDic = {}
    with open("self-mono.blast", 'r') as fin:
        for line in fin:
            tmpLst = line.rstrip().split('\t')
            qry, ref = tmpLst[:2]
            if ref not in homoDic:
                homoDic[ref] = []
            homoDic[ref].append(qry)
    #        homoDic[qry].append(ref)
    tmpDic = {'ref':[], 'qry':[]}
    for ref in homoDic:
        if len(homoDic[ref]) == 1:
            tmpDic['ref'].append(homoDic[ref][0])
            tmpDic['qry'].append(homoDic[ref][0])
        else:
            for e in it.combinations(homoDic[ref], 2):
                tmpDic['ref'].append(e[0])
                tmpDic['qry'].append(e[1])
    homo_df = pd.DataFrame.from_dict(tmpDic)
    homo_df.to_csv('homo_df.txt', header=True, sep='\t')
    #sys.exit()
    # 读取blast file
    faLenDB = read_faLength('xyz/merge.pep')
    self_blast_df = read_balst("self_self.blast")
    monoGeneDB = read_gff3(config_dic['monoGff3'])
    #filter_result['overlap'] = filter_result.apply(lambda x: check_overlap(x['qry'], x['ref'], monoGeneDB, mapping_region), axis=1)
    hapGeneDB = read_gff3(config_dic['selfGff3'])
    totalGeneDB = monoGeneDB
    totalGeneDB.update(hapGeneDB)
    syntenyDic = read_synteny('xyz/merge.collinearity', totalGeneDB)
    ## 添加覆盖度
    self_blast_df['cov'] = self_blast_df.apply(lambda x: 100*float(x['al'])/float(min(faLenDB[x['ref']], faLenDB[x['qry']])), axis=1)
    ## 添加坐标覆盖信息
    self_blast_df['overlap'] = self_blast_df.apply(lambda x: check_overlap(x['qry'], x['ref'], hapGeneDB, mapping_region), axis=1)
    ## 添加共线性
    self_blast_df['synteny'] = self_blast_df.apply(lambda x: check_synteny(x['qry'], x['ref'], syntenyDic), axis=1)
    self_blast_df = self_blast_df[['ref', 'qry', 'ident', 'cov', 'overlap', 'synteny']]
    #self_blast_df.to_csv('blast_all_info_database.txt',  header=True, sep="\t")
    ## 建立多重索引
    self_blast_df = self_blast_df.set_index(["ref", "qry"])
    self_blast_df = self_blast_df.sort_index()
    self_blast_df.to_csv('blast_all_info_database.txt',  header=True, sep="\t")
    time_print("read info db")
    infoDB = read_info_db('blast_all_info_database.txt')
    ## 添加到同源基因表里
    time_print("check ident")
    homo_df['ident'] = homo_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'ident', 0), axis=1)
    time_print("check cov")
    homo_df['cov'] = homo_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'cov', 0), axis=1)
    time_print("check overlap")
    homo_df['overlap'] = homo_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'overlap', False), axis=1)
    time_print("check synteny")
    homo_df['synteny'] = homo_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'synteny', False), axis=1)
    ## 添加亚基因组同源基因信息
    qry2ref = {}
    for ref in homoDic:
        for qry in homoDic[ref]:
            qry2ref[qry] = ref
    homo_df['homolog'] = homo_df.apply(lambda x: qry2ref[x['ref']], axis=1)
    ## 筛选同源基因表
    identityCutoff = float(config_dic['monoCutoff'][0])
    coverageCutoff = float(config_dic['monoCutoff'][1])
    filter_result = homo_df[ (homo_df['ident']> identityCutoff) & (homo_df['cov']>coverageCutoff)]
    backbone_result = filter_result[filter_result['overlap'] & filter_result['synteny']]
    homo_df.to_csv('unFilterHomolog.txt',  header=True, sep="\t")
    backbone_result.to_csv('backbone_result.txt',  header=True, sep="\t")
    #sys.exit()
    return backbone_result, syntenyDic, infoDB
    # return backbone_result


def rescue_gene(config_dic, backbone_result, infoDB):
    f = pysam.FastaFile(config_dic['selfCds'])
    totalGene = list(f.references)
    anchorGene = list(backbone_result['qry'])
    anchorGene.extend(list(backbone_result['ref']))
    ag = {}
    for g in anchorGene:
        if g not in ag:
            ag[g] = ''
    #anchorGene.extend(list(backbone_result['qry']))
    ## 输出未挂载上的基因
    time_print("output unanchor genes!")
    unAnchor = {}
    for g in totalGene:
        if g not in ag:
            unAnchor[g] = ''
    tmpDic = {'ref':[], 'qry':[]}
    tmpDic2 = {}
    time_print("check unAnchor-pairs")
    for blast_pair in infoDB.keys():
        if frozenset(blast_pair) in tmpDic2:
            continue
        else:
            tmpDic2[frozenset(blast_pair)] = ""
        ref, qry = blast_pair
        if ref in unAnchor and qry in unAnchor:
            tmpDic['ref'].append(ref)
            tmpDic['qry'].append(qry)
    unAnchor_df = pd.DataFrame.from_dict(tmpDic)
    unAnchor_df.to_csv('unAnchor_df.txt', header=True, sep='\t')
    #### add info   
    time_print("check ident")
    unAnchor_df['ident'] = unAnchor_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'ident', 0), axis=1)
    time_print("check cov")
    unAnchor_df['cov'] = unAnchor_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'cov', 0), axis=1)
    time_print("check overlap")
    unAnchor_df['overlap'] = unAnchor_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'overlap', False), axis=1)
    time_print("check synteny")
    unAnchor_df['synteny'] = unAnchor_df.apply(lambda x: check_homo_info(x['ref'], x['qry'], infoDB, 'synteny', False), axis=1)
    identityCutoff = float(config_dic['selfCutoff'][0])
    coverageCutoff = float(config_dic['selfCutoff'][1])
    ## 筛选高相似度的基因对
    filter_result = unAnchor_df[ (unAnchor_df['ident']> identityCutoff) & (unAnchor_df['cov']>coverageCutoff)]
    filter_result = filter_result[['ref', 'qry', 'ident', 'cov', 'overlap', 'synteny']]
    #filter_result.to_csv("unanchor-self.df", header=True, sep="\t")
    time_print('get ref gene of each cluster uniq gene.')
    homologDic = build_simiGene_network(filter_result)
    ## add info
    time_print('check homolog')
    filter_result['homolog'] = filter_result.apply(lambda x: homologDic[x['ref']], axis=1)
    filter_result.to_csv("unanchor.df", header=True, sep="\t")
    filter_result = filter_result[filter_result['overlap'] & filter_result['synteny']]
    ## 合并三个表
    finalAlleleTable = pd.concat([backbone_result, filter_result])
    return finalAlleleTable


def read_subgenome(sgFile):
    mono2sg2chrDic = {}
    chr2sgDic = {}
    #print(sgFile)
    with open(sgFile, 'r') as IN:
        for line in IN:
            tmpLst = line.split('\t')
            sg, chrLst = tmpLst[0], tmpLst[1].split(',')
     #       print(sg, chrLst)
            test_chr = chrLst[0][-1]
            for chrn in chrLst:
                chr2sgDic[chrn] = sg
    return chr2sgDic
    # return mono2sg2chrDic[monochr][sg] = chrLst

def rebuild_allele_table(finalAlleleTable, config_dic):
    ##  gene chr pos SS-A SS-B Rec
    chr2sgDic = read_subgenome(config_dic['subgenome'])
    monogeneDB = read_gff3(config_dic['monoGff3'])
    hapGeneDB = read_gff3(config_dic['selfGff3'])
    tandemDic = read_tandem('xyz/merge.tandem')
    sgdb = []
    for chr in chr2sgDic:
        sg = chr2sgDic[chr]
        if sg not in sgdb:
            sgdb.append(sg)
    outdb = []
    for i in range(len(sgdb)):
        outdb.append([])
    #outdic = dict(zip(sgdb,outdb))
    #monochrLst = list(chr2sgDic.keys())
    flipAlleleDict = {}
    for row in finalAlleleTable.itertuples():
        refG, qryG1, qryG2 = getattr(row, 'homolog'), getattr(row, 'ref'), getattr(row, 'qry')
        for qryG in [qryG1, qryG2]:
            qryC = hapGeneDB[qryG][0][0]
            try:
                qrysg = chr2sgDic[qryC]
            except:
                continue
            if refG not in flipAlleleDict:
                flipAlleleDict[refG] = {}
                for sg in sgdb:
                    flipAlleleDict[refG][sg] = {} 
            if qryC not in flipAlleleDict[refG][qrysg]:
                flipAlleleDict[refG][qrysg][qryC] = []
            if qryG in flipAlleleDict[refG][qrysg][qryC]:
                continue
            flipAlleleDict[refG][qrysg][qryC].append(qryG)
# Chr0719467      LAmono_Chr01    3344792 NA      NA      YZ-So-Chr07C
    for refG in flipAlleleDict:
        for sg in sgdb:
            if bool(flipAlleleDict[refG][sg]) == False:
                continue
            for chrn in flipAlleleDict[refG][sg]:
                qglst = flipAlleleDict[refG][sg][chrn]
                if len(qglst)==1:
                    continue
                else:
                    outqLst = []
                    for qg in qglst:
                        if qg in tandemDic:
                            qgT = qg + "-T"
                            outqLst.append(qgT)
                            qglst.remove(qg)
                    alleleQG = qglst[0]
                    outqLst.insert(0,alleleQG)
                    qglst.remove(alleleQG)
                    if len(qglst)>0:
                        qglst = list(map(lambda qg: qg + '-P', qglst))
                        outqLst.extend(qglst)
                    flipAlleleDict[refG][sg][chrn] = ['|'.join(outqLst)]
                    #print('|'.join(outqLst))

    ## 开始输出
    ### 给等位基因排个序
    allgeneDB = monogeneDB
    allgeneDB.update(hapGeneDB)
    tmpgeneDic = {}
    for refG in flipAlleleDict.keys():
        refC, refPos = allgeneDB[refG][0][:2]
        if refC not in tmpgeneDic:
            tmpgeneDic[refC] = [(refG, refPos)]
        else:
            tmpgeneDic[refC].append((refG, refPos))
    for refC in tmpgeneDic:
        tmpgeneDic[refC] = sorted(tmpgeneDic[refC], key=lambda x: x[1])
    ### 开始输出了
    with open("Allele.finnal.table",'w') as fout:
        ## 先确定亚基因组顺序
        header = "#REF\tCHR\tPOS\t{}\n".format('\t'.join(sgdb))
        fout.write(header)
        for refC in tmpgeneDic:
            for refG in tmpgeneDic[refC]:
                refG = refG[0]
                geneStr = []
                for sg in sgdb:
                    if bool(flipAlleleDict[refG][sg]):
                        try:
                            sgGeneStr = ",".join(list(flipAlleleDict[refG][sg].values()))
                        except:
                            sgGeneStr = ",".join(list(map(lambda x: ','.join(x), flipAlleleDict[refG][sg].values())))                      
                    else:
                        sgGeneStr = 'NA'
                    geneStr.append(sgGeneStr)
                fout.write("{}\t{}\t{}\t{}\n".format(refG, allgeneDB[refG][0][0], allgeneDB[refG][0][1], "\t".join(geneStr)))
    ## 如果按照祖先基因组输出的话，可能会存在有没有可能存在亚基因组间的同源交换，导致有些基因的等位基因弥散到另一条染色体上去了
    ## 需要区分挂载上的基因的属于哪个单倍型
    ## 同一个单倍型中如果存在多个基因，第一步先去掉tandem repeat
    ## 如果还有其他的单倍型，就比较和mono的位置，最近的是allele，其余的是paralog



def pipe(config, make_config_flag):
    if make_config_flag:
        make_conf()
    else:
        pass
    time_print("Run with config file '{}'.".format(config))
    time_print("Read config file.")
    config_dic = read_conf(config)
    time_print("Start to build allele table backbone.")
    backbone_result, syntenyDic, infoDB = back_bone(config_dic)
    with open('new_syntenyDic.txt', 'w') as fout:
        for i in syntenyDic:
            fout.write("{}\t{}\n".format(i,",".join(syntenyDic[i])))
   # sys.exit()
    #backbone_result = pd.read_csv("backbone_result.txt", sep="\t")
    with open('new_syntenyDic.txt', 'r') as fin:
        syntenyDic = {}
        for line in fin:
            tmppLst = line.rstrip().split('\t')
            gene, blockN = tmppLst
            blockN = blockN.split(',')
            syntenyDic[gene] = blockN
    time_print("Start to rescue unanchor genes.")
    finnalAllele = rescue_gene(config_dic, backbone_result, infoDB)
    time_print("Start to rebuild allele table and output.")
    finnalAllele.to_csv('finnalAllele.txt',  header=True, sep="\t")
    rebuild_allele_table(finnalAllele, config_dic)
    time_print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is the script for identifying allelic gene in an autoallopolyploid genome")
    parser.add_argument('-i', '--input_config', default=None,
                        help='<filepath>  The config file.')
    parser.add_argument('-c', '--make_config', action="store_true",
                        help='generate a un-edited config file.')
    args = parser.parse_args()
    pipe(args.input_config, args.make_config)





    
    
    


                

            
