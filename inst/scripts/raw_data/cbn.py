#!/usr/bin/env python

"""
  cbn.py -- conjunctive Bayesian networks
  Algorithms for data handling.

  Version : 0.1.00
  Author  : Niko Beerenwinkel and Moritz Gerstung
    
  Copyright (C)  2011 
  Niko Beerenwinkel                 Moritz Gerstung
  niko.beerenwinkel@bsse.ethz.ch    moritz.gerstung@bsse.ethz.ch
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
"""


import sys
import numpy as np

def read_tumors(name):
    '''Reads the tumor data.
    Input: File name
    Output: 2 dictionaries Tumor->Genes, Genes->Tumors
    '''
    file = open(name)
    header = file.next().split("\t")
    tu = header.index('Tumor')
    ge = header.index('Gene')
    try:
        ty = header.index('Mutation Type')
    except:
        ty = 0
    data = {}
    genes = {}
    for line in file:
        ls = line.split("\t")
        if ls[ty] != 'Synonymous': 
            tumor, gene = ls[tu].upper(), ls[ge]
            if gene != 'NONE':
                if tumor in data.keys():
                    data[tumor].append(gene)
                else:
                    data[tumor] = [gene]
                if gene in genes.keys():
                    genes[gene].append(tumor)
                else:
                    genes[gene] = [tumor]
            else:
                data[tumor] = []
    file.close()
    return data, genes

def read_coregroups(name = '../rawData/CoreGroup2Gene.txt'):
    '''Returns two dictionaries: CORE2gene and gene2CORE.
    '''
    file = open(name)
    line = file.next()
    CORE2gene = {}
    gene2CORE = {}
    for line in file:
        core, gene = line.strip().split("\t")
        if core == 'Invasion & Apoptosis': #Append to both
            core = 'Invasion'
            if core in CORE2gene.keys():
                CORE2gene[core].append(gene)
            else:
                CORE2gene[core] = [gene]
            core = 'Apoptosis'
        if core in CORE2gene.keys():
            CORE2gene[core].append(gene)
        else:
            CORE2gene[core] = [gene]
        if gene in gene2CORE.keys():
            gene2CORE[gene].append(core)
        else:
            gene2CORE[gene] = [core]
    file.close()
    return CORE2gene, gene2CORE

def tumor2core(tumors, gene2CORE):
    '''Input: dictionary 'tumors' with tumor to gene mapping and 'gene2CORE' dictionary with gene to CORE group mapping.
    Output: dictionary with tumor to CORE group mapping.
    '''
    tumor2CORE = {}
    for t in tumors:
        tumor2CORE[t] = []
        for g in tumors[t]:
            if g in gene2CORE.keys():
                for c in gene2CORE[g]:
                    if c not in tumor2CORE[t]:
                        tumor2CORE[t].append(c)
    return tumor2CORE

def randomize_tumors(tumors):
    """Permutes mutations across tumors.
    Input: Tumor dictionary
    Output: Permuted dictionary 
    """
    rumors = {}
    keys = tumors.keys()
    k = len(keys)
    for key in tumors:
        rumors[key] = []
    for key in tumors:
        l = len(tumors[key])
        r = np.random.randint(0,k,l)
        for i in range(l):
            rumors[keys[r[i]]].append(tumors[key][i])
    return rumors
            
def tumors2pat(tumors, genes = None, sort=False):
    '''Input: tumor (dictionary)
    Output: Binary pattern of mutations. Rows: tumors, cols: genes (alphabetical)
    '''
    if not genes:
        genes = set()
        for t in tumors.values():
            for g in t:
                genes.add(g)
        genes = list(genes)
    genes.sort()
    pat = np.zeros([len(tumors), len(genes)], dtype = 'int')
    tum = tumors.keys()
    if sort:
        tum.sort()
    for i in range(len(tum)):
        for j in range(len(genes)):
            pat[i,j] = int(genes[j] in tumors[tum[i]])
    return pat

def write_pat(pat, file = None):
    '''Prints pattern pat to file, otherwise to stdout
    '''
    if file:
        saveout = sys.stdout
        fh = open(file, 'w')
        sys.stdout = fh
    m, n = pat.shape[0], pat.shape[1]
    print m, n+1
    for i in range(m):
        print 1,
        for j in range(n):
            print int(pat[i,j]),
        print "\n",
    if file:
        fh.close()
        sys.stdout = saveout

def read_pat(file):
    '''Read pattern from file.
    '''
    fileh = open(file, 'r')
    fileh.next()
    pat = []
    for line in fileh:
        pat.append(map(int, line.split()[1:]))
    pat = np.matrix(pat)
    fileh.close()
    return pat

def read_poset(filen):
    '''Read poset from file.
    '''
    file = open(filen, 'r')
    m  = int(file.next())
    poset = np.zeros([m,m])
    for line in file:
        if line.strip() == '0':
            break
        i, j = map(int, line.split())
        poset[i-1,j-1] = 1
    file.close()
    return poset

def write_poset(poset, filen):
    '''Write Poset to file.
    '''
    file = open(filen, 'w')
    m = poset.shape[1]
    print >>file, m
    for i in range(m):
        for j in range(m):
            if poset[i,j]:
                print >>file, i+1, j+1
    print >>file, 0
    file.close()

def pos2rel(poset):
    '''Return list of relations.
    '''
    posrel = []
    i = 1
    for line in poset:
        j = 1
        for col in line:
            if col == 1:
                posrel.append([i, j])
            j += 1
        i+=1
    return posrel

def rel2pos(rel, m):
    '''Return poset as matrix.
    '''
    poset = np.zeros([m,m])
    for r in rel:
        poset[r[0] -1, r[1]-1] += 1
    return poset

def linear_poset(pat):
    '''Input: pattern
    Output: linear poset
    '''
    n, m = pat.shape
    sorted = np.array(pat).mean(0).argsort()[::-1]
    #print sorted
    poset = np.zeros([m,m])
    s = sorted[0]
    for t in sorted[1:]:
        poset[s,t] = 1
        s = t
    return poset

def write_label(label, file):
    """Writes LABEL to FILE.
    """
    file = open(file, "w")
    for l in label:
        print >>file, l
    file.close()

if __name__=='__main__':

    data, genes = read_tumors()
    #print data

    fav_genes = ['TP53', 'APC', 'KRAS', 'EVC2', 'FBXW7', 'PIK3CA', 'EPHA3', 'TCF7L2']

    gene_count = {}
    for gene in genes:
        gene_count[gene] = 0
        for tumor in data:
            gene_count[gene] = gene_count[gene] + int(gene in data[tumor])

    #print gene_count

    print len(data), len(fav_genes)+1
    for tumor in data.keys():
        print 1, 
        for gene in fav_genes:
            print int(gene in data[tumor]), 
        print "\n",

def poset_dist(poset,posetHat):
    #fpos = inPoset(posetHat,poset)
    #fneg = inPoset(poset,posetHat)
    #return fpos, fneg
    allrel,allrelHat = BFS(poset), BFS(posetHat)
    fpos, fneg = len(allrelHat), len(allrel)
    tpos = 0
    for rel in allrel:
        if rel in allrelHat:
            fpos -= 1
            tpos += 1

    for rel in allrelHat:
        if rel in allrel:
            fneg -= 1
    
    l = float(max(1,len(allrel)))
    return [fpos,tpos,fneg,l]

def BFS(relations):
    '''Breadth-first search of poset. Returns list with all relations.
    '''
    allrel = list()
    queue = list()
    for rel1 in relations:
        queue.append(rel1[0])
        while len(queue) > 0:
            mut = queue.pop(0)
            for rel2 in relations:
                if rel2[0] == mut:
                    queue.append(rel2[1])
                    if [rel1[0],rel2[1]] not in allrel:
                        allrel.append([rel1[0],rel2[1]])
    return allrel

def read_primary_type(name):
    """Read annotation and if primary tumor
    """
    file = open(name)
    g = {}
    file.next()
    for line in file:
        s = line.strip().split("\t")
        try:
            if s[11] == "Prevalence":
                g[s[0].upper()] = int(s[5]=="Yes")
        except:
            print s
    return g

def lambda_to_fitness(name, T0 = 3650, N = 1e6, mu=1e-7):
    """Input: 'experiment.lambda'
    Output: np.array of fitness values
    """
    file = open(name)
    fitness = []
    for line in file:
        l = float(line.strip())
        f = (1 + np.log(mu))**2 / (2 * T0/l + np.log(N) + 1 + np.log(mu)) 
        fitness.append(f)
    return np.array(fitness)

def bootstrap(pat):
    '''Input: pattern
    Output: bootstrapped pattern
    '''
    n, m = pat.shape
    bat = pat[np.random.randint(n, size=n),:]
    return bat

def randomize(pat):
    '''Input: pattern
    Output: permuted pattern (columnwise)
    '''
    n, m = pat.shape
    pet = np.array(pat)
    for j in range(m):
        pet[:,j] = pat[np.random.permutation(n),j]
    return pet
