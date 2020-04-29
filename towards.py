import numpy as np
import sys
import random


def target(a, b):
    return 2*a*a + b - 57

def guess_initial_set(num_guesses, num_numbers):
    return np.random.randint(0, high=15, size=(num_guesses, num_numbers))

def fitness(data):
    res = []
    for l in data:
        grade = target(l[0], l[1])
        if grade == 0:
            print ("Found optimal solution a={} b={}".format(l[0], l[1]))
            sys.exit(0)
        res.append(1.0 / np.abs(grade))
    return res

def natural_selection(FP):
    res = []
    s = 0
    FPT = [0.0]
    for p in FP:
        FPT.append(p + s)
        s += p
    FPT.append(1.0)
    Pr = sorted(np.random.uniform(0.0, 1.0, len(FP)))
    Pr = Pr
    iFP = 0
    while len(Pr):
        if (Pr[0] >= FPT[iFP] and Pr[0] < FPT[iFP + 1]):
            yield(iFP)
            Pr = Pr[1:]
        else:
            iFP += 1
    return

def split(word):
    return [char for char in word]

def get_bin_chromosome(p):
    res = []
    for g in bin(p[0])[2:].zfill(4) + bin(p[1])[2:].zfill(4):
        res.append(g == '1')
    return res

def to01string(s):
    res = []
    for b in s:
        res.append('1' if b else '0')
    return res

def get_int_pair_chromosome(s):
    c1 = to01string(s[:4])
    c2 = to01string(s[4:])
    return np.array((int(''.join(c1),2), int(''.join(c2), 2)))

def breed(parent1, parent2):
    child = []
    childP1 = []
    childP2 = []

    geneA = int(random.random() * len(parent1))
    geneB = int(random.random() * len(parent1))

    startGene = min(geneA, geneB)
    endGene = max(geneA, geneB)

    return parent1[:startGene] + parent2[startGene:endGene] + parent1[endGene:]

def breedPopulation(matingpool, eliteSize):
    children = []
    length = len(matingpool) - eliteSize
    pool = random.sample(matingpool, len(matingpool))

    for i in range(0,eliteSize):
        children.append(matingpool[i])

    for i in range(0, length):
        child = breed(pool[i], pool[len(matingpool)-i-1])
        children.append(child)
    return children

def sort_parents(parents):
    return sorted(parents, key=lambda x:x[1], reverse=True)

def mutate(children, rate):
    l = len(children[0])
    n_children = len(children)
    tot_num_genes = n_children * l
    num_mutations = int(np.rint(rate * tot_num_genes))
    num_indxs = np.random.randint(0, high=tot_num_genes - 1, size=num_mutations)
    for i in num_indxs:
        #print (i)
        i_child = int(i / l)
        i_gene = i % l
        #print ("mutating {},{}".format(i_child, i_gene))
        children[i_child][i_gene] = not children[i_child][i_gene]

def main():
    crms = guess_initial_set(6, 2)
    for i in range(1024):
        print (i)
        fit = fitness(crms)
        tot = sum(fit)
        FP = [ g / tot for g in fit]
        #print (FP)
        ns = list(natural_selection(FP))
        parents = [(crms[i], fit[i]) for i in ns]
        #print (parents)
        parents = sort_parents(parents)
        par_chromo = [ get_bin_chromosome(p[0]) for p in parents]
        eliteSize = 2
        children = breedPopulation(par_chromo, eliteSize)
        mutate(children, 0.1)
        #print (children)
        crms = [get_int_pair_chromosome(c) for c in children]
        #print (crms)
        #print (crms)
    print (crms)
    print ([target(c[0], c[1]) for c in crms])


main()
