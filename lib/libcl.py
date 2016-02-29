#!/usr/bin/env python

import copy

def initialize(pool):
    cl = [[pool[i]] for i in range(len(pool))]
    return cl

def dist(m,n,cl,tm):
    d = 0.0
    for x in cl[m]:
        for y in cl[n]:
            d += tm[(x,y)]
    return d/len(cl[m])/len(cl[n])

def spread(m,cl,tm):
    sp = 0.0
    for i in cl[m]:
        for j in cl[m]:
            if i <= j:
                continue
            sp += tm[(i,j)]
    n = len(cl[m])
    if n == 1:
        return 0.0
    else:
        return sp/(n*(n-1)/2)

def AvSp(cl,tm):
    sum = 0.0
    for m in range(len(cl)):
        sum += spread(m,cl,tm)
    return sum/float(len(cl))

def penalty(AvSp_s):
    min_AvSp = min(AvSp_s)
    max_AvSp = max(AvSp_s)
    nn = len(AvSp_s)-1

    norm = []
    for i,AvSp_i in enumerate(AvSp_s):
        norm.append( (i, -i-1,nn/(max_AvSp-min_AvSp)*(AvSp_i-min_AvSp)+1 + (nn-i+1)))
    return norm

def find_nearest(cl,tm):
    d = []
    n = len(cl)
    for i in range(n-1):
        for j in range(i+1,n):
            d.append((i,j,dist(i,j,cl,tm)))
    d.sort(key=lambda x:x[2])
    return d[0]

def merge_cl(cl,tm):
    new = []
    n = len(cl)
    nearest = find_nearest(cl,tm)
    merged = cl[nearest[0]]
    merged.extend(cl[nearest[1]])
    new.append(merged)
    for i in range(n):
        if i == nearest[0] or i == nearest[1]:
            continue
        new.append(cl[i])
    return new

def find_cntr(cl,tm):
    new = []
    for X in cl:
        n = len(X)
        d = [[i,0.0] for i in range(n)]
        for i in range(n-1):
            for j in range(i+1,n):
                d[i][1] += tm[(X[i],X[j])]
                d[j][1] += tm[(X[i],X[j])]
        d.sort(key=lambda x:x[1])
        #
        C = []
        for x in d:
            C.append(X[x[0]])
        new.append(C)
    return new

def NMRclust(pool,tm):
    N = len(pool)
    cl = initialize(pool)
    cl_bank = []
    AvSp_s = []
    for i in range(N-1):
        cl = merge_cl(cl,tm)
        c_cl = copy.deepcopy(cl)
        cl_bank.append(c_cl)
        AvSp_s.append(AvSp(cl,tm))
    P = penalty(AvSp_s)
    P.sort(key=lambda x:x[2])
    P0 = P[0]
    #
    n_cl = N+P0[1]
    i_cl = P0[0]

    cl = find_cntr(cl_bank[i_cl],tm)
    cl.sort(key=lambda x:len(x),reverse=True)
    return n_cl,cl
