#!/usr/bin/env python
import os, sys, copy, math
from subprocess import *

def calc_dist2(a,b,c, d,e,f):
    return (a-d)**2 + (b-e)**2 + (c-f)**2 

def file2dict(pdb_file, n_res_rec):
    PDBrec = {}
    PDBlig = {}
    for line in file(pdb_file):
        if not line.startswith('ATOM'):
            continue
        
        i_res = int(line[22:26])
        atom_name = line[12:16]
        if i_res <= n_res_rec:
            if not i_res in PDBrec:
                PDBrec[i_res] = {}
            PDBrec[i_res][atom_name] = line[:-1]
        else:
            if not i_res in PDBlig:
                PDBlig[i_res] = {}
            PDBlig[i_res][atom_name] = line[:-1]

    return PDBrec, PDBlig

def fn_rmsd(PDB1, PDB2):
    res_s = PDB1.keys()
    res_s.sort()
    dist_sum2 = 0.0
    res_num = 0
    for i_res in res_s:
        try:
            x1 = float(PDB1[i_res][' CA '][30:38])
            y1 = float(PDB1[i_res][' CA '][38:46])
            z1 = float(PDB1[i_res][' CA '][46:54])
            x2 = float(PDB2[i_res][' CA '][30:38])
            y2 = float(PDB2[i_res][' CA '][38:46])
            z2 = float(PDB2[i_res][' CA '][46:54])
            dist2 = calc_dist2(x1,y1,z1, x2,y2,z2)
        except:
            dist2 = 99999.9
        dist_sum2 += dist2
        res_num += 1

    dist_sum2 = dist_sum2 / float(res_num)
    rmsd = math.sqrt(dist_sum2)
    return rmsd

def calc_pair_rmsd(curr_list, n_res_rec):
    pdb_s = []
    for line in curr_list:
        if line.startswith('#'):
            continue
        linesp = line.strip().split()
        pdb_file = linesp[0]
        PDBrec, PDBlig = file2dict(pdb_file, n_res_rec)
        pdb_s.append(PDBlig)

    len_p = len(pdb_s)
    rmsd_list = []
    #initialize_rmsd_list
    for i in range(len_p):
        sub_list = [0.0 for j in range(len_p)]
        rmsd_list.append(sub_list)

    #calc_average_rmsd_of_rmsd_list
    for i in range(len_p-1):
        for j in range((i+1),len_p):
            PDBlig1 = pdb_s[i]
            PDBlig2 = pdb_s[j]
            rmsd = fn_rmsd(PDBlig1,PDBlig2)
            rmsd_list[i][j] = rmsd
            rmsd_list[j][i] = rmsd

    return rmsd_list

min_len_clus = 2
max_len_clus = 50

class Data:
    def __init__(self):
        rmsd = None
        score = None
        pdb_name = None

def len_cmp(a,b):
    return cmp(len(b), len(a))

def parse_list(input_list):
    data_s = []
    idx = 0
    # for line in file(list_file):
    for line in input_list:
        if line.startswith('#'):
            continue
        linesp = line.strip().split()
        pdb_name = linesp[0]
        data = Data()
        data.pdb_name = os.path.abspath(pdb_name)
        idx += 1
        data.score = idx
        data_s.append(data)
    return data_s

def find_min_dist(cluster, d_mat):
    n_clus = len(cluster)
    min_val = 9999999.9

    for i_clus in range(n_clus-1):
        for j_clus in range(i_clus+1, n_clus):
            if (len(cluster[i_clus])+len(cluster[j_clus])) > max_len_clus:
                dist = d_mat[i_clus][j_clus] + 9999.9
            else:
                dist = d_mat[i_clus][j_clus]
            dist = d_mat[i_clus][j_clus]
            if dist < min_val:
                min_val = d_mat[i_clus][j_clus]
                iclus_min = i_clus
                jclus_min = j_clus
    
    return iclus_min, jclus_min

def merge_cluster(cluster, i_clus, j_clus):
    new_cluster = []
    list = []
    for index in cluster[i_clus]:
        list.append(index)
    for index in cluster[j_clus]:
        list.append(index)
    new_cluster.append(list)
    for k_clus in range(len(cluster)):
        if (k_clus == i_clus) or (k_clus == j_clus):
            continue
        new_cluster.append(cluster[k_clus])

    return new_cluster

def dist_cluster(clus1, clus2, pair_rmsds):
    len1 = len(clus1)
    len2 = len(clus2)
    D = 0.0
    for i in range(len1):
        i_bank = clus1[i]
        for j in range(len2):
            j_bank = clus2[j]
            D += pair_rmsds[i_bank][j_bank]
    D = D / float(len1*len2)
    return D

def update_d_mat(cluster, pair_rmsds):
    n_clus = len(cluster)
    d_mat = []
    for i_clus in range(n_clus):
        list = [0.0 for j_clus in range(n_clus)]
        d_mat.append(list)
    
    for i_clus in range(n_clus-1):
        for j_clus in range(i_clus+1, n_clus):
            D = dist_cluster(cluster[i_clus], cluster[j_clus], pair_rmsds)
            d_mat[i_clus][j_clus] = D
            d_mat[j_clus][i_clus] = D
    return d_mat

def fn_spread(pair_rmsds, clus):
    n_bank = len(clus)
    dist_sum = 0.0
    for i in range(n_bank-1):
        for j in range(i+1, n_bank):
            i_bank = clus[i]
            j_bank = clus[j]
            dist_sum += pair_rmsds[i_bank][j_bank]
    
    n_bank = float(n_bank)
    spread = 2.0*dist_sum / (n_bank*(n_bank-1.0))

    return spread

def avg_spread(pair_rmsds, cluster):
    cnum = 0
    avg_spread = 0.0
    for clus in cluster:
        if (len(clus) > 1):
            cnum += 1
            spread = fn_spread(pair_rmsds, clus)
            avg_spread += spread

    avg_spread = avg_spread / float(cnum)
    
    return avg_spread

def normalize_spread(avsp_s, n_bank):
    avsp_norm_s = []
    min_avsp = min(avsp_s)
    max_avsp = max(avsp_s)
    inv_avsp = 1.0 / float(max_avsp - min_avsp)
    n_bank2 = n_bank - 2
    for avsp in avsp_s:
        avsp_norm = inv_avsp * n_bank2 * (avsp - min_avsp) + 1.0
        avsp_norm_s.append(avsp_norm)
    return avsp_norm_s

def fn_penalty(avsp_norm_s, n_step):
    penalty_s = []
    n_clus = n_step
    for i_step in range(n_step):
        penalty = avsp_norm_s[i_step] + n_clus
        penalty_s.append(penalty)
        n_clus += (-1)
    return penalty_s

def check_cluster(cluster):
    min_len = 999
    max_len = -999
    for clus in cluster:
        len_clus = len(clus)
        if len_clus > max_len:
            max_len = len_clus
        if len_clus < min_len:
            min_len = len_clus
    
    return min_len, max_len

def reform_cluster(cluster, pair_rmsds, merge, split):
    if merge:
        i_clusmin = 0  #Search smallest cluster
        n_clus = len(cluster)
        len_clusmin = len(cluster[0])
        for i_clus in range(1, n_clus):
            len_clus = len(cluster[i_clus])
            if len_clus < len_clusmin:
                len_clusmin = len_clus
                i_clusmin = i_clus
        
        min_val = 999999.9 # Merge cluster
        d_mat = update_d_mat(cluster, pair_rmsds)
        for j_clus in range(1, n_clus):
            if i_clusmin == j_clus:
                continue
            dist = d_mat[i_clusmin][j_clus]
            if dist < min_val:
                min_val = dist
                j_clusmin = j_clus
        cluster = merge_cluster(cluster, i_clusmin, j_clusmin)
        new_cluster = cluster

    if split:
        i_clusmax = 0  #Search largest cluster
        n_clus = len(cluster)
        len_clusmax = len(cluster[0])
        for i_clus in range(1, n_clus):
            len_clus = len(cluster[i_clus])
            if len_clus > len_clusmax:
                len_clusmax = len_clus
                i_clusmax = i_clus
        input_cluster = []
        for member in cluster[i_clusmax]:
            input_cluster.append([member])
        split_s = nmrclust(input_cluster, pair_rmsds, 2) #Split cluster
        new_cluster = []
        for i_clus, clus in enumerate(cluster):
            if i_clus == i_clusmax:
                continue
            new_cluster.append(clus)
        for clus in split_s:
            new_cluster.append(clus)

    return new_cluster
        
def nmrclust(cluster, pair_rmsds, n_clus):
    n_bank = len(cluster)
    n_step = n_bank - 1
    d_mat = update_d_mat(cluster, pair_rmsds)
    
    avsp_s = []
    cluster_list = []
    for i_step in range(n_step):
        i_clus, j_clus = find_min_dist(cluster, d_mat)
        cluster = merge_cluster(cluster, i_clus, j_clus)
        d_mat = update_d_mat(cluster, pair_rmsds)
        avsp = avg_spread(pair_rmsds, cluster)
        avsp_s.append(avsp)
        cluster_list.append(cluster)
    
    if n_clus > 0:
        return cluster_list[n_step-n_clus]
        
    avsp_norm_s = normalize_spread(avsp_s, n_bank)
    penalty_s = fn_penalty(avsp_norm_s, n_step)

    for i_step in range(0, n_step):
        for clus in cluster_list[i_step]:
            if len(clus) < min_len_clus:
                penalty_s[i_step] = 99999999.9
                break

    n_clus = min_nclus
    min_val = penalty_s[(n_step-min_nclus)]

    for i_clus in range(min_nclus,max_nclus+1):
        if penalty_s[(n_step-i_clus)] < min_val:
            min_val = penalty_s[(n_step-i_clus)]
            n_clus = i_clus

    return cluster_list[(n_step-n_clus)]

def clustering(pair_rmsds):
    n_bank = len(pair_rmsds)
    cluster = []
    for i_bank in range(n_bank):
        cluster.append([i_bank])
    cluster = nmrclust(cluster, pair_rmsds, 0)
    cluster.sort(len_cmp)

    for i in range(10):
        min_len, max_len = check_cluster(cluster)
        merge = False
        split = False
        if min_len < min_len_clus:
            merge = True
            split = True
        if max_len > max_len_clus:
            split = True
        
        if i >= 5:
            split = False

        if (not merge) and (not split):
            break
        else:
            cluster = reform_cluster(cluster, pair_rmsds, merge, split)
            cluster.sort(len_cmp)
    
    while(1):
        if len(cluster) <= 10:
            break
        split = False
        merge = True
        cluster = reform_cluster(cluster, pair_rmsds, merge, split)
        cluster.sort(len_cmp)

    return cluster

def sort_cluster(cluster_s, data_s):
    new_cluster_s = []
    for clus in cluster_s:
        score_s = []
        for index in clus:
            score = data_s[index].score
            score_s.append([score, index])
        score_s.sort()
        new_list = []
        for list in score_s:
            new_list.append(list[1])
        new_cluster_s.append(new_list)
    new_cluster_s.sort(len_cmp)

    return new_cluster_s

def write_cluster(cluster_s, data_s, out_file):
    curr = os.getcwd()
    w_file = file(out_file, 'w')
    for i_clus, clus in enumerate(cluster_s):
        len_clus = len(clus)
        w_file.write('#Cluster: %4d   Length: %4d\n'%(i_clus+1, len_clus))
        for index in clus:
            pdb_file = '%s'%data_s[index].pdb_name
            w_file.write('%s\n'%pdb_file)

    w_file.close()

def parse_org_cluster(clust_file):
    clust_list = []
    for line in file(clust_file):
        if not line.startswith('#'):
            continue
        linesp = line.strip().split()
        clust_list.append(int(linesp[3]))
    
    min_len_clus = min(clust_list)
    max_len_clus = max(clust_list)
    min_nclus = min(len(clust_list), 10)
    max_nclus = min(len(clust_list), 10)
    return min_len_clus, max_len_clus, min_nclus, max_nclus

def make_cluster(org_cluster, curr_list, new_clus_file, n_res_rec):
    global min_nclus
    global max_nclus
    min_len_clus, max_len_clus, min_nclus, max_nclus = parse_org_cluster(org_cluster)

    data_s = parse_list(curr_list)
    pair_rmsds = calc_pair_rmsd(curr_list, n_res_rec)
    cluster_s = clustering(pair_rmsds)
    write_cluster(cluster_s, data_s, new_clus_file)

