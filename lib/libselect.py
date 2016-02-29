#!/usr/bin/python
import os, sys, math
sys.path.insert(0,"/home/emiliar/.python_program")
def error_msg():
    print '[*.py] [pdblist_file] [pairrmsd_file] [cut_rmsd] or ( [-all] )'
    return

def len_cmp(a,b):
    return cmp(len(b), len(a))

def parse_rmsdfile(rmsd_file):
    file_s = []
    line_s = []
    for line in file(rmsd_file):
        if line.startswith('#') or line.startswith('!'):
            print line.rstrip()
            continue
        linesp = line.strip().split()
        file_name = linesp[0]
        
        file_s.append(file_name)
        line_s.append(line[:-1])

    return file_s, line_s

def parse_pairfile(pairrmsd_file):
    pair_rmsds = []
    for line in file(pairrmsd_file):
        if line.startswith('#'):
            continue
        list = []
        linesp = line.strip().split()
        for i in range(1, len(linesp)):
            list.append(float(linesp[i]))
        pair_rmsds.append(list)
    return pair_rmsds

def get_cluster(pair_rmsds, use_dict, cut_rmsd, pair_fnats, cut_fnat):
    neigh_list = []
    for i in use_dict:
        if use_dict[i]:
            continue
        
        sub_list = [i]
        for j, rmsd in enumerate(pair_rmsds[i]):
            if use_dict[j]:
                continue
            if i == j:
                continue
            fnat = pair_fnats[i][j]
            if (rmsd <= cut_rmsd) or (fnat >= cut_fnat):
                sub_list.append(j)
        neigh_list.append(sub_list)
    
    neigh_list.sort(len_cmp)
    return neigh_list[0]

def check_use(use_dict):
    all_use = True
    for i in use_dict:
        if not use_dict[i]:
            all_use = False
            break
    return all_use

def fn_clust(input_list, pair_rmsds, cut_rmsd, pair_fnats, cut_fnat):
    use_dict = {}
    for i in input_list:
        use_dict[i] = False

    cluster_s = []
    while(1):
        sub_clus = get_cluster(pair_rmsds, use_dict, cut_rmsd, pair_fnats, cut_fnat)
        sub_clus.sort()
        cluster_s.append(sub_clus)

        for i in sub_clus:
            use_dict[i] = True
        all_use = check_use(use_dict)
        if all_use:
            break
        
    return cluster_s

def clustering(rmsdlist_file, pair_rmsds, pair_fnats, cut_rmsd=5.0, cut_fnat=0.5):
    file_s, line_s = parse_rmsdfile(rmsdlist_file)
    input_list = [i for i in range(0, len(pair_rmsds))]
    cluster_s = fn_clust(input_list, pair_rmsds, cut_rmsd, pair_fnats, cut_fnat)
    return cluster_s, line_s

