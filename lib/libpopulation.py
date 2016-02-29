#!/usr/bin/env python
import os, sys
no_chi_list = ['GLY', 'PRO', 'ALA']
prior_dict = {'ARG':3, 'ASN':9, 'ASP':11, 'CYS':16, 'GLN':4, 'GLU':2, 'HIS':13, 'ILE':12, 'LEU':8, 'LYS':5, 'MET':7, 'PHE':15, 'SER':1, 'THR':6, 'TRP':17, 'TYR':14, 'VAL':10}
# This priority was determined by Chi1 angle transitions of each amino acids.
# Guharoy,M. et al., Side-chain rotamer transitions at protein-protein interfaces, Proteins, 2010, 3219-3225. (Table1).

def file2lib(pdb_file):
    PDB = {}
    num2name = {}
    for line in (pdb_file):
        i_res = int(line[22:26])
        if not i_res in PDB:
            PDB[i_res] = {}
        atom_name = line[12:16]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        PDB[i_res][atom_name] = line

        res_name = line[17:20]
        num2name[i_res] = res_name

    return PDB, num2name

def dist_func2(x1,y1,z1, x2,y2,z2):
    return ( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )

def check_inter(res_rec, res_lig, len_cut2):
    is_inter = False
    for atm1 in res_rec:
        if atm1[1:2] == 'H':
            continue
        pos1 = res_rec[atm1]
        for atm2 in res_lig:
            if atm2[1:2] == 'H':
                continue
            pos2 = res_lig[atm2]

            dis2 = dis_func2(pos1, pos2)
            if dis2 < len_cut2:
                is_inter = True
                break
        if is_inter:
            break
    return is_inter

def fn_search_new(PDB_rec, PDB_lig, len_cut):
    len_cut2 = len_cut**2
    inter_s = []
    rec_inter = []
    lig_inter = []
    for i_res in PDB_rec:
        res_name1 = PDB_rec[i_res][' CA '][17:20]
        if res_name1 == 'GLY':
            pos1 = PDB_rec[i_res][' CA ']
        else:
            pos1 = PDB_rec[i_res][' CB ']
        x1 = float(pos1[30:38])
        y1 = float(pos1[38:46])
        z1 = float(pos1[46:54])
        
        for j_res in PDB_lig:
            res_name2 = PDB_lig[j_res][' CA '][17:20]
            if res_name2 == 'GLY':
                pos2 = PDB_lig[j_res][' CA ']
            else:
                pos2 = PDB_lig[j_res][' CB ']
            x2 = float(pos2[30:38])
            y2 = float(pos2[38:46])
            z2 = float(pos2[46:54])

            dis2 = dist_func2(x1,y1,z1, x2,y2,z2)
            if dis2 > len_cut2:
                continue
            inter_s.append([i_res,j_res])
            if (not i_res in rec_inter) and (not res_name1 in no_chi_list):
                rec_inter.append(i_res)
            if (not j_res in lig_inter) and (not res_name2 in no_chi_list):
                lig_inter.append(j_res)

    rec_inter.sort()
    lig_inter.sort()

    return inter_s, rec_inter, lig_inter

def calc_interface(rec_file, lig_file, inter_cut):
    PDB_rec, num2name_rec = file2lib(rec_file)
    PDB_lig, num2name_lig = file2lib(lig_file)
    inter_s, rec_inter, lig_inter = fn_search_new(PDB_rec, PDB_lig, inter_cut)
    return inter_s, rec_inter, lig_inter, num2name_rec, num2name_lig

def population(list, num2name):
    pop_lib = {}
    for i_res in list:
        if not i_res in pop_lib:
            pop_lib[i_res] = 0
        pop_lib[i_res] += 1
    
    pop_list = []
    for i_res in pop_lib:
        res_name = num2name[i_res]
        prior = prior_dict[res_name]
        pop_list.append([-pop_lib[i_res], prior, i_res])
    pop_list.sort()
    return pop_list

