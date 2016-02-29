#!/usr/bin/python
import os, sys, math, glob
class PROTEIN:
    def __init__(self, pdb_file, n_res_rec):
        self.pdb_file = pdb_file
        self.n_res_rec = n_res_rec
        self.file2dict()
    def file2dict(self):
        PDB_rec = {}
        PDB_lig = {}
        for line in file(self.pdb_file):
            if not line.startswith('ATOM'):
                continue
            if line[13] == 'H':
                continue
            atom_name = line[12:16]
            i_res = int(line[22:26])
            if i_res <= self.n_res_rec:
                if not i_res in PDB_rec:
                    PDB_rec[i_res] = {}
                PDB_rec[i_res][atom_name] = line[:-1]
            else:
                if not i_res in PDB_lig:
                    PDB_lig[i_res] = {}
                PDB_lig[i_res][atom_name] = line[:-1]
        self.PDB_rec = PDB_rec
        self.PDB_lig = PDB_lig

def dist_func2(x1,y1,z1, x2,y2,z2):
    return ( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )

def fn_CA_rmsd(PDB1, PDB2):
    res_s = PDB1.keys()
    res_s.sort()
    dist_sum2 = 0.0
    res_num = 0
    for i_res in res_s:
        x1 = float(PDB1[i_res][' CA '][30:38])
        y1 = float(PDB1[i_res][' CA '][38:46])
        z1 = float(PDB1[i_res][' CA '][46:54])
        x2 = float(PDB2[i_res][' CA '][30:38])
        y2 = float(PDB2[i_res][' CA '][38:46])
        z2 = float(PDB2[i_res][' CA '][46:54])
        dist2 = dist_func2(x1,y1,z1, x2,y2,z2)
        dist_sum2 += dist2
        res_num += 1

    dist_sum2 = dist_sum2 / float(res_num)
    rmsd = math.sqrt(dist_sum2)
    return rmsd

def fn_search_CA(PDB_rec, PDB_lig, len_cut):
    len_cut2 = len_cut**2
    inter_CA = []
    for i_res in PDB_rec:
        x1 = float(PDB_rec[i_res][' CA '][30:38])
        y1 = float(PDB_rec[i_res][' CA '][38:46])
        z1 = float(PDB_rec[i_res][' CA '][46:54])
        for j_res in PDB_lig:
            x2 = float(PDB_lig[j_res][' CA '][30:38])
            y2 = float(PDB_lig[j_res][' CA '][38:46])
            z2 = float(PDB_lig[j_res][' CA '][46:54])
            dis2 = dist_func2(x1,y1,z1, x2,y2,z2)
            if dis2 > len_cut2:
                continue
            inter_CA.append([i_res,j_res])
    
    return inter_CA

def check_inter(res_rec, res_lig, len_cut2):
    is_inter = False
    for atm1 in res_rec:
        x1 = float(res_rec[atm1][30:38])
        y1 = float(res_rec[atm1][38:46])
        z1 = float(res_rec[atm1][46:54])
        for atm2 in res_lig:
            x2 = float(res_lig[atm2][30:38])
            y2 = float(res_lig[atm2][38:46])
            z2 = float(res_lig[atm2][46:54])
            dis2 = dist_func2(x1,y1,z1, x2,y2,z2)
            if dis2 < len_cut2:
                is_inter = True
                break
        if is_inter:
            break
    return is_inter


def fn_search_all(PDB_rec, PDB_lig, len_cut, inter_CA):
    len_cut2 = len_cut**2
    inter_s = []
    for line in inter_CA:
        i_res = line[0]
        j_res = line[1]
        is_inter = check_inter(PDB_rec[i_res], PDB_lig[j_res], len_cut2)
        if is_inter:
            inter_s.append([i_res, j_res])
        
    return inter_s

def fn_fnonnat(real_inters, pred_inters):
    noreal_num = 0
    pred_num = 0
    for inter in pred_inters:
        pred_num += 1
        if not inter in real_inters:
            noreal_num += 1
    
    if pred_num <= 0:
        fnonnat = 1.0
    else:
        fnonnat = noreal_num / float(pred_num)
    return fnonnat

def fn_fnat(real_inters, pred_inters):
    real_num = 0
    pred_num = 0
    for inter in real_inters:
        real_num += 1
        if inter in pred_inters:
            pred_num += 1
    
    if real_num <= 0:
        fnat = 0.0
    else:
        fnat = pred_num / float(real_num)
    return fnat

def calc_fnat(rec_ref_dict, lig_ref_dict, rec_pred_dict, lig_pred_dict, inter_cut=5.0):
    rough_cut = inter_cut + 15.0
    inter_CA_ref = fn_search_CA(rec_ref_dict, lig_ref_dict, rough_cut)
    inter_CA_pred = fn_search_CA(rec_pred_dict, lig_pred_dict, rough_cut)
    inter_s_ref = fn_search_all(rec_ref_dict, lig_ref_dict, inter_cut, inter_CA_ref)
    inter_s_pred = fn_search_all(rec_pred_dict, lig_pred_dict, inter_cut, inter_CA_pred)
    fnat = fn_fnat(inter_s_ref, inter_s_pred)

    return fnat

def make_proteinlist(pdb_list_file, n_res_rec, top_line=999999999):
    protein_list = []
    i_pdb = 0
    for line in file(pdb_list_file):
        if line.startswith('#'):
            continue
        
        i_pdb += 1
        if i_pdb > top_line:
            break

        linesp = line.strip().split()
        pdb_file = linesp[0]
        protein = PROTEIN(pdb_file, n_res_rec)
        protein.en = float(linesp[1])
        protein_list.append(protein)
    
    return protein_list

def calc_pair_lrmsd(protein_list):
    len_p = len(protein_list)
    rmsd_list = []
    #initialize_rmsd_list
    for i in range(len_p):
        sub_list = [0.0 for j in range(len_p)]
        rmsd_list.append(sub_list)

    #calc_average_rmsd_of_rmsd_list
    for i in range(len_p-1):
        for j in range((i+1),len_p):
            PDB1 = protein_list[i].PDB_lig
            PDB2 = protein_list[j].PDB_lig
            rmsd = fn_CA_rmsd(PDB1, PDB2)
            rmsd_list[i][j] = rmsd
            rmsd_list[j][i] = rmsd

    return rmsd_list

def calc_pair_fnat(protein_list, pair_rmsd):
    len_p = len(protein_list)
    fnat_list = []
    #initialize_fnat_list
    for i in range(len_p):
        sub_list = [0.0 for j in range(len_p)]
        fnat_list.append(sub_list)

    #calc_average_fnat_of_fnat_list
    for i in range(len_p-1):
        for j in range((i+1),len_p):
            if pair_rmsd[i][j] <= 5.0:
                fnat_list[i][j] = 1.0
                fnat_list[j][i] = 1.0
            elif pair_rmsd[i][j] > 10.0:
                fnat_list[i][j] = 0.0
                fnat_list[j][i] = 0.0
            else:
                PDB1_rec = protein_list[i].PDB_rec
                PDB1_lig = protein_list[i].PDB_lig
                PDB2_rec = protein_list[j].PDB_rec
                PDB2_lig = protein_list[j].PDB_lig
                fnat1 = calc_fnat(PDB1_rec, PDB1_lig, PDB2_rec, PDB2_lig)
                fnat2 = calc_fnat(PDB2_rec, PDB2_lig, PDB1_rec, PDB1_lig)
                fnat_list[i][j] = (fnat1 + fnat2) * 0.5
                fnat_list[j][i] = fnat_list[i][j]
    return fnat_list

def write_out(pair_list, out_file):
    len_p = len(pair_list) 
    w_file = file(out_file, 'w')
    w_file.write('%-2s '%'#')
    for i in range(len_p):
        w_file.write('%5d '%(i+1))
    w_file.write('\n')
    for i in range(len_p):
        w_file.write('%2d '%(i+1))
        for j in range(len_p):
            w_file.write('%5.2f '%(pair_list[i][j]))
        w_file.write('\n')

    avg_val = 0.0
    tot_num = len_p * (len_p-1) / 2.0
    for i in range(len_p-1):
        for j in range((i+1),len_p):
            avg_val += pair_list[i][j]
    avg_val = avg_val / tot_num
    w_file.write('%s %8.3f\n'%('#AVG:', avg_val))
    w_file.close()

def make_pairlist(pdb_list_file, n_res_rec):
    protein_list = make_proteinlist(pdb_list_file, n_res_rec)
    pair_rmsd = calc_pair_lrmsd(protein_list) # calc_pairwise RMSD
    pair_fnat = calc_pair_fnat(protein_list, pair_rmsd) # calc_pairwise RMSD
    return pair_rmsd, pair_fnat

def main():
    pdb_list_file = sys.argv[1]
    n_res_rec = int(sys.argv[2])
    out_file1 = 'out1.dat'
    out_file2 = 'out2.dat'

    pair_rmsd, pair_fnat = make_pairlist(pdb_list_file, n_res_rec)
    write_out(pair_rmsd, out_file1)
    write_out(pair_fnat, out_file2) 

if __name__ == '__main__':
    main()

