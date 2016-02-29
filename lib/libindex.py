#!/usr/bin/env python
import sys, os
def file2list(pdb_file):
    list = []
    for line in file(pdb_file):
        list.append(line)
    return list

def reindexing(pdb_list_s):
    first_line = True
    index = 1
    chain_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    std_res_s = {'GLY':'G', 'ALA':'A', 'SER':'S', 'THR':'T', 'CYS':'C', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', \
            'PRO':'P', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'ASP':'D', 'GLU':'E', 'ASN':'N', 'GLN':'Q', 'LYS':'K', 'ARG':'R', \
            'HIS':'H', 'PTR':'Y', 'SEP':'S', 'TPO':'T', 'ACE':'A', 'NME':'N'}
    not_std_s = {'HSE':'HIS', 'HID':'HIS', 'HIE':'HIS'}
    c_num = 0

    i_atm = 0

    out_list = []
    for pdb_list in pdb_list_s:
        std_chain = '!'
        for line in pdb_list:
            if (not line.startswith('ATOM')) and (not line.startswith('TER')) and (not line.startswith('END')):
                out_list.append('%s\n'%line[:-1])
                continue

            res_name = line[17:20]
            if (not res_name in std_res_s) and (not res_name in not_std_s):
                continue
            if (res_name in not_std_s):
                res_name = not_std_s[res_name]

            i_atm += 1

            if first_line:
                std_chain = line[21]
                std_num = line[22:27]
                index = 1
                first_line = False
            
            chain = line[21]
            res_num = line[22:27]
            
            if line.startswith('TER'):
                std_chain = '!'
                continue

            if not chain == std_chain:
                out_list.append('TER\n')
                index += 1
                std_chain = chain
                std_num = res_num
                c_num += 1

            if not res_num == std_num:
                index += 1
                std_num = res_num

            out_list.append('%s%5d%s %3s %s%4d %s\n'%(line[0:6],i_atm,line[11:16], res_name, chain_list[c_num], index, line[27:54]))
    out_list.append('TER\n')
    return out_list

def reindexing_file(pdb_file_s, out_file):
    pdb_list_s = []
    for pdb_file in pdb_file_s:
        pdb_list = file2list(pdb_file)
        pdb_list_s.append(pdb_list)
    
    out_list = reindexing(pdb_list_s)
    w_file = file(out_file, 'w')
    for line in out_list:
        w_file.write(line)
    w_file.write('END\n')
    w_file.close()
