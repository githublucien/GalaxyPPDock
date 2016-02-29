#!/usr/bin/python
import sys, os
resref_dict = {'GLY':'GLY', 'ALA':'ALA', 'SER':'SER', 'THR':'THR', 'CYS':'CYS', 'VAL':'VAL', 'LEU':'LEU', 'ILE':'ILE', 'MET':'MET', 'PRO':'PRO', 'PHE':'PHE', 'TYR':'TYR', 'TRP':'TRP', 'ASP':'ASP', 'GLU':'GLU', 'ASN':'ASN', 'GLN':'GLN', 'LYS':'LYS', 'ARG':'ARG', 'HIS':'HIS', 'HID':'HIS', 'HIE':'HIS', 'MSE':'MET'}
GALAXYPPDOCK_HOME = os.environ['GALAXYPPDOCK_HOME']
EXEC_ADD_MISS = '%s/bin/add_miss'%GALAXYPPDOCK_HOME
    
def make_aaref_dict():
    aaref_dict = {}
    bb_s = [' C  ',' O  ']
    aaref_dict['GLY'] = [' N  ',' CA ']
    aaref_dict['ALA'] = [' N  ',' CA ',' CB ']
    aaref_dict['SER'] = [' N  ',' CA ',' CB ',' OG ']
    aaref_dict['THR'] = [' N  ',' CA ',' CB ',' OG1',' CG2']
    aaref_dict['CYS'] = [' N  ',' CA ',' CB ',' SG ']
    aaref_dict['VAL'] = [' N  ',' CA ',' CB ',' CG1',' CG2']
    aaref_dict['LEU'] = [' N  ',' CA ',' CB ',' CG ',' CD1',' CD2']
    aaref_dict['ILE'] = [' N  ',' CA ',' CB ',' CG1',' CG2',' CD1']
    aaref_dict['MET'] = [' N  ',' CA ',' CB ',' CG ',' SD ',' CE ']
    aaref_dict['PRO'] = [' N  ',' CA ',' CB ',' CG ',' CD ']
    aaref_dict['PHE'] = [' N  ',' CA ',' CB ',' CG ',' CD1',' CD2',' CE1',' CE2',' CZ ']
    aaref_dict['TYR'] = [' N  ',' CA ',' CB ',' CG ',' CD1',' CD2',' CE1',' CE2',' CZ ',' OH ']
    aaref_dict['TRP'] = [' N  ',' CA ',' CB ',' CG ',' CD1',' CD2',' NE1',' CE2',' CE3',' CZ2',' CZ3',' CH2']
    aaref_dict['ASP'] = [' N  ',' CA ',' CB ',' CG ',' OD1',' OD2']
    aaref_dict['GLU'] = [' N  ',' CA ',' CB ',' CG ',' CD ',' OE1',' OE2']
    aaref_dict['ASN'] = [' N  ',' CA ',' CB ',' CG ',' OD1',' ND2']
    aaref_dict['GLN'] = [' N  ',' CA ',' CB ',' CG ',' CD ',' OE1',' NE2']
    aaref_dict['HIS'] = [' N  ',' CA ',' CB ',' CG ',' ND1',' CD2',' CE1',' NE2']
    aaref_dict['LYS'] = [' N  ',' CA ',' CB ',' CG ',' CD ',' CE ',' NZ ']
    aaref_dict['ARG'] = [' N  ',' CA ',' CB ',' CG ',' CD ',' NE ',' CZ ',' NH1',' NH2']
    for res_name in aaref_dict:
        aaref_dict[res_name] += bb_s
    return aaref_dict

def file2dict(pdb_file):
    resnum_list = []
    PDB = {}
    for line in file(pdb_file):
        if (not line.startswith('ATOM')) and (not line.startswith('HETATM')):
            continue
        if not line[17:20] in resref_dict:
            continue
        res_name = resref_dict[line[17:20]]
        res_num = line[21:27]
        if not res_num in resnum_list:
            resnum_list.append(res_num)
            PDB[res_num] = {}
        atom_name = line[12:16]
        if atom_name in PDB[res_num]:
            continue
        PDB[res_num][atom_name] = line[0:54]
    return resnum_list, PDB 

def write_pdb(PDB, resnum_list, out_file):
    wrt_list = []
    aaref_dict = make_aaref_dict()
    i_res = 0
    i_atm = 0
    for res_num in resnum_list:
        if not ' N  ' in PDB[res_num]:
            continue
        if not ' CA ' in PDB[res_num]:
            continue
        if not ' C  ' in PDB[res_num]:
            continue
        if not ' O  ' in PDB[res_num]:
            continue
        res_name = resref_dict[PDB[res_num][' N  '][17:20]]
        i_res += 1
        for atom_name in aaref_dict[res_name]:
            if not atom_name in PDB[res_num]:
                continue
            i_atm += 1
            line = PDB[res_num][atom_name]
            wrt_line = '%-6s%5d %4s %3s %s%4d %s\n'%\
                    ('ATOM', i_atm, atom_name, res_name, line[21], i_res, line[27:54])
            wrt_list.append(wrt_line)
    
    w_file = file(out_file, 'w')
    for line in wrt_list:
        w_file.write(line)
    w_file.write('TER\n')
    w_file.close()

def write_infile(infile, infile_pdb):
    w_file = file(infile, 'w')
    w_file.write('data_directory %s/data/\n'%GALAXYPPDOCK_HOME)
    w_file.write('infile_pdb     %s\n'%infile_pdb)
    w_file.write('top_type       polarh\n')
    w_file.write('END\n')
    w_file.close()
    

def modify_pdb(infile_pdb, outfile_pdb):
    resnum_list, PDB = file2dict(infile_pdb)
    write_pdb(PDB, resnum_list, outfile_pdb)
    write_infile('add_miss.in', outfile_pdb)
    os.system('%s add_miss.in'%EXEC_ADD_MISS)
    os.system('mv %s.add %s'%(outfile_pdb, outfile_pdb))
    
