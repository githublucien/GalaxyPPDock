#!/usr/bin/env python
import sys, os, glob, math
pi = 3.14159265359
two_pi = 2.0*pi
rad2deg = 180.0 / pi
small_real = 0.00001

def file2list(pdb_file):
    list = []
    for line in file(pdb_file):
        if line.startswith('END'):
            continue
        line = line.rstrip()
        list.append(line)
    return list

class ZDOCK:
    def __init__(self):
        self.grid = None

def read_zdock_output(zdock_output, rec_file, lig_file):
    rec_pdb = file2list(rec_file)
    lig_pdb = file2list(lig_file)
    i_line = 0
    zdock_s = []
    for line in file(zdock_output):
        i_line += 1
        if i_line == 1:
            linesp = line.strip().split()
            n_grid = float(int(linesp[0])) / 2.0
            grid_step = float(linesp[1])
        if i_line == 3:
            linesp = line.strip().split()
            center2 = [float(linesp[1]), float(linesp[2]), float(linesp[3])]
        if i_line == 4:
            linesp = line.strip().split()
            center = [float(linesp[1]), float(linesp[2]), float(linesp[3])]
        if i_line >= 5:
            linesp = line.strip().split()
            zdock = ZDOCK()
            zdock.alpha = float(linesp[0])
            zdock.beta = float(linesp[1])
            zdock.gamma = float(linesp[2])
            zdock.nx = float(linesp[3])
            zdock.ny = float(linesp[4])
            zdock.nz = float(linesp[5])
            zdock_s.append(zdock)
    
    return rec_pdb, lig_pdb, center, center2, zdock_s, n_grid, grid_step

def euler_to_rotmat(angle):
    c1 = math.cos(angle[0])
    c2 = math.cos(angle[1])
    c3 = math.cos(angle[2])
    s1 = math.sin(angle[0])
    s2 = math.sin(angle[1])
    s3 = math.sin(angle[2])
    mat_rot = []
    mat_rot.append([c1*c3-c2*s1*s3, -c1*s3-c2*c3*s1, s1*s2])
    mat_rot.append([c3*s1+c1*c2*s3, -s1*s3+c1*c2*c3, -c1*s2])
    mat_rot.append([s2*s3, c3*s2, c2])
    return mat_rot

def rotate(pos, alpha, beta, gamma):
    mat_rot = euler_to_rotmat([alpha, beta, gamma])
    x = mat_rot[0][0] * pos[0] + mat_rot[0][1] * pos[1] + mat_rot[0][2] * pos[2]
    y = mat_rot[1][0] * pos[0] + mat_rot[1][1] * pos[1] + mat_rot[1][2] * pos[2]
    z = mat_rot[2][0] * pos[0] + mat_rot[2][1] * pos[1] + mat_rot[2][2] * pos[2]
    return [x,y,z]
    
def translate(coord, x_move, y_move, z_move):
    new_coord = [0.0, 0.0, 0.0]
    new_coord[0] = coord[0] + x_move
    new_coord[1] = coord[1] + y_move
    new_coord[2] = coord[2] + z_move
    return new_coord

def translate_pdb(pdb_list, x_move, y_move, z_move):
    pdb_list_new = []
    for line in pdb_list:
        if not line.startswith('ATOM'):
            pdb_list_new.append(line)
            continue
        coord = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        new_coord = translate(coord, x_move, y_move, z_move)
        new_line = '%s%8.3f%8.3f%8.3f'%(line[0:30], new_coord[0], new_coord[1], new_coord[2])
        pdb_list_new.append(new_line)
    return pdb_list_new

def transrot_pdb(pdb_list, alpha, beta, gamma, x_move=0.0, y_move=0.0, z_move=0.0):
    pdb_list_new = []
    for line in pdb_list:
        if not line.startswith('ATOM'):
            pdb_list_new.append(line)
            continue
        coord = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        new_coord = rotate(coord, alpha, beta, gamma)
        new_coord = translate(new_coord, x_move, y_move, z_move)
        new_line = '%s%8.3f%8.3f%8.3f'%(line[0:30], new_coord[0], new_coord[1], new_coord[2])
        pdb_list_new.append(new_line)
    return pdb_list_new

def build_pdbs(rec_pdb, lig_pdb, center, center2, zdock_s, n_grid, grid_step):
    lig_pdb2 = translate_pdb(lig_pdb, -center[0], -center[1], -center[2])
    pdb_s = []
    for zdock in zdock_s:
        if (zdock.nx > n_grid):
            x_move = grid_step*(-1.0*zdock.nx + 2.0*n_grid) + center2[0]
        else:
            x_move = grid_step*(-1.0*zdock.nx) + center2[0]
        if (zdock.ny > n_grid):
            y_move = grid_step*(-1.0*zdock.ny + 2.0*n_grid) + center2[1]
        else:
            y_move = grid_step*(-1.0*zdock.ny) + center2[1]
        if (zdock.nz > n_grid):
            z_move = grid_step*(-1.0*zdock.nz + 2.0*n_grid) + center2[2]
        else:
            z_move = grid_step*(-1.0*zdock.nz) + center2[2]
        lig_pdb = transrot_pdb(lig_pdb2, zdock.alpha, zdock.beta, zdock.gamma, x_move, y_move, z_move)
        pdb_s.append(lig_pdb)
    return pdb_s

def write_pdbs(pdb_s, out_file):
    w_file = file(out_file, 'w')
    i_pdb = 0
    for pdb in pdb_s:
        i_pdb += 1
        w_file.write('MODEL %s\n'%i_pdb)
        for line in pdb:
            w_file.write('%s\n'%line)
        w_file.write('ENDMDL\n')
    w_file.write('END\n')
    w_file.close()
    
def make_ligprotein(zdock_file, ligprotein_file, rec_file, lig_file):
    rec_pdb, lig_pdb, center, center2, zdock_s, n_grid, grid_step = read_zdock_output(zdock_file, rec_file, lig_file)
    pdb_s = build_pdbs(rec_pdb, lig_pdb, center, center2, zdock_s, n_grid, grid_step)
    write_pdbs(pdb_s, ligprotein_file)
