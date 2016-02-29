#!/usr/bin/env python

import os
import sys
from libcl     import NMRclust
from libbuild import make_ligprotein
from libindex import reindexing_file

to_one_letter = {'GLY':'G', 'ALA':'A', 'SER':'S', 'THR':'T', 'CYS':'C', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'PRO':'P', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'ASP':'D', 'GLU':'E', 'ASN':'N', 'GLN':'Q', 'LYS':'K', 'ARG':'R', 'HIS':'H'}

# GALAXY PPDOCK
#
PARAM_NFILT = 200
PARAM_NSELECT = 50
PARAM_CUT_RATIO = 0.6

GALAXYPPDOCK_HOME = os.environ['GALAXYPPDOCK_HOME']
PARAM_GALAXY_DATA_HOME="%s/data"%GALAXYPPDOCK_HOME
EXEC_GALAXY_PPDOCK = '%s/bin/energy_ligprotein'%GALAXYPPDOCK_HOME

TEMPL_PATH="%s/templates"%GALAXYPPDOCK_HOME
TEMPL_ZDOCK_ENERGY = '%s/ppdock_energy.in'%TEMPL_PATH
TEMPL_PPDOCK_OPT = '%s/ppdock_opt.in'%TEMPL_PATH
TEMPL_PPDOCK_CSA = '%s/ppdock_csa.in'%TEMPL_PATH
base_score = -0.1

def calc_n_res(pdb_file):
    n_res = 0
    for line in file(pdb_file):
        if not line.startswith('ATOM'):
            continue
        if not line[12:16] == ' CA ':
            continue
        n_res += 1
    return n_res

class Zdock:
    def __init__(self, home, rec, lig, zdock_out):
        self.home = home
        self.cwd = os.getcwd()
        #
        self.pdb_s = []
        self.pdb_s.append(rec)
        self.pdb_s.append(lig)
        self.receptor_fn = os.path.abspath(rec)
        self.ligand_fn = os.path.abspath(lig)
        self.n_res_rec = calc_n_res(self.receptor_fn)
        self.ligall_pdb = '%s/lig_all.pdb'%os.getcwd()
        self.out_fn = zdock_out
    def make(self):
        if not os.path.exists(self.ligall_pdb):
            self.comp_fn = 'ref.pdb'
            reindexing_file([self.receptor_fn, self.ligand_fn], self.comp_fn)
            make_ligprotein(self.out_fn, self.ligall_pdb, self.receptor_fn, self.ligand_fn)

class Score:
    def __init__(self,home, n_res_rec):
        self.home = home
        self.cwd = os.getcwd()
        self.n_res_rec = n_res_rec

    def prep_in(self,complex_s,in_f, pdb_file):
        self.in_f = in_f
        wrt = file(TEMPL_ZDOCK_ENERGY).read()
        wrt = wrt.replace("GALAXY_DATA_HOME",PARAM_GALAXY_DATA_HOME)
        wrt = wrt.replace("COMPLEX_S",os.path.relpath(complex_s))
        wrt = wrt.replace("N_RES_REC",'%s'%self.n_res_rec)
        wrt = wrt.replace('REFERENCE_PDB', 'ref.pdb')
        wrt = wrt.replace('LIGAND_PROTEINS', pdb_file)
        wrt = wrt.split("\n")
        fout = file(in_f,'wt')
        for line in wrt:
            fout.write('%s\n'%line)
        fout.close()

    def run(self,out_f):
        self.out_f = out_f
        if os.path.exists(out_f):
            return
        #
        cmd = []
        cmd.append(EXEC_GALAXY_PPDOCK)
        cmd.append(self.in_f)
        cmd.append(">")
        cmd.append(self.out_f)
        cmd = ' '.join(cmd)
        #
        os.system(cmd)

    def report(self,zdock_out, zdock_enr):
        self.score_s = []
        #
        ith = 0
        for line in file(zdock_out):
            x = line.strip().split()
            if len(x) < 7:
                continue
            ith += 1
            self.score_s.append([ith,-float(x[6])]) # 0, 1
        #
        read = False
        ith = 0
        for line in file(self.out_f):
            if not line.startswith('- ENERGY.2.0'):
                continue
            linesp = line.strip().split()
            dfire = float(linesp[13])
            elec = float(linesp[7])
            self.score_s[ith].append(dfire)
            self.score_s[ith].append(elec)
            ith += 1
        #
        self.score_s = dat_to_Z(self.score_s)
        self.score_s.sort(key=lambda x:x[-1])
        #
        report = []
        report.append('%4s %10s %12s %12s\n'%\
            ('# ID','ZDOCK','DFIRE','ELEC'))
        #
        is_selected = False
        self.selected = -1
        fmt = '%4d '+' %10.3f'+' %12.3f'*2+ ' ' + ' %9.3f'*4
        for i in range(ith):
            if not is_selected and self.score_s[i][-2] < self.score_s[i][-1]:
                is_selected = True
                report.append('%s\n'%fmt%tuple(self.score_s[i]))
                self.selected = self.score_s[i][0]-1
            else:
                report.append('%s\n'%fmt%tuple(self.score_s[i]))
        #
        if not os.path.exists(zdock_enr):
            w_file = file(zdock_enr, 'w')
            for line in report:
                w_file.write(line)
            w_file.close()
        return report

    def score_read(self, zdock_enr):
        self.score_s = []
        for line in file(zdock_enr):
            if line.startswith('#'):
                continue
            linesp = line.strip().split()
            j = int(linesp[0])
            score = float(linesp[7])
            self.score_s.append([j, score])

def dijsq(X,Y):
    d2 = 0.0
    for i in range(3):
        d2 += (X[i]-Y[i])**2
    return d2

def stat(X):
    n = float(len(X))
    X2 = [x**2 for x in X]
    m = sum(X)/n
    m2 = sum(X2)/n
    s = (m2-m**2)**0.5
    return m,s

def dat_to_Z(dat):
    n_comp = len(dat[0])
    X = [[] for i in range(1,n_comp)]
    for x in dat:
        for i in range(1,n_comp):
            X[i-1].append(x[i])
    #
    for x in X:
        m,s = stat(x)
        for i in range(len(x)):
            dat[i].append((x[i]-m)/s)
    #
    for i in range(len(dat)):
        sum_score = dat[i][-3] + dat[i][-2] + dat[i][-1]
        dat[i].append(sum_score)
    #
    return dat
    
def read_pdb_nres(pdb_fn):
    tnr = 0
    for line in file(pdb_fn):
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != 'CA':
            continue
        tnr += 1
    return tnr

def file2list(pdb_file):
    pdb_list = []
    for line in file(pdb_file):
        if not line.startswith('ATOM'):
            pdb_list.append(line)
            continue
        if line[12:16].strip() != 'CA':
            continue
        pdb_list.append(line[0:54])
    return pdb_list

def read_pdbs(monomer_s, pdb_file):
    pdb = file2list(pdb_file)
    R_s = []
    for fn in monomer_s:
        R = []
        read = False
        for line in pdb:
            if line.startswith('MODEL'):
                linesp = line.strip().split()
                if int(linesp[1])-1 == fn[0]:
                    read = True
                    continue
            if read and line.startswith('ENDMDL'):
                read = False
                break
            if not line.startswith("ATOM"):
                continue
            if not read:
                continue
            R.append([float(line[30+8*i:38+8*i]) for i in range(3)])
        R_s.append(R)
    return R_s

def update_pdb(pdb_fn):
    pdb = []
    fp = file(pdb_fn)
    for line in fp:
        if not line.startswith("ATOM"):
            continue
        pdb.append(line)
    fp.close()
    #
    fout = file(pdb_fn,'wt')
    fout.writelines(pdb)
    fout.write("TER\nEND\n")
    fout.close()

def CArmsd(X,Y):
    dev = 0.0
    for i in range(len(X)):
        dev += dijsq(X[i],Y[i])
    dev /= float(len(X))
    return dev**0.5

def build_cluster(zdock_home,score_s, pdb_file):
    monomer_s = []
    for i in range(PARAM_NFILT):
        j = score_s[i][0]-1
        monomer_s.append((j,j+1, score_s[i][-1]))
    #
    R_s = read_pdbs(monomer_s, pdb_file)
    dmat = {}
    for i in range(PARAM_NFILT):
        dmat[(monomer_s[i],monomer_s[i])] = 0.0
    for i in range(PARAM_NFILT-1):
        for j in range(i+1,PARAM_NFILT):
            rmsd = CArmsd(R_s[i],R_s[j])
            dmat[(monomer_s[i],monomer_s[j])] = rmsd
            dmat[(monomer_s[j],monomer_s[i])] = rmsd
    #
    n_cluster,cluster = NMRclust(monomer_s,dmat)
    #
    report = []
    for icl in range(n_cluster):
        cluster[icl].sort(key=lambda x:x[2])
        report.append("#Cluster: %4d  Length: %4d\n"%(icl+1,len(cluster[icl])))
        for i,x in enumerate(cluster[icl]):
            report.append("%s/ %4d %8.3f\n"%(zdock_home,x[1],x[2]))
    return cluster,report

def renumbering(cluster_s, n_bank):
    tot_num = 0
    tot_score = 0.0
    for list in cluster_s:
        score = list[0]
        clus = list[1]
        tot_score += abs(score)
    
    nseed_s = []
    seed_num = 0
    for list in cluster_s:
        score = abs(list[0])
        clus = list[1]
        nseed = min( int(n_bank* score/tot_score), len(clus))
        seed_num += nseed
        nseed_s.append(nseed)
    
    while(1):
        for i in range(len(nseed_s)):
            if seed_num >= n_bank:
                break
            clus = cluster_s[i][1]
            if nseed_s[i] < len(clus):
                nseed_s[i] += 1
                seed_num += 1
        if seed_num >= n_bank:
            break

    return nseed_s

def pick_bank(new_clusters, nseed_s):
    clus_final = []
    for j, list in enumerate(new_clusters):
        nseed = nseed_s[j]
        clus = list[1]
        clus_new = []
        for iseed, line in enumerate(clus):
            if iseed >= nseed:
                break
            clus_new.append(line)
        clus_final.append(clus_new)
    return clus_final

def select_cluster(zdock_home, cluster_file):
    cluster = []
    first_line = True
    for line in file(cluster_file):
        if line.startswith('#'):
            if not first_line:
                cluster.append(cl)
            first_line = False
            cl = []
            continue
        linesp = line.strip().split()
        pdb_file = linesp[1]
        i_pdb = int(pdb_file)
        score = min(float(linesp[2]), base_score)
        cl.append((i_pdb, pdb_file, score))
    cluster.append(cl)
        
    filt = []
    tot_score = 0.0
    for cl in cluster:
        sum = 0.0
        for x in cl:
            sum += x[2]
        filt.append((sum,cl))
        tot_score += sum
    filt.sort()

    selected = []
    sum = 0.0
    n_conf = 0
    for cl in filt:
        selected.append(cl)
        sum += cl[0]
        n_conf += len(cl)
        if (sum/tot_score) > PARAM_CUT_RATIO:
            break
    nseed_s = renumbering(selected, PARAM_NSELECT)
    clus_final = pick_bank(selected, nseed_s)

    report1 = []
    i_pdb = 0
    for icl, clus in enumerate(clus_final):
        report1.append("#Cluster: %4d  Length: %4d\n"%(icl+1,len(clus)))
        for list in clus:
            i_pdb += 1
            pdb_file1 = '%s %4d'%(zdock_home, int(list[1]))
            report1.append("%s %10.3f\n"%(pdb_file1, list[2]))
    return clus_final, report1

def prefix(fn):
    return '.'.join(fn.split(".")[:-1])

