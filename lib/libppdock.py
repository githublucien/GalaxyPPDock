#!/usr/bin/env python

import os
import sys
from libpopulation import population, calc_interface
from libcl2 import make_cluster
from libpair import make_pairlist
from libselect import clustering
import libindex

to_one_letter = {'GLY':'G', 'ALA':'A', 'SER':'S', 'THR':'T', 'CYS':'C', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'PRO':'P', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'ASP':'D', 'GLU':'E', 'ASN':'N', 'GLN':'Q', 'LYS':'K', 'ARG':'R', 'HIS':'H'}
std_seq = 'ACDEFGHIKLMNPQRSTVWY'

# GALAXY PPDOCK
#
PARAM_NSELECT = 50
GALAXYPPDOCK_HOME = os.environ['GALAXYPPDOCK_HOME']
PARAM_GALAXY_DATA_HOME="%s/data/"%GALAXYPPDOCK_HOME
TEMPL_PATH="%s/templates"%GALAXYPPDOCK_HOME
EXEC_GALAXY_PPDOCK_OPT = '%s/bin/ppdock_opt.mpi'%GALAXYPPDOCK_HOME
EXEC_GALAXY_PPDOCK_CSA = '%s/bin/ppdock_csa.mpi'%GALAXYPPDOCK_HOME
TEMPL_PPDOCK_OPT = '%s/ppdock_opt.in'%TEMPL_PATH
TEMPL_PPDOCK_CSA = '%s/ppdock_csa.in'%TEMPL_PATH
EXEC_MPI='mpiexec'
EXEC_MPI=os.environ['EXEC_MPI']

def pdb2seq(pdb_file):
    seq = ''
    for line in file(pdb_file):
        if not line.startswith('ATOM'):
            continue
        if not line[12:16] == ' CA ':
            continue
        aa = to_one_letter[line[17:20]]
        seq += aa
    return seq

def read_pssm(pssm_file, i_res):
    seq = ''
    list = []
    for line in file(pssm_file):
        if line.startswith('\n'):
            continue
        linest = line.strip()
        if linest == '':
            continue
        linesp = line.strip().split()
        try:
            res_num = int(linesp[0])
        except:
            continue
        
        is_read = True
        aa = linesp[1]
        if not aa in std_seq:
            continue
        if len(linesp) < 22:
            continue
        
        for i in range(2, 22):
            try:
                score = int(linesp[i])
            except:
                score = -99
                is_read = False
            
            if int(linesp[i]) in list:
                continue
        if not is_read:
            continue
        
        i_res += 1
        seq += aa
        wrt_line = '%5d'%i_res
        wrt_line += ' %s'%aa
        for i in range(2, 22):
            score = int(linesp[i])
            wrt_line += ' %3d'%score
        
        list.append(wrt_line)
    
    return seq, list

def prefix(fn):
    return '.'.join(fn.split(".")[:-1])

def parse_input(population_list, inter_num):
    res_s = []
    n_res = 0
    for j, line in enumerate(population_list):
        if j%8 == 0:
            if j != 0:
                res_s.append(dum)
            dum = []
        i_res = int(line[0])
        dum.append(i_res)
        n_res += 1
        if n_res >= inter_num:
            break
    res_s.append(dum)
    return res_s

class Opt:
    def __init__(self, home, n_res_rec, PARAM_N_PROC):
        self.home = home
        self.n_res_rec = n_res_rec
        self.cwd = os.getcwd()
        self.PARAM_N_PROC = PARAM_N_PROC

    def prep_pdb(self, rec_pdb, ligall_pdb, select_file):
        select_s = []
        for line in file(select_file):
            if line.startswith('#'):
                continue
            linesp = line.strip().split()
            select_s.append(int(linesp[1]))
        
        rec_list = libindex.file2list(rec_pdb)
        
        self.init_pdb = '%s/init.pdb'%os.getcwd()
        wrt_dict = {}
        read = False
        for line in file(ligall_pdb):
            if line.startswith('MODEL'):
                linesp = line.strip().split()
                num = int(linesp[1])
                if num in select_s:
                    read = True
                    pdb = []
                continue
            if line.startswith('ENDMDL') and read:
                wrt_dict[num] = pdb
                read = False
                continue
            if not read:
                continue
            if line.startswith('END'):
                continue
            pdb.append(line)

        i_model = 0
        self.init_pdb = 'init.pdb'
        w_file = file(self.init_pdb, 'w')
        for num in select_s:
            if not num in wrt_dict:
                continue
            i_model += 1
            lig_list = wrt_dict[num]
            comp = libindex.reindexing([rec_list, lig_list])
            w_file.write('MODEL %d\n'%i_model)
            for line in comp:
                w_file.write(line)
            w_file.write('ENDMDL\n')
        w_file.write('END\n')
        w_file.close()

    def get_population(self):
        pdb_s = []
        for line in file(self.init_pdb):
            if line.startswith('MODEL'):
                rec_pdb = []
                lig_pdb = []
                continue
            if line.startswith('ENDMDL'):
                pdb_s.append([rec_pdb, lig_pdb])
                continue
            if not line.startswith('ATOM'):
                continue
            i_res = int(line[22:26])
            if i_res <= self.n_res_rec:
                rec_pdb.append(line[0:54])
            else:
                lig_pdb.append(line[0:54])

        inter_cut = 10.0
        n_target = 0.0
        self.NinterR = 0.0
        self.NinterL = 0.0
            
        rec_inter_tot = []
        lig_inter_tot = []
        for comp in pdb_s:
            rec_pdb = comp[0]
            lig_pdb = comp[1]
            inter_s, rec_inter, lig_inter, num2name_rec, num2name_lig \
                    = calc_interface(rec_pdb, lig_pdb, inter_cut)
            
            rec_inter_tot += rec_inter
            lig_inter_tot += lig_inter
            self.NinterR += float(len(rec_inter))
            self.NinterL += float(len(lig_inter))
            n_target += 1.0
        self.NinterR = self.NinterR / n_target
        self.NinterL = self.NinterL / n_target

        r_list_rec = population(rec_inter_tot, num2name_rec)
        r_list_lig = population(lig_inter_tot, num2name_lig)
        rec_popul = []
        lig_popul = []
        for num in r_list_rec:
            rec_popul.append([ num[2], -num[0] ])
        for num in r_list_lig:
            lig_popul.append([ num[2], -num[0] ])
        
        self.pdb_s = pdb_s
        self.rec_popul = rec_popul
        self.lig_popul = lig_popul
        self.res_s_rec = parse_input(self.rec_popul, self.NinterR)
        self.res_s_lig = parse_input(self.lig_popul, self.NinterL)
    
    def make_comp_pssm(self, rec_pssm, lig_pssm, rec_pdb, lig_pdb):
        if rec_pssm == None or lig_pssm == None:
            self.comp_pssm = None
            return
        rec_seq, rec_list = read_pssm(rec_pssm, 0)
        lig_seq, lig_list = read_pssm(lig_pssm, len(rec_list))
        recpdb_seq = pdb2seq(rec_pdb)
        ligpdb_seq = pdb2seq(lig_pdb)
        if (not rec_seq == recpdb_seq) or (not lig_seq == ligpdb_seq):
            self.comp_pssm = None
            return
        
        comp_list = rec_list + lig_list
        self.comp_pssm = '%s/comp.pssm'%os.getcwd()
        w_file = file(self.comp_pssm, 'w')
        for line in comp_list:
            w_file.write('%s\n'%line)
        w_file.close()

    def prep_in(self, opt_input, rec_pssm, lig_pssm, rec_pdb, lig_pdb):
        self.make_comp_pssm(rec_pssm, lig_pssm, rec_pdb, lig_pdb)
        wrt = file(TEMPL_PPDOCK_OPT).read()
        wrt = wrt.replace("GALAXY_DATA_HOME",PARAM_GALAXY_DATA_HOME)
        wrt = wrt.replace("COMPLEX_S",self.init_pdb)
        wrt = wrt.replace("N_RES_REC",str(self.n_res_rec))
        if not self.comp_pssm == None:
            wrt = wrt.replace("PSSM_COMPLEX", self.comp_pssm)
        else:
            wrt = wrt.replace("blastfile", "!blastfile")
            wrt = wrt.replace("conserve_weight     3.0", "conserve_weight     0.0")
        wrt = wrt.split('\n')
        for list in self.res_s_rec:
            dum = 'USC               '
            for i_res in list:
                dum += ' %d'%i_res
            wrt.append(dum)
        for list in self.res_s_lig:
            dum = 'USC               '
            for i_res in list:
                dum += ' %d'%i_res
            wrt.append(dum)
        wrt.append('!')
        wrt.append('END')
        self.in_f = opt_input
        fout = file(self.in_f, 'wt')
        for line in wrt:
            if len(line) == 0:
                continue
            fout.write('%s\n'%line)
        fout.close()

    def run(self,out_f):
        self.out_f = out_f
        if os.path.exists(out_f):
            return
        #
        cmd = []
        cmd.append(EXEC_MPI)
        cmd.append('-n %d'%min(PARAM_NSELECT, self.PARAM_N_PROC))
        cmd.append(EXEC_GALAXY_PPDOCK_OPT)
        cmd.append(self.in_f)
        cmd.append(">")
        cmd.append(self.out_f)
        cmd = ' '.join(cmd)
        #
        os.system(cmd)

    def report(self, org_clust_file, new_clust_file):
        if os.path.exists(new_clust_file):
            return
        if not os.path.exists(self.out_f):
            return
        self.out_pdb = 'out.pdb'
        curr_list = []
        pwd = os.getcwd()
        i_model = 0
        read = False
        for line in file(self.out_pdb):
            if line.startswith('MODEL'):
                i_model += 1
                pdb_file = '%s/opt%03d.pdb'%(pwd, i_model)
                pdb_list = []
                read = True
                continue
            if line.startswith('ENDMDL'):
                w_file = file(pdb_file, 'w')
                for line in pdb_list:
                    w_file.write('%s'%line)
                w_file.write('END\n')
                w_file.close()
                curr_list.append(pdb_file)
                read = False
                continue
            if line.startswith('END'):
                continue
            if not read:
                continue
            pdb_list.append(line)
        
        make_cluster(org_clust_file, curr_list, new_clust_file, self.n_res_rec)

def split_pdb(pdb_file, n_res_rec):
    pdb_name = pdb_file.rstrip('.pdb')
    rec_file = '%s_rec.pdb'%pdb_name
    lig_file = '%s_lig.pdb'%pdb_name
    w_rec = file(rec_file, 'w')
    w_lig = file(lig_file, 'w')
    
    flag_lig = False
    for line in file(pdb_file):
        if line.startswith('END'):
            continue
        if line.startswith('ATOM'):
            i_res = int(line[22:26])
            if i_res > n_res_rec:
                flag_lig = True
        
        if not flag_lig:
            w_rec.write(line)
        else:
            w_lig.write(line)
    w_rec.write('END\n')
    w_lig.write('END\n')
    
        
class CSA:
    def __init__(self, home, n_res_rec, comp_pssm, PARAM_N_PROC):
        self.home = home
        self.cwd = os.getcwd()
        self.n_res_rec = n_res_rec
        self.comp_pssm = comp_pssm
        self.PARAM_N_PROC = PARAM_N_PROC

    def read_opt_in(self, opt_in):
        self.sc_lines = []
        for line in file(opt_in):
            if not line.startswith('USC'):
                continue
            line = line.rstrip()
            self.sc_lines.append(line)

    def read_clust_file(self, clust_file):
        pwd = os.getcwd()
        self.cluster_s = []
        i_pdb = 0
        is_first = True
        for line in file(clust_file):
            if line.startswith('#Cluster:'):
                if not is_first:
                    self.cluster_s.append(clus)
                clus = []
                is_first = False
                continue
            linesp = line.strip().split()
            pdb_file = linesp[0]
            i_pdb += 1
            new_pdb = '%s/r%03d.pdb'%(pwd, i_pdb)
            os.system('cp %s %s'%(pdb_file, new_pdb))
            clus.append((new_pdb))
        self.cluster_s.append(clus)
        
    def do_csa(self, csa_in, csa_log):
        wrt = file(TEMPL_PPDOCK_CSA).read()
        wrt = wrt.replace("GALAXY_DATA_HOME",PARAM_GALAXY_DATA_HOME)
        wrt = wrt.replace("N_RES_REC",str(self.n_res_rec))
        wrt = wrt.replace("N_BANK  N_BANK",'%s     %s'%(PARAM_NSELECT, PARAM_NSELECT))
        if not self.comp_pssm == None:
            wrt = wrt.replace("PSSM_COMPLEX", self.comp_pssm)
        else:
            wrt = wrt.replace("blastfile", "!blastfile")
            wrt = wrt.replace("conserve_weight          3.0", "conserve_weight          0.0")
        wrt = wrt.split('\n')
        for line in self.sc_lines:
            wrt.append(line)
        wrt.append('!')
        wrt.append('cluster_csa              yes')
        wrt.append('cluster_number           %d'%len(self.cluster_s))
        wrt.append('cluster_change_member    1')
        wrt.append('cluster_period           2')
        wrt.append('cluster_boundary         5   20')
        i_pdb = 0
        for j, clus in enumerate(self.cluster_s):
            write_prefix = True
            index_line = 0
            for index, pdb_file in enumerate(clus):
                if write_prefix:
                    clus_line = 'cluster_index           %d  '%(j+1)
                    write_prefix = False
                i_pdb += 1
                clus_line += '%2d '%i_pdb
                if index % 8 == 7 or index >= (len(clus)-1):
                    wrt.append(clus_line)
                    write_prefix = True
        wrt.append('!')
        wrt.append('END')
        fout = file(csa_in, 'wt')
        for line in wrt:
            if len(line) == 0:
                continue
            fout.write('%s\n'%line)
        fout.close()
        
        if os.path.exists(csa_log):
            return

        cmd = []
        cmd.append(EXEC_MPI)
        cmd.append('-n %d'%self.PARAM_N_PROC)
        cmd.append(EXEC_GALAXY_PPDOCK_CSA)
        cmd.append(csa_in)
        cmd.append(">")
        cmd.append(csa_log)
        cmd = ' '.join(cmd)

        os.system(cmd)

    def read_csa_log(self, csa_log_file):
        en_list = []
        for line in file(csa_log_file):
            if not line.startswith('- ENERGY.0'):
                continue
            linesp = line.strip().split()
            en = float(linesp[2])
            en_list.append(en)
        new_en_list = []
        cut_len = len(en_list) / 2
        for index, en in enumerate(en_list):
            if index < cut_len:
                continue
            new_en_list.append(en)
        return new_en_list

    def parse_result(self, csa_log_file, csa_output_file):
        curr = os.getcwd()
        en_list = self.read_csa_log(csa_log_file)
        result_list = []
        for i_pdb in range(len(en_list)):
            pdb_file = '%s/b%03d.pdb'%(curr, i_pdb+1)
            en = en_list[i_pdb]
            result_list.append([en, pdb_file])
        result_list.sort()

        w_file = file(csa_output_file, 'w')
        for sub_list in result_list:
            en = sub_list[0]
            pdb_file = sub_list[1]
            w_file.write('%s %12.3f\n'%(pdb_file, en))
        w_file.close()
        self.csa_output_file = csa_output_file
    
    def make_cluster(self, clust_file):
        pair_rmsd, pair_fnat = make_pairlist(self.csa_output_file, self.n_res_rec)
        self.cluster_s, self.line_s = clustering(self.csa_output_file, pair_rmsd, pair_fnat)
        w_file = file(clust_file, 'w')
        for i_clus, sub_clus in enumerate(self.cluster_s):
            w_file.write('#Cluster: %4d    Length: %4d\n'%(i_clus+1, len(sub_clus)))
            for i in sub_clus:
                w_file.write('%s\n'%self.line_s[i])
        w_file.close()

    def select10(self, select10_file):
        en_min = 999999999.9
        line_s = []
        for sub_clus in self.cluster_s:
            i = sub_clus[0]
            line_s.append(self.line_s[i])

        list = []
        for line in line_s:
            linesp = line.strip().split()
            pdb_file = linesp[0]
            en = float(linesp[1])
            list.append([pdb_file, en])
            if en < en_min:
                en_min = en
                pdb_min = pdb_file
        
        w_file = file(select10_file, 'w')
        n_pdb = 1
        w_file.write('%-80s %12.3f\n'%(pdb_min, en_min))
        for sub_list in list:
            pdb_file = sub_list[0]
            en = sub_list[1]
            if pdb_file == pdb_min:
                continue
            w_file.write('%-80s %12.3f\n'%(pdb_file, en))
            n_pdb += 1
            if n_pdb >= 10:
                break
        w_file.close()
    
    def finalize(self, select10_file):
        i_model = 0
        for line in file(select10_file):
            if line.startswith('#'):
                continue
            linesp = line.strip().split()
            pdb_file = linesp[0]
            i_model += 1
            new_file = 'model%02d.pdb'%i_model
            if os.path.exists(new_file):
                os.system('rm %s'%new_file)
            os.system('ln -s %s %s'%(pdb_file, new_file))
            split_pdb(new_file, self.n_res_rec)
            
