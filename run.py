#!/usr/bin/python

import os
import sys
sys.path.insert(0, '%s/lib'%os.environ['GALAXYPPDOCK_HOME'])
import libreadpdb
import libzdock
import libppdock

def report_msg():
    print "USAGE"
    print "./run.py [RECEPTOR_PDB] [LIGAND_PDB] [ZDOCK_OUTPUT_FILE] [-pssm] [RECEPTOR_PSSM_FILE] [LIGAND_PSSM_FILE] [-n] [2]"
    print
    print "EXAMPLE" 
    print "Go to examples/TA12/ directory"
    print "../../run.py receptor.pdb ligand.pdb zdock.out -pssm receptor.pssm ligand.pssm -n 2"

def main():
    pwd = os.getcwd()
    PARAM_N_PROC = 2

    #input files
    if len(sys.argv) < 4:
        report_msg()
        return
    rec_pdb_org = sys.argv[1]
    lig_pdb_org = sys.argv[2]
    zdock_out = sys.argv[3]
    if not os.path.exists(rec_pdb_org):
        report_msg()
        return
    if not os.path.exists(lig_pdb_org):
        report_msg()
        return
    if not os.path.exists(zdock_out):
        report_msg()
        return
    rec_pdb_org = os.path.abspath(rec_pdb_org)
    lig_pdb_org = os.path.abspath(lig_pdb_org)
    zdock_out = os.path.abspath(zdock_out)
    
    #read pssm
    rec_pssm = None
    lig_pssm = None
    for i_sys in range(len(sys.argv)):
        if sys.argv[i_sys] == '-pssm':
            if len(sys.argv) >= i_sys + 3:
                rec_pssm = sys.argv[i_sys+1]
                lig_pssm = sys.argv[i_sys+2]
                if not os.path.exists(rec_pssm):
                    rec_pssm = None
                if not os.path.exists(lig_pssm):
                    lig_pssm = None

    if not rec_pssm == None:
        rec_pssm = os.path.abspath(rec_pssm)
    if not lig_pssm == None:
        lig_pssm = os.path.abspath(lig_pssm)
    if rec_pssm == None or lig_pssm == None:
        rec_pssm = None
        lig_pssm = None
    
    #set number of cpu
    for i_sys in range(len(sys.argv)):
        if sys.argv[i_sys] == '-n' or sys.argv[i_sys] == '-np':
            if len(sys.argv) >= i_sys + 2:
                try:
                    n_proc = int(sys.argv[i_sys+1])
                except:
                    n_proc = 2
            PARAM_N_PROC = max(PARAM_N_PROC, n_proc)
                
    #initialize input_pdb
    rec_pdb = '%s/rec.pdb'%pwd
    lig_pdb = '%s/lig.pdb'%pwd
    libreadpdb.modify_pdb(rec_pdb_org, rec_pdb)
    libreadpdb.modify_pdb(lig_pdb_org, lig_pdb)
    
    #output files
    zdock_home = '%s/zdock'%pwd
    complex_s = '%s/complex_s'%zdock_home
    energy_in = '%s/energy.in'%zdock_home
    energy_log = '%s/energy.log'%zdock_home
    zdock_enr = '%s/zdock.enr'%zdock_home
    cluster_fn = '%s/zdock.cl'%zdock_home
    selected_fn = '%s/zdock.select'%zdock_home
    
    init_home = '%s/init'%pwd
    opt_input = '%s/opt.in'%init_home
    opt_output = '%s/opt.log'%init_home
    clust_output = '%s/opt.newclus'%init_home
    
    csa_home = '%s/csa'%pwd
    csa_input = '%s/csa.in'%csa_home
    csa_input_clus = '%s/csa.in.clus'%csa_home
    csa_log = '%s/csa.log'%csa_home
    csa_output = '%s/csa.output'%csa_home
    csa_clust = '%s/csa.clust'%csa_home
    csa_select10 = '%s/csa.select10'%csa_home

    if not os.path.exists(zdock_home): #ZDOCK START
        os.system('mkdir %s/'%zdock_home)
    os.chdir(zdock_home)

    zdock = libzdock.Zdock(zdock_home, rec_pdb, lig_pdb, zdock_out)
    zdock.make()

    score = libzdock.Score(zdock_home, zdock.n_res_rec)
    score.prep_in(complex_s, energy_in, zdock.ligall_pdb)
    score.run(energy_log)
    score.report(zdock_out, zdock_enr)

    score.score_read(zdock_enr)
    cluster,report = libzdock.build_cluster(zdock_home,score.score_s, zdock.ligall_pdb)
    if not os.path.exists(cluster_fn):
        fout = file(cluster_fn,'wt')
        fout.writelines(report)
        fout.close()

    selected, report = libzdock.select_cluster(zdock_home, cluster_fn)
    fout = file(selected_fn, 'wt')
    fout.writelines(report)
    fout.close()

    os.chdir(pwd)


    if not os.path.exists(init_home): #initial_optimize
        os.system('mkdir %s/'%init_home)
    os.chdir(init_home)

    opt = libppdock.Opt(init_home, zdock.n_res_rec, PARAM_N_PROC)
    opt.prep_pdb(rec_pdb, zdock.ligall_pdb, selected_fn)
    opt.get_population()
    opt.prep_in(opt_input, rec_pssm, lig_pssm, rec_pdb, lig_pdb)
    opt.run(opt_output)
    opt.report(selected_fn, clust_output)

    os.chdir(pwd)


    if not os.path.exists(csa_home): #run PPDock-CSA (Cluster-Guided)
        os.system('mkdir %s/'%csa_home)
    os.chdir(csa_home)

    csa = libppdock.CSA(csa_home, zdock.n_res_rec, opt.comp_pssm, PARAM_N_PROC)
    csa.read_opt_in(opt_input)
    csa.read_clust_file(clust_output)
    csa.do_csa(csa_input, csa_log)
    csa.parse_result(csa_log, csa_output)
    csa.make_cluster(csa_clust)
    csa.select10(csa_select10)
    
    os.chdir(pwd)
    
    os.system('ln -s %s'%csa_select10)
    csa.finalize(csa_select10)

if __name__=='__main__':
    main()
